#! /usr/bin/env/python

from __future__ import print_function
import argparse
import csv
import datetime
import glob
import multiprocessing
import numpy as np
import os.path
import pandas as pd
import sys
import matplotlib as mpl
import pybedtools
# to prevent DISPLAY weirdness when running in the cluster:
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from functools import partial
from operator import itemgetter
#from functools import reduce


# returns ith column from the given matrix
def get_column(matrix,i):
    f = itemgetter(i)
    return map(f,matrix)


# return whether the given interval is within a window of size window_size
# around the given mean
def is_in_window(motif_interval, atac_median, window_size):
    start = int(motif_interval.start)
    end = int(motif_interval.end)
    if (end >= (atac_median - window_size) and end <= (atac_median + window_size)) \
        or (start >= (atac_median - window_size) and start <= (atac_median + window_size)):
        return True
    else:
        return False


# this will be ran in parallel
def find_motifs_in_chrom(current_chrom, files):
    tf_motif_filename, atac_peaks_filename, radius = files
    H = radius          # in bps, the MD-score parameter (large window)
    h = radius * 0.1           # in bps, the MD-score parameter (small window)
    REPRESSOR_MARGIN = 500      # in bps, distance from the large window boundaries

    atac_df = pd.read_csv(atac_peaks_filename, header=None, comment='#', sep="\t", usecols=[0, 1, 2], \
                          names=['chrom', 'start', 'end'], na_filter=False, dtype={'chrom':'str', 'start':'int', 'end':'int'})
    atac_iter = atac_df[(atac_df.chrom == current_chrom)].itertuples()
    motif_df = pd.read_csv(tf_motif_filename, header=None, comment='#', sep="\t", usecols=[0, 1, 2], \
                           names=['chrom', 'start', 'end'], na_filter=False, dtype={'chrom':'str', 'start':'int', 'end':'int'})
    if len(motif_df) == 0:
        return None
    motif_iter = motif_df[(motif_df.chrom == current_chrom)].itertuples()
    last_motif = None
    g_h = 0
    g_H = 0
    total_motif_sites = 0
    tf_distances = []
    try:
        motif_region = next(motif_iter)   # get the first motif in the list
    except StopIteration:
        print('No motifs for chromosome ' + current_chrom + ' on file ' + tf_motif_filename)
        return None

    peak_count_overlapping_motif = 0
    for atac_peak in atac_iter:
        motifs_within_region = True

        atac_median = atac_peak.start + (atac_peak.end - atac_peak.start)/2
        # check the last motif, too, in case any of them falls within the region of
        # interest of two sequential ATAC-Seq peaks
        if last_motif:
            if is_in_window(last_motif, atac_median, h):
                g_h = g_h + 1
            if is_in_window(last_motif, atac_median, H):
                g_H = g_H + 1
                tf_median = last_motif.start + \
                            (last_motif.end - last_motif.start)/2
                tf_distances.append(atac_median - tf_median)
                try:
                    motif_region = next(motif_iter)   # get the next motif in the list
                    total_motif_sites += 1
                except StopIteration:
                    pass
                last_motif = motif_region
                if motif_region.start > (atac_median + H):
                    motifs_within_region = False

        # Move to the next putative motif sites until we get one past our evaluation window
        while motifs_within_region:
            total_motif_sites += 1
            # account for those within the smaller window (h)
            if is_in_window(motif_region, atac_median, h):
                g_h = g_h + 1
            # account for those within the larger window (H)
            if is_in_window(motif_region, atac_median, H):
                g_H = g_H + 1
                tf_median = motif_region.start + (motif_region.end - \
                            motif_region.start)/2
                tf_distances.append(atac_median - tf_median)

            # if we still haven't shifted past this peak...
            if motif_region.start <= (atac_median + H):
                try:
                    motif_region = next(motif_iter)   # get the next motif in the list
                    total_motif_sites += 1
                except StopIteration:
                    # No more TF motifs for this chromosome
                    break
                last_motif = motif_region
            else:
                motifs_within_region = False

    # Count any remaining TF motif sites after the last ATAC peak
    while(len(motif_region) > 0):
        try:
            motif_region = next(motif_iter)   # this gets the next motif in the list
        except StopIteration:
            break

    return [tf_distances, g_h, g_H, total_motif_sites]


def get_md_score(tf_motif_filename, mp_threads, atac_peaks_filename, CHROMOSOMES, radius):

    HISTOGRAM_BINS = 150
    pool = multiprocessing.Pool(mp_threads)
    results = pool.map(partial( find_motifs_in_chrom, \
                                files=[tf_motif_filename, atac_peaks_filename, radius]), \
                       CHROMOSOMES)
    pool.close()
    pool.join()

    results_matching_motif = [x for x in results if x is not None]
    if len(results_matching_motif) > 0:
        sums = np.sum(results_matching_motif, axis=0)
        overall_g_h = sums[1]
        overall_g_H = sums[2]
        overall_motif_sites = sums[3]

        # Calculate the heatmap for this motif's barcode
        #tf_distances = reduce(lambda a, b: [*a, *b], [x[0] for x in results if x is not None])
        tf_distances = sums[0]
        heatmap, xedges = np.histogram(tf_distances, bins=HISTOGRAM_BINS)
        str_heatmap = np.char.mod('%d', heatmap)
        # TODO: Use the motif sequences to generate a logo for each motif, based
        #       only on the overlapping ATAC-Seq peaks
        if overall_g_H >= 0:
            return [float(overall_g_h + 1)/(overall_g_H + 1), \
                    (overall_g_h + 1), \
                    (overall_g_H + 1), \
                    (overall_motif_sites + 1), \
                    ';'.join(str_heatmap)]
    else:
        return None

def main():
    parser = argparse.ArgumentParser(description='This script analyzes ATAC-Seq and GRO-Seq data and produces various plots for further data analysis.', epilog='IMPORTANT: Please ensure that ALL bed files used with this script are sorted by the same criteria.')
    parser.add_argument('-e', '--atac-peaks', dest='atac_peaks_filename', \
                        help='Full path to the ATAC-Seq broadPeak file.', \
                        default='', required=True)
    parser.add_argument('-m', '--motif-path', dest='tf_motif_path', \
                        help='Path to the location of the motif sites for the desired reference genome (i.e., "/usr/local/motifs/human/hg19/*").', \
                        default='', required=True)                        
    parser.add_argument('-r', '--radius', dest='radius', \
                        help='Radius around BED regions for which to scan for motifs. Default = 1500', \
                        default=1500, required=False, type=int)
    parser.add_argument('-g', '--genome', dest='genome', \
                        help='Genome to which the organism is mapped (e.g. hg38, mm10)', \
                        default='', required=False)
    parser.add_argument('-c', '--chromosomes', dest='chromosomes', \
                        help='Chromosome size file. See README for details in generating this file.', \
                        default='', required=False)    
    parser.add_argument('-t', '--threads', dest='mp_threads', \
                        help='Number of CPUs to use for multiprocessing of MD-score calculations. Depends on your hardware architecture.', \
                        default='1', required=False)
    parser.add_argument('-o', '--output', dest='output_dir', \
                        help='Path to where scores file will be saved. Save output will be your peak file rootname + _md_scores.txt.', \
                        default='', required=True)                        
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.3.0')                        

    args = parser.parse_args()

    #evaluation_radius = 750   # in bps
    ATAC_WIDTH_THRESHOLD = 5000   # in bp

    print('Starting --- ' + str(datetime.datetime.now()))
    atac_peaks_file = open(args.atac_peaks_filename)
    output_prefix = os.path.splitext(os.path.basename(args.atac_peaks_filename))[0]
    atac_csv_reader = csv.reader(atac_peaks_file, delimiter='\t')
    atac_line = next(atac_csv_reader)
    atac_peak_count = 0
    
    if os.path.exists(args.output_dir) is False:
        raise NameError("The specified output directory does not exist. Please specifiy a valid directory.")
    
    # Get chromosomes for multiprocessing
    if args.genome:
        chr_size_file = pybedtools.chromsizes(args.genome)
        unique_chr = list(chr_size_file.keys())[0:]
        CHROMOSOMES = [word for word in unique_chr if len(word) <= 6]
    elif args.chromosomes:
        chr_df = pd.read_csv(args.chromosomes, header=None, comment='#', sep="\t", usecols=[0,1], \
                                  names=['chrom', 'size'], na_filter=False, dtype={'chrom':'str', 'size':'int'})
        chr_list = list(chr_df['chrom'])
        CHROMOSOMES = [word for word in chr_list if len(word) <= 6]
    else:
        raise ValueError("Either -g/--genome or -c/--chromosomes argument must be specified.") 

    # skip the BedGraph headers
    while(atac_line[0][0] == '#'):
        atac_line = next(atac_csv_reader)

    # Determine the evaluation radius from the data itself
    all_widths = []
    while(atac_line):  # count the rest of ATAC peaks
        try:
            new_peak_width = int(atac_line[2]) - int(atac_line[1])
            if new_peak_width < ATAC_WIDTH_THRESHOLD:
                all_widths.append(new_peak_width)
            atac_line = next(atac_csv_reader)
        except StopIteration:
            break
        except IndexError:
            print("\nThere was an error with the ATAC-seq peaks file.\nPlease verify it conforms with a BedGraph-like format\n(tab-separated columns, any other lines commented with a leading '#')")
            sys.exit(1)
    atac_peak_mean = np.mean(all_widths)
    atac_peak_std = np.std(all_widths)
    evaluation_radius = (int(atac_peak_mean) + 2 * int(atac_peak_std)) / 2
    print('ATAC mean width: %d bp (std: %d bp). Determined an evaluation radius of %d bp' % \
          (atac_peak_mean, atac_peak_std, evaluation_radius))

    # Start reading from the top of the ATAC-Seq BedGraph again
    atac_peaks_file.seek(0)
    atac_csv_reader = csv.reader(atac_peaks_file, delimiter='\t')
    atac_line = next(atac_csv_reader)
    while(atac_line[0][0] == '#'):
        atac_line = next(atac_csv_reader)

    motif_stats = []
    sorted_motif_stats = []
    tf_motif_path = args.tf_motif_path + '/*'
    motif_filenames = glob.glob(tf_motif_path)
    motif_count = len(motif_filenames)
    print("Processing motif files in %s" % tf_motif_path)
    for filename in motif_filenames:
        filename_no_path = filename.split('/')[-1]
        if os.path.getsize(filename) > 0 and \
           os.path.basename(filename).endswith(tuple(['.bed', '.BedGraph', '.txt'])):
            [md_score, small_window, large_window, motif_site_count, heat] = get_md_score(filename, int(args.mp_threads), args.atac_peaks_filename, CHROMOSOMES, args.radius)
            print('The MD-score for the provided BED regions vs %s is %.6f' % (filename_no_path, md_score))
            motif_stats.append({ 'motif_file': filename_no_path, \
                                 'md_score': md_score, \
                                 'small_window': small_window, \
                                 'large_window': large_window, \
                                 'motif_site_count': motif_site_count, \
                                 'heat': heat })

    # sort the stats dictionary by MD-score, descending order
    sorted_motif_stats = sorted(motif_stats, key=itemgetter('md_score'), reverse=True)

    md_score_fp = open("%s/%s_md_scores.txt" % (args.output_dir, output_prefix), 'w')
    for stat in sorted_motif_stats:
        md_score_fp.write("%s,%s,%s,%s,%s,%s\n" % \
                          (stat['motif_file'], stat['md_score'], stat['small_window'], \
                           stat['large_window'], stat['motif_site_count'], stat['heat']))
    md_score_fp.close()

    print('All done --- %s' % str(datetime.datetime.now()))
    sys.exit(0)


if __name__=='__main__':
    main()
