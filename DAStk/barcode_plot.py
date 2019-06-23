from __future__ import print_function
import argparse
import datetime
import numpy as np
import sys
import matplotlib as mpl
# to prevent display-related issues
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib import cm
from adjustText import adjust_text
from scipy.stats import norm
from argparse import RawTextHelpFormatter

def main():
    parser = argparse.ArgumentParser(description='This script generates barcodes using the output from DAStk.', \
                        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-md', '--scores-1', dest='scores_1', metavar='ASSAY_1', \
                    help='First MD score file generated from DAStk.', required=True, type=str)
    parser.add_argument('-MD', '--scores-2', dest='scores_2', metavar='ASSAY_1', \
                    help='Second MD score file generated from DAStk. Not required if argument \'-s/--single\' is specified.', required=False, type=str)
    parser.add_argument('-a', '--assay-1', dest='assay_1', metavar='ASSAY_1', \
                    help='Assay name of first MD score file.', required=True, type=str)
    parser.add_argument('-A', '--assay-2', dest='assay_2', metavar='ASSAY_2', \
                    help='Assay name of second MD score file. Will be title of barcode plot. Not required if argument \'-s/--single\' is specified.', required=False, type=str)
    parser.add_argument('-tf', '--transcription-factor', dest='relevant_tf', metavar='TF', \
                    help='Transcription factor you would like to plot. Should be full prefix before the .bed extension printed in the MD score output (e.g. JUND_HUMAN.H11MO.0.A).', required=True, type=str)
    parser.add_argument('-o', '--output', dest='output_dir', metavar='OUTPUT_DIR', \
                    help='Path to directory where plot will be saved.', required=True, type=str)
    parser.add_argument('-s', '--single', dest='single', action='store_true', \
                    help='Generate a single barcode rather than a side-by-side comparison.', default=False, required=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.3.0')       
    args = parser.parse_args()

    if not args.single:
        if args.scores_2 is None or args.assay_2 is None:
            parser.error('If argument single is not specified, you must provide a second MD score file and assay name.')
            print(args)
            return 0
    
    control = args.scores_1
    perturbation = args.scores_2
    relevant_tf = args.relevant_tf
    HISTOGRAM_BINS = 150
    assay_1 = args.assay_1
    assay_2 = args.assay_2
    output_dir = args.output_dir
    
    print('Parsing MD score file...' + str(datetime.datetime.now()))

    if args.single:
        control_mds = {}
        control_nr_peaks = {}
        control_barcodes = {}
        control_total_motifs = {}
        labels = []
        control_fd = open('%s' % assay_1)
        for line in control_fd:
            line_chunks = line.split(',')
            if '.bed' in line_chunks[0]:
                control_mds[line_chunks[0][:-4]] = float(line_chunks[1])
                labels.append(line_chunks[0][:-4])
                control_nr_peaks[line_chunks[0][:-4]] = round(float(line_chunks[3]))
                control_total_motifs[line_chunks[0][:-4]] = int(line_chunks[4])
                control_barcodes[line_chunks[0][:-4]] = line_chunks[5]
        
    else:    
        control_mds = {}
        control_nr_peaks = {}
        control_barcodes = {}
        control_total_motifs = {}
        labels = []
        control_fd = open('%s' % assay_1)
        for line in control_fd:
            line_chunks = line.split(',')
            if '.bed' in line_chunks[0]:
                control_mds[line_chunks[0][:-4]] = float(line_chunks[1])
                labels.append(line_chunks[0][:-4])
                control_nr_peaks[line_chunks[0][:-4]] = round(float(line_chunks[3]))
                control_total_motifs[line_chunks[0][:-4]] = int(line_chunks[4])
                control_barcodes[line_chunks[0][:-4]] = line_chunks[5]
        perturbation_mds = {}
        perturbation_nr_peaks = {}
        perturbation_barcodes = {}
        perturbation_total_motifs = {}
        perturbation_fd = open('%s' % assay_2)
        for line in perturbation_fd:
            line_chunks = line.split(',')
            if '.bed' in line_chunks[0]:
                assert(line_chunks[0][:-4] in labels)
                perturbation_mds[line_chunks[0][:-4]] = float(line_chunks[1])
                perturbation_nr_peaks[line_chunks[0][:-4]] = round(float(line_chunks[3]))
                perturbation_total_motifs[line_chunks[0][:-4]] = int(line_chunks[4])
                perturbation_barcodes[line_chunks[0][:-4]] = line_chunks[5]
                    
    print('Plotting barcode(s)... ' + str(datetime.datetime.now()))            
            
    if args.single:            
        plt.clf()
        plt.title('Barcode plots for %s' % relevant_tf)
        fig, (ax0) = plt.subplots(ncols=1)
        
        control_bc_data = np.array(control_barcode[relevant_tf].split(';'))
        # Condense the barcode to half the bins for prettier display
        control_bc_data = control_bc_data.astype(int)
        heat_m = np.nan * np.empty(shape=(int(HISTOGRAM_BINS/4), HISTOGRAM_BINS))
        for row in range(int(HISTOGRAM_BINS/4)):
            heat_m[row] = control_bc_data
            ax0.matshow(heat_m, cmap=cm.YlGnBu)
        ax0.axis('off')
        ax0.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (control_nr_peaks[relevant_tf], control_mds[relevant_tf]), ha='center', size=12, zorder=0)
        ax0.text(HISTOGRAM_BINS/2, -5, assay_1, ha='center', size=12, zorder=0)
        plt.savefig('%s/%s_%s.png' % (output_dir, relevant_tf ,assay_1) , dpi=600)
        
    else:
        plt.clf()
        plt.title('Barcode plots for %s' % relevant_tf)
        fig, (ax0, ax1) = plt.subplots(ncols=2)
        
        control_bc_data = np.array(control_barcode[relevant_tf].split(';'))
        # Condense the barcode to half the bins for prettier display
        control_bc_data = control_bc_data.astype(int)
        heat_m = np.nan * np.empty(shape=(int(HISTOGRAM_BINS/4), HISTOGRAM_BINS))
        for row in range(int(HISTOGRAM_BINS/4)):
            heat_m[row] = control_bc_data
            ax0.matshow(heat_m, cmap=cm.YlGnBu)
        ax0.axis('off')
        ax0.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (control_nr_peaks[relevant_tf], control_mds[relevant_tf]), ha='center', size=12, zorder=0)
        ax0.text(HISTOGRAM_BINS/2, -5, assay_1, ha='center', size=12, zorder=0)
        
        perturbation_bc_data = np.array(perturbation_barcode[relevant_tf].split(';'))
        perturbation_bc_data = perturbation_bc_data.astype(int)
        heat_m = np.nan * np.empty(shape=(int(HISTOGRAM_BINS/4), HISTOGRAM_BINS))
        for row in range(int(HISTOGRAM_BINS/4)):
            heat_m[row] = perturbation_bc_data
            ax1.matshow(heat_m, cmap=cm.YlGnBu)
        ax1.axis('off')
        ax1.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (perturbation_nr_peaks[relevant_tf], perturbation_mds[relevant_tf]), ha='center', size=12, zorder=0)
        ax1.text(HISTOGRAM_BINS/2, -5, assay_2, ha='center', size=12, zorder=0)
        
        plt.savefig('%s/%s_%s_vs_%s.png' % (output_dir, relevant_tf ,assay_1, assay_2) , dpi=600)
        
    print('Done! ' + str(datetime.datetime.now()))   
        
if __name__=='__main__':
    main()