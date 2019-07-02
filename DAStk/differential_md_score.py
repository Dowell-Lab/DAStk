#! /usr/bin/env/python

from __future__ import print_function
import argparse
import datetime
import numpy as np
import sys
import os
import matplotlib as mpl
# to prevent display-related issues
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib import cm
from adjustText import adjust_text
from scipy.stats import norm
from argparse import RawTextHelpFormatter
from sklearn.utils import resample
from operator import itemgetter
from sklearn.utils import resample
from sklearn.tree import DecisionTreeClassifier
import concurrent.futures

# Usage:
#
# $ python differential_md_score.py -1 control -2 tamoxifen -p 0.000000000001 -o /outdir
#


def get_differential_md_scores(diff_params):
    (label, p1, p2, n1, n2, control_barcode, perturbation_barcode, is_chip, pval_cutoff) = diff_params
    control = p1
    perturbation = p2
    nr_peaks = n1 + n2
    perturbation_bc_array = np.array(perturbation_barcode.split(';'))
    control_bc_array = np.array(control_barcode.split(';'))

    if n1 <= 70:
        p1 = .1
    if n2 <= 70:
        p2 = .1
        
    fold_change = p2 - p1        

    if n1 > 70:
        control_bc_boot = control_bc_array.astype(int)
        # configure bootstrap
        values = control_bc_boot
        control_n_iterations = 1099
        control_n_size = int(len(control_bc_boot) * 0.5)
        # run bootstrap
        stats = list()
        for i in range(control_n_iterations):
            # prepare train and test sets
            train = resample(values, n_samples=control_n_size, replace=True)
            a = np.var(train)
            stats.append(a)
        control_bootstrap = abs(np.log(np.median(stats))) / 10

    if n2 > 70:
        perturbation_bc_boot = perturbation_bc_array.astype(int)
        # configure bootstrap
        values = perturbation_bc_boot
        perturbation_n_iterations = 1099
        perturbation_n_size = int(len(perturbation_bc_boot) * 0.5)
        # run bootstrap
        stats = list()
        for i in range(perturbation_n_iterations):
            # prepare train and test sets
            train = resample(values, n_samples=perturbation_n_size, replace=True)
            a = np.var(train)
            stats.append(a)
        perturbation_bootstrap = abs(np.log(np.median(stats))) / 10

    if (is_chip):
        if (n1 <= 70) & (n2 > 70):
            z_value = (p1 - p2) / np.sqrt(perturbation_bootstrap / n2)
        elif (n2 <= 70) & (n1 > 70):
            z_value = (p1 - p2) / np.sqrt(control_bootstrap / n1)
        elif (n1 <= 70) & (n2 <= 70):
            z_value = (p1 - p2)
            print('There were not enough regions to calculate the differential MDS for motif %s.' % label)
        else:
            z_value = (p1 - p2) / np.sqrt((control_bootstrap / n1) + (perturbation_bootstrap  / n2))

    else:
        if (n1 <= 70) & (n2 > 70):
            x1 = p1 * n1
            x2 = p2 * n2
            pooled = (x1 + x2)/(n1 + n2)
            z_value = (p1 - p2) / np.sqrt(pooled * (1 - pooled) * ((1/n1) + (1/n2)))
        elif (n2 <= 70) & (n1 > 70):
            x1 = p1 * n1
            x2 = p2 * n2
            pooled = (x1 + x2)/(n1 + n2)
            z_value = (p1 - p2) / np.sqrt(pooled * (1 - pooled) * ((1/n1) + (1/n2)))
        elif (n1 <= 70) & (n2 <= 70):
            z_value = (p1 - p2)
            print('There were not enough regions to calculate the differential MDS for motif %s.' % label)
        else:
            z_value = (p1 - p2) / np.sqrt((control_bootstrap / n1) + (perturbation_bootstrap  / n2))

    p_value = norm.sf(abs(z_value))*2

    # Get colors for MD plot based on p-value

    if p_value < (pval_cutoff / 10) and p2 > p1:
        color = '#c64e50'
        size = 70
    elif p_value < pval_cutoff and p2 > p1:
        color = 'maroon'
        size = 70
    elif p_value < (pval_cutoff / 10) and p2 < p1:
        color = 'darkviolet'
        size = 70
    elif p_value < pval_cutoff and p2 < p1:
        color = 'purple'
        size = 70
    else:
        color = '#4e74ae'
        size = 70

    differential_stats = { 'motif_name': label, \
                         'p_value': p_value, \
                         'control_peaks': int(n1), \
                         'perturbation_peaks': int(n2), \
                         'control_md_score': round(p1, 3), \
                         'perturbation_md_score': round(p2, 3), \
                         'color': color, \
                         'size': size, \
                         'fold_change': round(fold_change, 3), \
                         'control': control, \
                         'perturbation': perturbation, \
                         'nr_peaks': int(nr_peaks), \
                         'label': label}

    print('Done with %s.' % label)
    return differential_stats


def main():
    parser = argparse.ArgumentParser(description='This script produces an MA plot of TFs from ATAC-Seq data, for DMSO vs. treatment conditions.', \
        epilog="IMPORTANT: Please ensure that ALL files used with this script are sorted by the same criteria.\n\nExample:\nFor your files:\n     * mcf7_DMSO_md_scores.txt\n     * mcf7_Nutlin_md_scores.txt\n... you can use the following arguments to generate an MA plot with barcodes at a p-value cutoff of 1e-4:\n\n$ python differential_md_score.py -x mcf7 -1 DMSO -2 Nutlin -p 0.0001 -b\n\n", formatter_class=RawTextHelpFormatter)
    parser.add_argument('-p', '--p-value', dest='p_value', metavar='P_VALUE', \
                        help='p-value cutoff to define which motifs to label in the MA plot. Defaults to 0.00001.', default=0.00001, required=False)
    parser.add_argument('-1', '--assay-1', dest='assay_1', metavar='ASSAY_1', \
                        help='Control file generated from process_atac ending in the extension "md_scores.txt" (e.g. "DMSO", "control", "wildtype").', required=True)
    parser.add_argument('-2', '--assay-2', dest='assay_2', metavar='ASSAY_2', \
                        help='Perturbation file generated from process_atac ending in the extension "md_scores.txt" (e.g., "doxycyclin", "p53_knockout").', required=True)
    parser.add_argument('-m', '--label-1', dest='label_1', metavar='LABEL_1', \
                        help='Label for the MA plot title corresponding to assay 1', required=False)
    parser.add_argument('-n', '--label-2', dest='label_2', metavar='LABEL_2', \
                        help='Label for the MA plot title corresponding to assay 2', required=False)
    parser.add_argument('-b', '--barcodes', dest='gen_barcode', action='store_true', \
                        help='Generate a barcode plot for each significant motif', default=False, required=False)
    parser.add_argument('-o', '--output', dest='output_dir', \
                        help='Path to where output files will be saved.', \
                        default='', required=True)
    parser.add_argument('-t', '--threads', dest='threads', metavar='THREADS', \
                        help='Number of threads for multi-processing. Defaults to 1.', default=1, required=False)
    parser.add_argument('-c', '--chip', dest='chip', action='store_true', \
                        help='If the input is ChIP data, it may be useful to specify this flag as it will change the variance calulation because a large difference in sites between control and treatment will be expected.', default=False, required=False)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.3.0')   
    args = parser.parse_args()

    HISTOGRAM_BINS = 150
    assay_1_prefix = os.path.splitext(os.path.basename(args.assay_1))[0]
    if assay_1_prefix.endswith('_md_scores'):
        assay_1_prefix = assay_1_prefix[:-10]
    assay_2_prefix = os.path.splitext(os.path.basename(args.assay_2))[0]
    if assay_2_prefix.endswith('_md_scores'):
        assay_2_prefix = assay_2_prefix[:-10]    
    P_VALUE_CUTOFF = float(args.p_value)
    threads = int(args.threads)

    print('Starting --- ' + str(datetime.datetime.now()))

    control_mds = {}
    control_nr_peaks = {}
    control_barcodes = {}
    control_total_motifs = {}
    labels = []
    control_fd = open('%s' % args.assay_1)
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
    perturbation_fd = open('%s' % args.assay_2)
    for line in perturbation_fd:
        line_chunks = line.split(',')
        if '.bed' in line_chunks[0]:
            assert(line_chunks[0][:-4] in labels)
            perturbation_mds[line_chunks[0][:-4]] = float(line_chunks[1])
            perturbation_nr_peaks[line_chunks[0][:-4]] = round(float(line_chunks[3]))
            perturbation_total_motifs[line_chunks[0][:-4]] = int(line_chunks[4])
            perturbation_barcodes[line_chunks[0][:-4]] = line_chunks[5]


    with concurrent.futures.ProcessPoolExecutor(threads) as executor:
        jobs = [executor.submit(get_differential_md_scores, 
                                [   label, 
                                    float(control_mds[label]), \
                                    float(perturbation_mds[label]), \
                                    float(control_nr_peaks[label]), \
                                    float(perturbation_nr_peaks[label]), \
                                    control_barcodes[label], \
                                    perturbation_barcodes[label], \
                                    args.chip, \
                                    P_VALUE_CUTOFF ]
                               )
                for label in labels]
        differential_results = [r.result() for r in jobs]

    sorted_differential_stats = sorted(differential_results, key=itemgetter('p_value'))

    differential_stats_file = open("%s/%s_vs_%s_differential_md_scores.txt" \
                                % (args.output_dir, assay_1_prefix, assay_2_prefix), 'w')
    for stat in sorted_differential_stats:
        differential_stats_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
                        (stat['motif_name'], stat['p_value'], stat['control_peaks'], \
                        stat['perturbation_peaks'], stat['control_md_score'], stat['perturbation_md_score'], \
                        stat['fold_change']))
    differential_stats_file.close()

    nr_peaks = [np.log2(d['nr_peaks']) for d in sorted_differential_stats]
    fold_change = [d['fold_change'] for d in sorted_differential_stats]
    p_values = [d['p_value'] for d in sorted_differential_stats]
    colors = [d['color'] for d in sorted_differential_stats]
    sizes = [d['size'] for d in sorted_differential_stats]
    control = [d['control'] for d in sorted_differential_stats]
    perturbation = [d['perturbation'] for d in sorted_differential_stats]
    labels = [d['motif_name'] for d in sorted_differential_stats]

    np_control, np_perturbation, np_labels = np.array(control), np.array(perturbation), np.array(labels)

    # MA (mean/average) plot
    most_relevant_tfs = []
    plt.clf()
    fig, ax = plt.subplots()
    ax.scatter(nr_peaks, fold_change, s=sizes, edgecolor='white', linewidth=0.5, color=colors)
    texts = []

    for x, y, text, p_value in zip(nr_peaks, fold_change, np_labels, p_values):
        if p_value < P_VALUE_CUTOFF:
            print('%s (%.3f, p-value = %.2E)' % (text, y, p_value))
            label_color = 'maroon'
            if p_value < (P_VALUE_CUTOFF/10):
                label_color = '#c64e50'
            if y < 0:
                label_color = 'purple'
                if p_value < (P_VALUE_CUTOFF/10):
                    label_color = 'darkviolet'
            short_text = text.replace('HO_', '')
            short_text = short_text.replace('_HUMAN.H10MO', '')
            short_text = short_text.replace('_MOUSE.H10MO', '')
            short_text = short_text.split('_M', 1)[0]              
            texts.append(ax.text(x, y, u'%s' % short_text, fontsize=8, color=label_color))
            most_relevant_tfs.append(text)
    #adjust_text(texts, force_points=1, on_basemap=True, expand_points=(5,5), expand_text=(3,3), arrowprops=dict(arrowstyle="-", lw=1, color='grey', alpha=0.5))
    adjust_text(texts, force_points=1, expand_points=(2,2), expand_text=(2,2), arrowprops=dict(arrowstyle="-", lw=1, color='black', alpha=0.8))
    
    #Check to make sure plot titles are not too large... will skew plot image otherwise
    
    if len(assay_1_prefix) <= 19:
        label_1_str = assay_1_prefix
    else:
        label_1_str = assay_1_prefix[:19]
        
    if args.label_1:
        label_1_str = args.label_1
        
    if len(assay_2_prefix) <= 19:
        label_2_str = assay_2_prefix
    else:
        label_2_str = assay_2_prefix[:19]
        
    if args.label_2:
        label_2_str = args.label_2

    plt.title(u'MA for %s vs. %s MD-scores\n(p-value cutoff: %.2E)' % \
                (label_1_str, label_2_str, P_VALUE_CUTOFF), fontsize=12)
    plt.xlabel(u'$\log_2$(Sum #peaks overlapping 3kbp window)', fontsize=14)
    plt.ylabel(u'${\Delta}$ MD-score', fontsize=14)
    plt.xlim(np.min(nr_peaks), np.max(nr_peaks) + 1)
    plt.xscale('log',basex=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x',reset=False,which='both',length=5,width=1)
    y_bound = max(np.abs(np.min(fold_change)), np.max(fold_change)) + 0.01
    plt.ylim(-1 * y_bound, y_bound)
    plt.tight_layout()
    plt.savefig('%s/MA_%s_to_%s_md_score.png' % (args.output_dir, assay_1_prefix, assay_2_prefix), dpi=600)

    if args.gen_barcode:
        # Generate barcodes for each relevant TF, for both conditions
        print('Generating barcode plots for significant motifs...')
        np.seterr(invalid='ignore')
        for relevant_tf in most_relevant_tfs:
            plt.clf()
            fig, (ax0, ax1) = plt.subplots(ncols=2)
            fig.suptitle('Barcode plots for %s' % relevant_tf)

            control_bc_data = np.array(control_barcodes[relevant_tf].split(';'))
            # Condense the barcode to half the bins for prettier display
            control_bc_data = control_bc_data.astype(int)
            heat_m = np.nan * np.empty(shape=(int(HISTOGRAM_BINS/4), HISTOGRAM_BINS))
            for row in range(int(HISTOGRAM_BINS/4)):
                heat_m[row] = control_bc_data
                ax0.matshow(heat_m, cmap=cm.YlGnBu)
            ax0.axis('off')
            ax0.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (control_nr_peaks[relevant_tf], control_mds[relevant_tf]), ha='center', size=18, zorder=0)
            ax0.text(HISTOGRAM_BINS/2, -10, label_1_str, ha='center', size=18, zorder=0)

            perturbation_bc_data = np.array(perturbation_barcodes[relevant_tf].split(';'))
            perturbation_bc_data = perturbation_bc_data.astype(float)
            heat_m = np.nan * np.empty(shape=(int(HISTOGRAM_BINS/4), HISTOGRAM_BINS))
            for row in range(int(HISTOGRAM_BINS/4)):
                heat_m[row] = perturbation_bc_data
                ax1.matshow(heat_m, cmap=cm.YlGnBu)
            ax1.axis('off')
            ax1.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (perturbation_nr_peaks[relevant_tf], perturbation_mds[relevant_tf]), ha='center', size=18, zorder=0)
            ax1.text(HISTOGRAM_BINS/2, -10, label_2_str, ha='center', size=18, zorder=0)

            plt.tight_layout()
            plt.savefig('%s/%s_barcode_%s_vs_%s.png' % (args.output_dir, relevant_tf, assay_1_prefix, assay_2_prefix), dpi=600)

        print('All done --- ' + str(datetime.datetime.now()))
        sys.exit(0)


if __name__=='__main__':
    main()
