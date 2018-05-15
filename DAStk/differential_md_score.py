#! /usr/bin/env/python

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

# Usage:
#
# $ python differential_md_score.py -x brg1fl -1 control -2 tamoxifen -p 0.000000000001
#

def main():
    parser = argparse.ArgumentParser(description='This script produces an MA plot of TFs from ATAC-Seq data, for DMSO vs. treatment conditions.', \
        epilog="IMPORTANT: Please ensure that ALL files used with this script are sorted by the same criteria.\n\nExample:\nFor your files:\n     * mcf7_DMSO_md_scores.txt\n     * mcf7_Nutlin_md_scores.txt\n... you can use the following arguments to generate an MA plot with barcodes at a p-value cutoff of 1e-4:\n\n$ python differential_md_score.py -x mcf7 -1 DMSO -2 Nutlin -p 0.0001 -b\n\n", formatter_class=RawTextHelpFormatter)
    parser.add_argument('-x', '--prefix', dest='output_prefix', metavar='CELL_TYPE', \
                        help='Cell type (k562, imr90, etc), or any other appropriate output file prefix', required=True)
    parser.add_argument('-p', '--p-value', dest='p_value', metavar='P_VALUE', \
                        help='p-value cutoff to define which motifs to label in the MA plot. Defaults to 0.00001.', default=0.00001, required=False)
    parser.add_argument('-1', '--assay-1', dest='assay_1', metavar='ASSAY_1', \
                        help='Conditions label for the reference/control assay of the differential pair (e.g. "DMSO", "control", "wildtype"). Used to find the proper file with the calculated MD-scores and on the plot labels.', required=True)
    parser.add_argument('-2', '--assay-2', dest='assay_2', metavar='ASSAY_2', \
                        help='Conditions label for the perturbation assay of the differential pair (e.g., "doxycyclin", "p53_knockout"). Used to find the proper file with the calculated MD-scores and on the plot labels.', required=True)
    parser.add_argument('-b', '--barcodes', dest='gen_barcode', action='store_true', \
                        help='Generate a barcode plot for each significant motif', default=False, required=False)
    args = parser.parse_args()

    HISTOGRAM_BINS = 100
    P_VALUE_CUTOFF = float(args.p_value)

    print('Starting --- ' + str(datetime.datetime.now()))

    control_mds = {}
    control_nr_peaks = {}
    control_barcode = {}
    labels = []
    control_fd = open('%s_%s_md_scores.txt' % (args.output_prefix, args.assay_1))
    for line in control_fd:
        line_chunks = line.split(',')
        if '.bed' in line_chunks[0]:
            control_mds[line_chunks[0][:-4]] = float(line_chunks[1])
            labels.append(line_chunks[0][:-4])
            control_nr_peaks[line_chunks[0][:-4]] = int(line_chunks[3])
            control_barcode[line_chunks[0][:-4]] = line_chunks[5]
    perturbation_mds = {}
    perturbation_nr_peaks = {}
    perturbation_barcode = {}
    perturbation_fd = open('%s_%s_md_scores.txt' % (args.output_prefix, args.assay_2))
    for line in perturbation_fd:
        line_chunks = line.split(',')
        if '.bed' in line_chunks[0]:
            assert(line_chunks[0][:-4] in labels)
            perturbation_mds[line_chunks[0][:-4]] = float(line_chunks[1])
            perturbation_nr_peaks[line_chunks[0][:-4]] = int(line_chunks[3])
            perturbation_barcode[line_chunks[0][:-4]] = line_chunks[5]

    print('Done gathering data, ready to plot --- ' + str(datetime.datetime.now()))

    control = []
    perturbation = []
    nr_peaks = []
    fold_change = []
    p_values = []
    colors = []
    sizes = []
    for label in sorted(labels):
        if control_mds[label] == 0 and perturbation_mds[label] == 0:
            # Skip this motif, it will only skew the graph and won't provide any useful info
            continue

        control.append(float(control_mds[label]))
        perturbation.append(float(perturbation_mds[label]))
        nr_peaks.append(np.log2(float(control_nr_peaks[label]) + float(perturbation_nr_peaks[label])))
        fold_change.append(float(perturbation_mds[label]) - float(control_mds[label]))
        p1 = float(control_mds[label])
        p2 = float(perturbation_mds[label])
        n1 = float(control_nr_peaks[label])
        n2 = float(perturbation_nr_peaks[label])
        x1 = p1 * n1
        x2 = p2 * n2
        if n1 == 0:
            print('%s had an MD-score of 0 in %s' % (label, args.assay_1))
            n1 = 1
        if n2 == 0:
            print('%s had an MD-score of 0 in %s' % (label, args.assay_2))
            n2 = 1
        pooled = (x1 + x2)/(n1 + n2)
        z_value = (p1 - p2) / np.sqrt(pooled * (1 - pooled) * ((1/n1) + (1/n2)))
        p_value = norm.sf(abs(z_value))*2
        p_values.append(p_value)
        if p_value < (P_VALUE_CUTOFF / 10) and float(perturbation_mds[label]) > float(control_mds[label]):
            colors.append('#c64e50')
            sizes.append(70)
        elif p_value < P_VALUE_CUTOFF and float(perturbation_mds[label]) > float(control_mds[label]):
            colors.append('maroon')
            sizes.append(70)
        elif p_value < (P_VALUE_CUTOFF / 10) and float(perturbation_mds[label]) < float(control_mds[label]):
            colors.append('darkviolet')
            sizes.append(70)
        elif p_value < P_VALUE_CUTOFF and float(perturbation_mds[label]) < float(control_mds[label]):
            colors.append('purple')
            sizes.append(70)
        else:
            colors.append('#4e74ae')
            sizes.append(70)


    np_control, np_perturbation, np_labels = np.array(control), np.array(perturbation), np.array(sorted(labels))

    # MA (mean/average) plot
    most_relevant_tfs = []
    plt.clf()
    fig, ax = plt.subplots()
    ax.scatter(nr_peaks, fold_change, s=sizes, edgecolor='white', linewidth=0.5, color=colors)
    texts = []
    for x, y, text, p_value in zip(nr_peaks, fold_change, np_labels, p_values):
        if p_value < P_VALUE_CUTOFF:
        #if (y > 0.02 or y < -0.02) and x > 12:
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
            texts.append(ax.text(x, y, u'%s' % short_text, fontsize=8, color=label_color))
            most_relevant_tfs.append(text)
    #adjust_text(texts, force_points=1, on_basemap=True, expand_points=(5,5), expand_text=(3,3), arrowprops=dict(arrowstyle="-", lw=1, color='grey', alpha=0.5))
    adjust_text(texts, force_points=1, expand_points=(2,2), expand_text=(2,2), arrowprops=dict(arrowstyle="-", lw=1, color='black', alpha=0.8))
    plt.title(u'MA for %s vs. %s MD-scores\n(%s, p-value cutoff: %.2E)' % \
              (args.assay_1, args.assay_2, args.output_prefix, P_VALUE_CUTOFF), fontsize=12)
    plt.xlabel(u'$\log_2$(Sum #peaks overlapping 3kbp window)', fontsize=14)
    plt.ylabel(u'${\Delta}$ MD-score', fontsize=14)
    plt.xlim(np.min(nr_peaks), np.max(nr_peaks) + 1)
    plt.xscale('log',basex=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x',reset=False,which='both',length=5,width=1)
    y_bound = max(np.abs(np.min(fold_change)), np.max(fold_change)) + 0.01
    plt.ylim(-1 * y_bound, y_bound)
    plt.savefig('%s_MA_%s_to_%s_md_score.png' % (args.output_prefix, args.assay_1, args.assay_2), dpi=600)

    if args.gen_barcode:
        # Generate barcodes for each relevant TF, for both conditions
        print('Generating barcode plots for significant motifs...')
        for relevant_tf in most_relevant_tfs:
            plt.clf()
            plt.title('Barcode plots for %s' % relevant_tf)
            fig, (ax0, ax1) = plt.subplots(ncols=2)

            control_bc_data = np.array(control_barcode[relevant_tf].split(';'))
            # Condense the barcode to half the bins for prettier display
            control_bc_data = control_bc_data.astype(int)
            heat_m = np.nan * np.empty(shape=(int(HISTOGRAM_BINS/4), HISTOGRAM_BINS))
            for row in range(HISTOGRAM_BINS/4):
                heat_m[row] = control_bc_data
                ax0.matshow(heat_m, cmap=cm.YlGnBu)
            ax0.axis('off')
            ax0.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (control_nr_peaks[relevant_tf], control_mds[relevant_tf]), ha='center', size=18, zorder=0)
            ax0.text(HISTOGRAM_BINS/2, -10, args.assay_1, ha='center', size=18, zorder=0)

            perturbation_bc_data = np.array(perturbation_barcode[relevant_tf].split(';'))
            perturbation_bc_data = perturbation_bc_data.astype(int)
            heat_m = np.nan * np.empty(shape=(int(HISTOGRAM_BINS/4), HISTOGRAM_BINS))
            for row in range(HISTOGRAM_BINS/4):
                heat_m[row] = perturbation_bc_data
                ax1.matshow(heat_m, cmap=cm.YlGnBu)
            ax1.axis('off')
            ax1.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (perturbation_nr_peaks[relevant_tf], perturbation_mds[relevant_tf]), ha='center', size=18, zorder=0)
            ax1.text(HISTOGRAM_BINS/2, -10, args.assay_2, ha='center', size=18, zorder=0)

            plt.savefig('%s_%s_barcode_%s_vs_%s.png' % (args.output_prefix, relevant_tf, args.assay_1, args.assay_2), dpi=600)

    print('All done --- ' + str(datetime.datetime.now()))
    sys.exit(0)


if __name__=='__main__':
    main()
