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

#def main():
parser = argparse.ArgumentParser(description='This script produces an MA plot of TFs from ATAC-Seq data, for DMSO vs. treatment conditions.', \
    epilog="IMPORTANT: Please ensure that ALL files used with this script are sorted by the same criteria.\n\nExample:\nFor your files:\n     * mcf7_DMSO_md_scores.txt\n     * mcf7_Nutlin_md_scores.txt\n... you can use the following arguments to generate an MA plot with barcodes at a p-value cutoff of 1e-4:\n\n$ python differential_md_score.py -x mcf7 -1 DMSO -2 Nutlin -p 0.0001 -b\n\n", formatter_class=RawTextHelpFormatter)
parser.add_argument('-p', '--p-value', dest='p_value', metavar='P_VALUE', \
                    help='p-value cutoff to define which motifs to label in the MA plot. Defaults to 0.00001.', default=0.00001, required=False)
parser.add_argument('-1', '--assay-1', dest='assay_1', metavar='ASSAY_1', \
                    help='Conditions label for the reference/control assay of the differential pair (e.g. "DMSO", "control", "wildtype"). This will be the rootname before _md_scores.txt generated from process_atac. Used to find the proper file with the calculated MD-scores and on the plot labels.', required=True)
parser.add_argument('-2', '--assay-2', dest='assay_2', metavar='ASSAY_2', \
                    help='Conditions label for the perturbation assay of the differential pair (e.g., "doxycyclin", "p53_knockout"). This will be the rootname before _md_scores.txt generated from process_atac. Used to find the proper file with the calculated MD-scores and on the plot labels.', required=True)
parser.add_argument('-b', '--barcodes', dest='gen_barcode', action='store_true', \
                    help='Generate a barcode plot for each significant motif', default=False, required=False)
parser.add_argument('-o', '--output', dest='output_dir', \
                    help='Path to where output files will be saved.', \
                    default='', required=True)
parser.add_argument('-n', '--threads', dest='threads', metavar='THREADS', \
                    help='Number of threads for multi-processing. Defaults to 1.', default=1, required=False)    
args = parser.parse_args()

HISTOGRAM_BINS = 150
P_VALUE_CUTOFF = float(args.p_value)
assay_1_prefix = os.path.splitext(os.path.basename(args.assay_1))[0]
assay_2_prefix = os.path.splitext(os.path.basename(args.assay_2))[0]
threads = int(args.threads)

print('Starting --- ' + str(datetime.datetime.now()))

control_mds = {}
control_nr_peaks = {}
control_barcode = {}
labels = []
control_fd = open('%s' % args.assay_1)
for line in control_fd:
    line_chunks = line.split(',')
    if '.bed' in line_chunks[0]:
        control_mds[line_chunks[0][:-4]] = float(line_chunks[1])
        labels.append(line_chunks[0][:-4])
        control_nr_peaks[line_chunks[0][:-4]] = float(line_chunks[3])
        control_barcode[line_chunks[0][:-4]] = line_chunks[4]
perturbation_mds = {}
perturbation_nr_peaks = {}
perturbation_barcode = {}
perturbation_fd = open('%s' % args.assay_2)
for line in perturbation_fd:
    line_chunks = line.split(',')
    if '.bed' in line_chunks[0]:
        assert(line_chunks[0][:-4] in labels)
        perturbation_mds[line_chunks[0][:-4]] = float(line_chunks[1])
        perturbation_nr_peaks[line_chunks[0][:-4]] = int(line_chunks[3])
        perturbation_barcode[line_chunks[0][:-4]] = line_chunks[4]
        
def get_differential_md_scores(label):        
    control = float(control_mds[label])
    perturbation = float(perturbation_mds[label])
    nr_peaks = np.log2(float(control_nr_peaks[label]) + float(perturbation_nr_peaks[label]))
    fold_change = float(perturbation_mds[label]) - float(control_mds[label])
    p1 = float(abs(np.log10(control_mds[label])))
    p2 = float(abs(np.log10(perturbation_mds[label])))
    n1 = float(control_nr_peaks[label])
    n2 = float(perturbation_nr_peaks[label])
    if n1 <= 70:
        print('%s does not have a sufficient number of calls in %s. Setting MD score to 0.1 (background) for differntial analysis.' % (label, assay_1_prefix))
        p1 = .1
    if n2 <= 70:
        print('%s does not have a sufficient number of calls in %s. Setting MD score to 0.1 (background) for differntial analysis.' % (label, assay_2_prefix))
        p2 = .1
    if n1 >= 70:
        for line in control_barcode:
            control_bc_array = np.array(control_barcode[line].split(';'))
            control_bc_boot = control_bc_array.astype(int)
            #print(control_bc_boot)
            
            # configure bootstrap
            values = control_bc_boot
            control_n_iterations = 599
            control_n_size = int(len(control_bc_boot) * 0.1)
            # run bootstrap
            stats = list()
            for i in range(control_n_iterations):
                # prepare train and test sets
                train = resample(values, n_samples=control_n_size, replace=True)
                a = np.var(train)
                stats.append(a)
            control_bootstrap = np.median(stats)
            print(control_bootstrap)
                
    if n2 >=70:            
        for line in perturbation_barcode:
            perturbation_bc_array = np.array(perturbation_barcode[line].split(';'))
            perturbation_bc_boot = perturbation_bc_array.astype(int)
            
            # configure bootstrap
            values = perturbation_bc_boot
            perturbation_n_iterations = 599
            perturbation_n_size = int(len(perturbation_bc_boot) * 0.1)
            # run bootstrap
            stats = list()
            for i in range(perturbation_n_iterations):
                # prepare train and test sets
                train = resample(values, n_samples=perturbation_n_size, replace=True)
                a = np.var(train)
                stats.append(a)
            perturbation_bootstrap = np.median(stats)
            print(perturbation_bootstrap)
            
    if n1 <= 70:
        z_value = (abs(np.log10(p1)) - abs(np.log10(p2))) / np.sqrt((perturbation_bootstrap/perturbation_n_iterations))
    elif n2 <= 70:
        z_value = (abs(np.log10(p1)) - abs(np.log10(p2))) / np.sqrt((control_bootstrap/control_n_iterations))
    else: 
        z_value = (abs(np.log10(p1)) - abs(np.log10(p2))) / np.sqrt((control_bootstrap / n1) + (perturbation_bootstrap  / n2))
    p_value = norm.sf(abs(z_value))*2
    p_values.append(p_value)
    
    if p_value < (P_VALUE_CUTOFF / 10) and float(perturbation_mds[label]) > float(control_mds[label]):
        color = '#c64e50'
        size = 70
    elif p_value < P_VALUE_CUTOFF and float(perturbation_mds[label]) > float(control_mds[label]):
        color = 'maroon'
        size = 70
    elif p_value < (P_VALUE_CUTOFF / 10) and float(perturbation_mds[label]) < float(control_mds[label]):
        color = 'darkviolet'
        size = 70
    elif p_value < P_VALUE_CUTOFF and float(perturbation_mds[label]) < float(control_mds[label]):
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
                         'fold_change': fold_change, \
                         'control': control, \
                         'perturbation': perturbation, \
                         'nr_peaks': nr_peaks, \
                         'label': label}
    
    print('Done with %s.' % label)
    return differential_stats

with concurrent.futures.ProcessPoolExecutor(threads) as executor:
    jobs = [executor.submit(get_differential_md_scores, label) 
               for label in labels]
    differential_results = [r.result() for r in jobs]
    
sorted_differential_stats = sorted(differential_results, key=itemgetter('p_value'))

differential_stats_file = open("%s/%s_vs_%s_differential_md_scores.txt" \
                               % (args.output_dir, assay_1_prefix, assay_2_prefix), 'w')
for stat in sorted_differential_stats:
    differential_stats_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % \
                      (stat['motif_name'], stat['p_value'], stat['control_peaks'], \
                       stat['perturbation_peaks'], stat['control_md_score'], stat['perturbation_md_score']))   
differential_stats_file.close()

nr_peaks = [d['nr_peaks'] for d in sorted_differential_stats]
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
label_1_str = assay_1_prefix
label_2_str = assay_2_prefix

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
plt.savefig('%s/MA_%s_to_%s_md_score.png' % (args.output_dir, assay_1_prefix, assay_2_prefix), dpi=600)

#if args.gen_barcode:
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
    for row in range(int(HISTOGRAM_BINS/4)):
        heat_m[row] = control_bc_data
        ax0.matshow(heat_m, cmap=cm.YlGnBu)
    ax0.axis('off')
    ax0.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (control_nr_peaks[relevant_tf], control_mds[relevant_tf]), ha='center', size=18, zorder=0)
    ax0.text(HISTOGRAM_BINS/2, -10, assay_1_prefix, ha='center', size=18, zorder=0)

    perturbation_bc_data = np.array(perturbation_barcode[relevant_tf].split(';'))
    perturbation_bc_data = perturbation_bc_data.astype(float)
    heat_m = np.nan * np.empty(shape=(int(HISTOGRAM_BINS/4), HISTOGRAM_BINS))
    for row in range(int(HISTOGRAM_BINS/4)):
        heat_m[row] = perturbation_bc_data
        ax1.matshow(heat_m, cmap=cm.YlGnBu)
    ax1.axis('off')
    ax1.text(HISTOGRAM_BINS/2, HISTOGRAM_BINS/2, 'N(total) = %d\nMD-score = %.3f' % (perturbation_nr_peaks[relevant_tf], perturbation_mds[relevant_tf]), ha='center', size=18, zorder=0)
    ax1.text(HISTOGRAM_BINS/2, -10, assay_2_prefix, ha='center', size=18, zorder=0)

    plt.savefig('%s/%s_barcode_%s_vs_%s.png' % (args.output_dir, relevant_tf, assay_1_prefix, assay_2_prefix), dpi=600)

print('All done --- ' + str(datetime.datetime.now()))
sys.exit(0)


#if __name__=='__main__':
#    main()
