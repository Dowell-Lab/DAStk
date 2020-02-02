from __future__ import print_function
import argparse
import datetime
import numpy as np
import sys
import matplotlib as mpl
# to prevent display-related issues
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
plt.ioff()
from matplotlib import cm
from adjustText import adjust_text
from scipy.stats import norm
from argparse import RawTextHelpFormatter

def main():
    parser = argparse.ArgumentParser(description='This script generates barcodes using the output from DAStk.', \
                        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-s', '--stats', dest='stats', metavar='STATS_FILE', \
                    help='First MD score file generated from DAStk.', required=True, type=str)
    parser.add_argument('-m', '--label-1', dest='label_1', metavar='LABEL_1', \
                        help='Label for the MA plot title corresponding to assay 1', required=True)
    parser.add_argument('-n', '--label-2', dest='label_2', metavar='LABEL_2', \
                        help='Label for the MA plot title corresponding to assay 2', required=True)
    parser.add_argument('-w', '--window', dest='window_size', metavar='WINDOW', \
                        help='Label for the MA plot title corresponding to window size (str). Default = \'3kb\'', type=str, default='3kb', required=False)    
    #parser.add_argument('-tf', '--transcription-factor', nargs='+', dest='relevant_tfs', metavar='TFs', \
    #                help='Transcription factor you would like to plot. Should be full prefix before the .bed extension printed in the MD score output (e.g. JUND_HUMAN.H11MO.0.A).', required=False)      
    parser.add_argument('-p', '--p-value', dest='p_value', metavar='P_VALUE', \
                        help='p-value cutoff to define which motifs to label in the MA plot. Defaults to 0.00001.', default=0.00001, required=False)
    parser.add_argument('-l', '--label_p-value', dest='label_by_p_value', action='store_true', \
                        help='Label all TFs falling below the specified p-value cutoff.', required=False)
    parser.add_argument('-o', '--output', dest='output_dir', metavar='OUTPUT_DIR', \
                    help='Path to directory where plot will be saved.', required=True, type=str)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0')       
    args = parser.parse_args()
    
    output_dir = args.output_dir
    P_VALUE_CUTOFF = float(args.p_value)
    
    print('Parsing MD score stats...' + str(datetime.datetime.now()))
    
    nr_peaks = []
    delta_md = []
    p_values = []
    labels = []
    sizes=[]
    colors=[]
    stats_file = open('%s' % args.stats)
    for line in stats_file:
        line_chunks = line.split('\t')
        if (line_chunks[7] != '0.0\n') & (line_chunks[7] != '0\n') :
            labels.append(line_chunks[0])
            p_value = float(line_chunks[1])
            p_values.append(float(line_chunks[1]))
            nr_peaks.append(np.log2(int(line_chunks[3]) + int(line_chunks[4])))
            p1 = float(line_chunks[5])
            p2 = float(line_chunks[6])
            delta_md.append(float(line_chunks[7]))
        
            if p_value < (P_VALUE_CUTOFF / 10) and p2 > p1:
                colors.append('#c64e50')
                sizes.append(70)
            elif p_value < P_VALUE_CUTOFF and p2 > p1:
                colors.append('maroon')
                sizes.append(70)
            elif p_value < (P_VALUE_CUTOFF / 10) and p2 < p1:
                colors.append('darkviolet')
                sizes.append(70)
            elif p_value < P_VALUE_CUTOFF and p2 < p1:
                colors.append('purple')
                sizes.append(70)
            else:
                colors.append('#4e74ae')
                sizes.append(70)
                    
    # MA (mean/average) plot
    #most_relevant_tfs = []    
    plt.clf()
    fig, ax = plt.subplots()
    ax.scatter(nr_peaks, delta_md, s=sizes, edgecolor='white', linewidth=0.5, color=colors)
    #relevant_tfs = [args.relevant_tfs]
    texts = []

    for x, y, text, p_value, color in zip(nr_peaks, delta_md, labels, p_values, colors):
        #if args.relevant_tfs:
        #    short_text = text.replace('HO_', '')
        #    short_text = short_text.replace('_HUMAN.H10MO', '')
        #    short_text = short_text.split('_M', 1)[0]
        #    texts.append(ax.text(x, y, u'%s' % short_text, fontsize=8, color=color))            
            
        if p_value < P_VALUE_CUTOFF:
            print('%s (%.3f, p-value = %.2E)' % (text, y, p_value))
            if args.label_by_p_value:
                short_text = text.replace('HO_', '')
                short_text = short_text.replace('_HUMAN.H10MO', '')
                short_text = short_text.split('_M', 1)[0]
                short_text = short_text.split('_', 1)[0]
                texts.append(ax.text(x, y, u'%s' % short_text, fontsize=8, color=color))
    adjust_text(texts, force_points=1, expand_points=(2,2), expand_text=(2,2), arrowprops=dict(arrowstyle="-", lw=1, color='black', alpha=0.8))

    label_1_str = args.label_1
    label_2_str = args.label_2
    plt.title(u'MA for %s vs. %s MD-scores\n(p-value cutoff: %.2E)' % \
                (label_1_str, label_2_str, P_VALUE_CUTOFF), fontsize=12)
    plt.xlabel(u'$\log_2$(Sum #peaks overlapping ' + args.window_size + ' window)', fontsize=14)
    plt.ylabel(u'${\Delta}$ MD-score', fontsize=14)
    plt.xlim(np.min(nr_peaks), np.max(nr_peaks) + 1)
    plt.xscale('log',basex=2)
    loc = plticker.MultipleLocator(base=2.0) # this locator puts ticks at regular intervals
    ax.xaxis.set_major_locator(loc)     
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x',reset=False,which='both',length=4,width=1)
    y_bound = max(np.abs(np.min(delta_md)), np.max(delta_md)) + 0.01
    plt.ylim(-1 * y_bound, y_bound)
    plt.tight_layout()
    plt.savefig('%s/MA_%s_to_%s_md_score.png' % (args.output_dir, label_1_str, label_2_str), dpi=600)          

if __name__=='__main__':
    main()    
