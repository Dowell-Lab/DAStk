#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import numpy as np
import io
import argparse
import os
import sys
import datetime
import upsetplot
import matplotlib.colors as mc
import colorsys
from upsetplot import UpSet
from collections import defaultdict
from functools import reduce
from itertools import chain
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted, venn2
from io import StringIO
from argparse import RawTextHelpFormatter

######## For adjusting darkness/lightness in venn diagrams

def adjust_lightness(color, amount=0.5):
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

######## For merging multiple dictionaries on a matching key

def merge(d1, d2, key='motif_key'):
    r = defaultdict(list)

    for k, v in chain(d1.items(), d2.items()):
        if k != key:
            r[k].extend(v if isinstance(v, list) else [v])

    return {**r, key: d1[key]}

####### Function to "explode" Pandas dataframes. As of v.0.25.3, this is now available as a function, but is not yet super fine-tunable

def explode(df, lst_cols, fill_value='', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values    
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:        
        res = res.reset_index(drop=True)
    return res

########## List of arguments for plotting

def main():

    parser = argparse.ArgumentParser(description='Venn Diagram Generator\n\n===============================================================================================================\n\nGenerates a Venn Diagram for lists of significant TFs coming out of the DAStk differential MD score stats file results.', usage='%(prog)s --stats <stats files> --rootname <rootname> --output /path/to/out' , formatter_class=RawTextHelpFormatter)
    
    required = parser.add_argument_group('Required Arguments')
    optional = parser.add_argument_group('Optional Arguments')
    
    required.add_argument('-s', '--stats', dest='stats_files', nargs='+', \
                        help='Full path to stats files (minimum of two) separated by a space. Limit of 9. EITHER results from differential_md_score OR process_atac may be used, but file types should not be combined. File type will be determined based on presense or absense of \'differential\' from the filename (included as an extension by default in the differential_md_score stats file output).', required=True)

    required.add_argument('-r', '--rootname', dest='rootname', metavar='<ROOTNAME>', \
                        help='Rootname for saving plots.', required=True)
    
    required.add_argument('-o', '--output', dest='output', metavar='/path/to/out/dir', \
                        help='Directory where plots will be saves.', required=True)
    
    optional.add_argument('-uw', '--unweighted', dest='unweighted', action='store_true', \
                        help='Produce an unweighted (vs. propotional, default) venn diagram. Default = False', default=False, required=False)
    
    optional.add_argument('-sig', '--significant', dest='significant', action='store_true', \
                        help='Intersect motifs which are significant (differential MD score stats files only). Default = False.', default=False, required=False)
    
    optional.add_argument('-e', '--enriched', dest='enriched', action='store_true', \
                        help='Intersect motifs which are enriched (either differentially or raw MD scores). Default = False', default=False, required=False)
    
    optional.add_argument('-d', '--depleted', dest='depleted', action='store_true', \
                        help='Intersect motifs which are depleted (either differentially or raw MD scores). Default = False', default=False, required=False)    
    
    optional.add_argument('-p', '--pvalue', dest='pval_threshold', metavar='<PVALUE>', \
                       help='p-value threshold for significant values which will be plotted. Default = pe-6.', default=.000001, required=False, type=float)
    
    optional.add_argument('-m', '--md-score-threshold', dest='md_score_threshold', metavar='<PVALUE>', \
                       help='Differential MD score threshold for values which will be plotted (positive OR negative). Default = 0.1.', default=0.1, required=False, type=float)
    
    optional.add_argument('-md', '--depleted-threshold', dest='depleted_threshold', metavar='<PVALUE>', \
                       help='Threshold for depletion raw MD score files. Default = 0.08.', default=0.08, required=False, type=float)
    
    optional.add_argument('-me', '--enriched-threshold', dest='enriched_threshold', metavar='<PVALUE>', \
                       help='Threshold for enrichment raw MD score files. Default = 0.2.', default=0.2, required=False, type=float)      
    
    optional.add_argument('-l', '--labels', dest='plot_labels', metavar='<PLOT_LABELS>', nargs='+',  \
                       help='Plot labels for files provided. Number of labels provided must match the number of files provided.', required=False, type=str)
    
    optional.add_argument('-c', '--colors', dest='colors', metavar='<COLORS>', nargs='+',  \
                       help='Hex colors for files provided. Number of labels provided must match the number of files provided for venn diagrams (2 or 3 files) or equal 2 for the upset catplots (4 or more files).', required=False, type=str)       

    args = parser.parse_args()
    
################################### SET GLOBAL VARIABLES AND ERROR CHECK ######################################
    PVAL_THRESH = args.pval_threshold
    DIFF_MD_THRESH = args.md_score_threshold
    ENRICHED_THRESH = args.enriched_threshold
    DEPLETED_THRESH = args.depleted_threshold
    
    pd.options.mode.chained_assignment = None
    
    ### Number of files and provided files have to match if less than 12 files are provided and plotting is not disabled
    if (len(args.stats_files) != len(args.plot_labels)) & (len(args.stats_files) <= 12):
        raise ValueError("Number of stats files provied and number of labels provided do not match -- please adjust these values such that they are equal. \n" + str(datetime.datetime.now()))
    
    try:
        for file in args.stats_files:
            if os.stat(file).st_size < 0:
                raise ValueError("One or more of the stats files provided is empty. Please check your specified files.")
            else:
                continue
    except OSError:
        raise ValueError("Stats file(s) missing ... exiting. \n" + str(datetime.datetime.now()))
        
###############################################################################################################
## Identify which stats files are being provided and parse according to either p-value or enriched / depleted
# Will determine file type based on presence / absense of 'differential' from file name
###############################################################################################################
    print("Beginning processing..." + str(datetime.datetime.now()))
    motif_lists = []
    data_dict_list = []
    for file in args.stats_files :
        if ('differential' not in file):
            print('Running intersecctions on raw MD score files.')
            file_columns = ['motif', 'md_score', 'h', 'H', 'total_motif_hits']
            stats_df = pd.read_csv(file, header=None, usecols=range(5), \
                             names=file_columns)
            stats_df['motif_key'] = stats_df['motif'].str.split('.').str[0]
            if args.enriched &~ args.depleted:
                stats_df_sig =  stats_df[(stats_df['md_score'] > ENRICHED_THRESH) & (stats_df['H'] >= 70 )]
            elif args.depleted &~ args.enriched:
                stats_df_sig =  stats_df[(stats_df['md_score'] < DEPLETED_THRESH) & (stats_df['H'] >= 70 )]
            elif args.enriched & args.depleted:
                stats_df_sig =  stats_df[((stats_df['md_score'] > ENRICHED_THRESH) | (stats_df['md_score'] < DEPLETED_THRESH)) & (stats_df['H'] >= 70 )]       
            else:
                raise ValueError("User must specify enriched, depleted, or both for the intersection type.")
                
            motif_lists.append(stats_df_sig['motif_key'].tolist())
            data_dict_list.append(stats_df.drop(columns=['motif']).to_dict('records'))
            explode_columns = ['md_score', 'h', 'H', 'total_motif_hits']
            
        else:
            print('Running intersecctions on differential MD score files.')
            file_columns = ['motif', 'pval', 'total_motif_hits', 'ctrl_hits', 'perturb_hits', 'ctrl_md', 'perturb_md', 'md_score']
            stats_df = pd.read_csv(file, sep='\t', header=None, \
                             names=file_columns)
            stats_df['H'] = (stats_df['ctrl_hits'] + stats_df['perturb_hits']) / 2
            stats_df['motif_key'] = stats_df['motif'].str.split('.').str[0]            
            if args.significant:
                stats_df_sig =  stats_df[(stats_df['pval'] <= PVAL_THRESH) & (stats_df['H'] >= 70 )] # Set the p-value threshold for what we're calling as significant for each set              
            elif args.enriched &~ args.depleted:
                stats_df_sig =  stats_df[(stats_df['md_score'] > DIFF_MD_THRESH) & (stats_df['H'] >= 70 )] # Set the differential MD vlaue threshold                 
            elif args.depleted &~ args.enriched:
                stats_df_sig =  stats_df[(stats_df['md_score'] < -DIFF_MD_THRESH) & (stats_df['H'] >= 70 )] # Set the differential MD vlaue threshold                      
            elif args.enriched & args.depleted:
                stats_df_sig =  stats_df[(stats_df['md_score'] < -DIFF_MD_THRESH) | (stats_df['md_score'] > DIFF_MD_THRESH) & (stats_df['H'] >= 70 )] # Set the differential MD vlaue threshold
            else:
                raise ValueError("User must specify pval, enriched, or depleted for the intersection type. Both enriched and depleted may also be specified.")
            
            motif_lists.append(stats_df_sig['motif_key'].tolist())
            data_dict_list.append(stats_df.drop(columns=['motif']).to_dict('records'))
            explode_columns = ['pval', 'total_motif_hits', 'ctrl_hits', 'perturb_hits', 'ctrl_md', 'perturb_md', 'md_score']
            
    unique_motifs_set = set(x for l in motif_lists for x in l)
    unique_motifs = list(unique_motifs_set)
    improved_data_dict_list = [num for elem in data_dict_list for num in elem]
    common = defaultdict(list)
    for d in improved_data_dict_list:
        common[d['motif_key']].append(d)
    
    results_df = pd.DataFrame([reduce(merge, value) for value in common.values()])
    
#############################################################################################            
############################### VENN DIAGRAMS ################################################
#############################################################################################    
                
########### Parse each file list into a variable for plotting in venn2/venn3

    df = pd.DataFrame(motif_lists, index=args.plot_labels).transpose()
    df.to_csv('%s/%s_motifs.txt' % (args.output, args.rootname), sep='\t', index=None)
        
    if len(args.stats_files) == 2:
        if args.colors:
            COLORS = (adjust_lightness(args.colors[0], 0.8), adjust_lightness(args.colors[1], 0.8))
        else:
            COLORS = (adjust_lightness('#154360', 0.8), adjust_lightness('#C70039', 0.8))
        A = set(motif_lists[0])
        B = set(motif_lists[1])
        venn_list = [A,B]

        if args.unweighted:
            v = venn2_unweighted(venn_list, args.plot_labels, set_colors=COLORS)
            plt.tight_layout()
            plt.savefig('%s/%s_venn2_unweighted.png' % (args.output, args.rootname), dpi=600)  
        else:
            v = venn2(venn_list, args.plot_labels, set_colors=COLORS)
            plt.tight_layout()
            plt.savefig('%s/%s_venn2.png' % (args.output, args.rootname), dpi=600)

########### Venn3 settings if number of files == 3
        
    elif len(args.stats_files) == 3:
        if args.colors:
            COLORS = (adjust_lightness(args.colors[0], 0.8), adjust_lightness(args.colors[1], 0.8), adjust_lightness(args.colors[2], 0.8))
        else:        
            COLORS = (adjust_lightness('#154360', 1.2), adjust_lightness('#FFC300', 1.2), adjust_lightness('#a60711', 1.2))
        A = set(motif_lists[0])
        B = set(motif_lists[1])
        C = set(motif_lists[2])
        venn_list = [A,B,C]
        
        if args.unweighted:
            v2 = venn3_unweighted(venn_list, args.plot_labels, set_colors=COLORS)
            plt.tight_layout()
            plt.savefig('%s/%s_venn3_unweighted.png' % (args.output, args.rootname), dpi=600)
        else:
            v = venn3(venn_list, args.plot_labels, set_colors=COLORS)
            plt.tight_layout()
            plt.savefig('%s/%s_venn3.png' % (args.output, args.rootname), dpi=600)
            
#############################################################################################            
############################### USPET PLOT ##################################################
#############################################################################################

########### Set plot lists for each stats file if <= 12 files provided
    
    elif (len(args.stats_files) > 3) & (len(args.stats_files) <= 12):
        a, b, c, d, e, f, g, h, i, j, k, l = ([] for i in range(12))       
        if len(args.stats_files) == 4:
            lists = [a,b,c,d]
        elif len(args.stats_files) == 5:
            lists = [a,b,c,d,e]
        elif len(args.stats_files) == 6:
            lists = [a,b,c,d,e,f]
        elif len(args.stats_files) == 7:
            lists = [a,b,c,d,e,f,g]
        elif len(args.stats_files) == 8:
            lists = [a,b,c,d,e,f,g,h]                    
        elif len(args.stats_files) == 9:
            lists = [a,b,c,d,e,f,g,h,i]
        elif len(args.stats_files) == 10:
            lists = [a,b,c,d,e,f,g,h,i,j]
        elif len(args.stats_files) == 11:
            lists = [a,b,c,d,e,f,g,h,i,j,k]
        elif len(args.stats_files) == 12:
            lists = [a,b,c,d,e,f,g,h,i,j,k,l]            
                 
        for bin, motif_list in zip(lists, motif_lists):
            for motif in unique_motifs:
                if motif in motif_list:
                    bin.append('True')
                else:
                    bin.append('False')
                
        data_dictionary = {'motif_key': unique_motifs, **{args.plot_labels[i]: lists[i] for i in range(len(args.plot_labels))}}  
                
        if args.colors:
            COLORS = [args.colors[0], args.colors[1], args.colors[2]]
        else:
            COLORS = ['#2e506e', '#6e2e50', '#506e2e']
                
        index_df = pd.DataFrame(data_dictionary)
        for col in args.plot_labels:
            index_df[col] = index_df[col].map({'False':False, 'True':True})
        
        motif_df = results_df[results_df['motif_key'].isin(unique_motifs)].drop_duplicates(subset = 'motif_key')
        #test = explode(motif_df, explode_columns, fill_value='', preserve_index=True)
        
        full_plot_df = pd.merge(motif_df, index_df,
                      on=['motif_key']).set_index(args.plot_labels)
        full_plot_df.to_csv('%s/%s_upset_data.txt' % (args.output, args.rootname), sep='\t')
        full_plot_df['Motif Hits'] = np.log2([np.array(x).mean() for x in full_plot_df['total_motif_hits'].values])        
        full_plot_df['Motif Hits (3kb)'] = np.log2([np.array(x).mean() for x in full_plot_df['H'].values])
        full_plot_df['MD Score'] = [np.array(x).mean() for x in full_plot_df['md_score'].values]
        upset = UpSet(full_plot_df, subset_size='count', intersection_plot_elements=3, sort_by='cardinality', show_counts=True)
        upset.add_catplot(value='Motif Hits', kind='strip', color=COLORS[0])
        upset.add_catplot(value='Motif Hits (3kb)', kind='strip', color=COLORS[1])
        upset.add_catplot(value='MD Score', kind='strip', color=COLORS[2])              
        upset.plot()
        plt.savefig('%s/%s_upset.png' % (args.output, args.rootname), dpi=600)
        
    else:
        results_df.to_csv('%s/%s_common_motif_data.txt' % (args.output, args.rootname), sep='\t')
        print('More than 12 files were provided. No plots will be generated, but a common motif dataframe will be saved.')
        
    print('Done! ' + str(datetime.datetime.now()))
    
    ### If we ever decide to add motif names inside of the venn diagrams
    
    #v.get_label_by_id('100').set_text('\n'.join(A-B-C))
    #v.get_label_by_id('110').set_text('\n'.join(A&B-C))
    #v.get_label_by_id('011').set_text('\n'.join(B&C-A))
    #v.get_label_by_id('001').set_text('\n'.join(C-A-B))
    #v.get_label_by_id('010').set_text('')
    #plt.annotate(',\n'.join(B-A-C), xy=v.get_label_by_id('010').get_position() +index_col=None
    #             np.array([0, 0.2]), xytext=(-20,40), ha='center',
    #             textcoords='offset points', 
    #             bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
    #             arrowprops=dict(arrowstyle='->',              
    #                             connectionstyle='arc',color='gray'))
    #for text in v.set_labels:
    #    text.set_fontsize(6)
    #for text in v.subset_labels:
    #    text.set_fontsize(4)    
    #plt.rcParams["figure.figsize"] = [20,12]
    #plt.tight_layout()
    #plt.savefig('venn_test.png', dpi=600)        

if __name__ == '__main__':
    main()