from __future__ import print_function
from itertools import combinations
import networkx as nx
import pandas as pd
import argparse
from os import path

#
# Usage example:
# 
# $ python3 tf_change_explanations.py -p 0.01 -d some_directory/Untreated_vs_Treatment_differential_md_scores.txt -o untreated_vs_treated_explanation.txt
#

basedir = path.split(__file__)[0]
KNOWLEDGE_GRAPH = "%s/public_knowledge/prot_reactome_interactions.pkl" % basedir
TF_UNIPROT_MAP = "%s/public_knowledge/human_TFs_to_uniprot.txt" % basedir
TF_HOCOMOCO_UNIPROT_MAP = "%s/public_knowledge/HOCOMOCOv11_to_uniprot.txt" % basedir
ONTO_LABELS = "%s/public_knowledge/all_labels.tsv" % basedir

PATH_SEARCH_DEPTH = 2

EXTRA_NODES = []
EXTRA_NODES_LABEL = []
UNINTERESTING_RELATIONS = []


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--p-value', dest='p_val', help="P-value cutoff to determine which TFs to include from DAStk's output (default=0.05).", required=False, default=0.05)
    parser.add_argument('-d', '--dastk-results', dest='dastk_results', help="Results file from DAStk (*differential_md_scores.txt) used to find relations between the most significant TF changes in activity.", required=True)
    parser.add_argument('-o', '--output', dest='output_filename', help='Output filename for the report.', required=True)
    parser.add_argument('-u', '--uninteresting-nodes', dest='uninteresting_nodes', help="File listing ontology concept URIs to ignore during the pathway searches, because they are uninformative for TFs (e.g.\"binds to DNA\") or they don't apply to the current study, just to minimize noise. One URI per line, can optionally add a description after a TAB to track the label of the ignored intersecting concepts. See the example file for a format guide.", required=False)
    parser.add_argument('-e', '--extra-concepts', dest='extra_concepts', help="File listing extra ontology concepts to include in the pathway searches, that are relevant to this study. This is a two-column (TAB-separated) list, first the ontology URI and second a label you'd like to use in the report.", required=False)
    args = parser.parse_args()

    G = nx.read_gpickle(KNOWLEDGE_GRAPH)

    if args.extra_concepts:
        with open(args.extra_concepts, 'r') as fd:
            for line in fd.readlines():
                if len(line) > 0:
                    chunks = line.split('\t')
                    EXTRA_NODES.append(chunks[0])
                    EXTRA_NODES_LABEL.append(chunks[1][:-1])

    if args.uninteresting_nodes:
        with open(args.uninteresting_nodes, 'r') as fd:
            for line in fd.readlines():
                if len(line) > 0:
                    chunks = line.split('\t')
                    UNINTERESTING_RELATIONS.append(chunks[0])

    # cache the Uniprot mappings
    uniprot_tf_ids = dict()
    for mapping_file in [TF_UNIPROT_MAP, TF_HOCOMOCO_UNIPROT_MAP]:
        with open(mapping_file, 'r') as fd:
            for line in fd.readlines():
                chunks = line.split('\t')
                uniprot_tf_ids[chunks[0]] = chunks[1]

    # Gather the list of significantly changing TFs
    significant_tfs = []
    significant_tfs_key = []
    for extra_node in EXTRA_NODES:
        significant_tfs.append(extra_node)
    for extra_key in EXTRA_NODES_LABEL:
        significant_tfs_key.append(extra_key)
    with open(args.dastk_results, 'r') as fd:
        for line in fd.readlines():
            chunks = line.split('\t')
            tf_name = chunks[0].split('_')[0]
            p_value = float(chunks[1])
            if p_value < float(args.p_val) and tf_name not in significant_tfs_key:
                significant_tfs.append("<http://purl.obolibrary.org/obo/PR_%s>" % uniprot_tf_ids[tf_name])
                significant_tfs_key.append(tf_name)

    labels_df = pd.read_csv(ONTO_LABELS, sep="\t", na_filter=False, names=['concept', 'label'])
    idx = 0
    for extra_concept in EXTRA_NODES:
        labels_df.loc[labels_df.index.max() + 1] = [extra_concept, EXTRA_NODES_LABEL[idx]]
        idx += 1

    output_fd = open(args.output_filename, 'w')
    output_fd.write("Transcription factors displaying a significant difference in activity:\n%s" % ", ".join(sorted(significant_tfs_key)))
    #print(sorted(significant_tfs_key))

    common_intersections = dict()
    rel_list = dict()
    rel_list[1000] = []

    for pair in combinations(significant_tfs, 2):
        #print(pair)
        try:
            for path in nx.all_simple_paths(G, source=pair[0], target=pair[1], cutoff=PATH_SEARCH_DEPTH):
                idx_1 = significant_tfs.index(pair[0])
                idx_2 = significant_tfs.index(pair[1])
                if len(path) == 2:
                    relation_concept = G.edges[pair[0], pair[1]]['label']
                    relation = labels_df[(labels_df.concept==relation_concept)].label.values[0]
                    if 0 in rel_list:
                        rel_list[0].append("%s %s %s\n" % (significant_tfs_key[idx_1], relation, significant_tfs_key[idx_2]))
                    else:
                        rel_list[0] = ["%s %s %s\n" % (significant_tfs_key[idx_1], relation, significant_tfs_key[idx_2])]
                elif len(path) == 3:
                    if path[1] in UNINTERESTING_RELATIONS:
                        continue
                    relation_1 = G.edges[path[0], path[1]]['label']
                    relation_1_label = labels_df[(labels_df.concept==relation_1)].label.values[0]
                    relation_2 = G.edges[path[1], path[2]]['label']
                    relation_2_label = labels_df[(labels_df.concept==relation_2)].label.values[0]
                    intersecting_concept = labels_df[(labels_df.concept==path[1])].label.values[0]
                    if intersecting_concept in common_intersections:
                        if relation_1_label in common_intersections[intersecting_concept]:
                            common_intersections[intersecting_concept][relation_1_label].append(significant_tfs_key[idx_1])
                        else:
                            common_intersections[intersecting_concept][relation_1_label] = [significant_tfs_key[idx_1]]
                        if relation_2_label in common_intersections[intersecting_concept]:
                            common_intersections[intersecting_concept][relation_2_label].append(significant_tfs_key[idx_2])
                        else:
                            common_intersections[intersecting_concept][relation_2_label] = [significant_tfs_key[idx_2]]
                    else:
                        common_intersections[intersecting_concept] = dict()
                        common_intersections[intersecting_concept][relation_1_label] = [significant_tfs_key[idx_1]]
                        if relation_1 == relation_2:
                            common_intersections[intersecting_concept][relation_1_label].append(significant_tfs_key[idx_2])
                        else:
                            common_intersections[intersecting_concept][relation_2_label] = [significant_tfs_key[idx_2]]
                else:
                    if any([x in path for x in UNINTERESTING_RELATIONS]):
                        continue
                    last_node = None
                    path_desc = ''
                    for i in range(len(path) - 1):
                        if last_node:
                            node1 = last_node
                            path_desc += ", and "
                        else:
                            node1 = labels_df[(labels_df.concept==path[i])].label.values[0]
                        node2 = labels_df[(labels_df.concept==path[i+1])].label.values[0]
                        rel_concept = G.edges[path[i], path[i+1]]['label']
                        rel = labels_df[(labels_df.concept==rel_concept)].label.values[0]
                        if i == len(path) - 2:
                            path_desc += "%s %s %s" % (node2, rel, node1)
                        else:
                            path_desc += "%s %s %s" % (node1, rel, node2)
                        last_node = node2
                    rel_list[1000].append(path_desc + '\n')
        except nx.exception.NetworkXNoPath as e:
            #print("No paths")
            pass
        except nx.exception.NodeNotFound as e:
            #print("\nWarning: No node for %s or %s" % (pair[0], pair[1]))
            pass
        except IndexError as e:
            print("\nWarning: Problem adding this path (likely missing label, or deprecated ontology concept):")
            print(path)

    for intersection in common_intersections.keys():
        for relation in common_intersections[intersection].keys():
            tfs_involved = sorted(list(set(common_intersections[intersection][relation])))
            if len(tfs_involved) == 1:
                explanation = "%s %s %s\n\n" % (tfs_involved[0], relation, intersection)
            elif len(tfs_involved) == 2:
                explanation = "%s and %s: %s %s\n\n" % (", ".join(tfs_involved[:-1]), tfs_involved[-1], relation, intersection)
            else:
                explanation = "%s, and %s: %s %s\n\n" % (", ".join(tfs_involved[:-1]), tfs_involved[-1], relation, intersection)
            if len(tfs_involved) in rel_list:
                rel_list[len(tfs_involved)].append(explanation)
            else:
                rel_list[len(tfs_involved)] = [explanation]


    # Output results
    output_fd.write("\n\nHere's what we know about these TFs presenting significant activity changes (p=%.2E):\n-----------------------------------------------------------------------------------------\n" % float(args.p_val))
    if 0 in rel_list:
        output_fd.write("\nDirect interactions between each of these TFs:\n--------------------------------\n")
        for line in rel_list[0]:
            output_fd.write(line)

    output_fd.write("\nOther ways these TFs are related:\n---------------------------------\n")
    for tf_count in sorted(rel_list.keys(), reverse=True):
        if tf_count > 0:
            for line in rel_list[tf_count]:
                output_fd.write(line)

    output_fd.close()


if __name__=='__main__':
    main()    

