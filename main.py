#!/bin/usr/env python3

from phylogenetic_pipeline import cluster_sequences, calculate_gene_trees_ml, get_protein_id_to_sequence, \
    calculate_constenus_tree, calculate_supertree, select_best_supertree, visualize_newick, bootstrap_supertree
import argparse
from pathlib import Path




parser = argparse.ArgumentParser(description='Phylogenetic pipeline')
parser.add_argument('-cluster', type=Path,
                    help='Provide path to file with accessions of genomes to be clustered')
parser.add_argument('-approach', choices={'one_to_one', 'allow_paralogs', 'bootstrap'},
                    help='Select approach')
parser.add_argument('-consensus_tree', nargs=2,
                    help='Calculate consensus tree, provide path to dir containing PhyML ML trees and name for file')
parser.add_argument('-supertree', nargs=1,
                    help='Calculate supertree, provide path to dir containing PhyML ML trees')
parser.add_argument('-supertree_bootstrap', nargs=2,
                    help='Calculate supertree for bootstrap trees, provide path to dir containing PhyML ML trees and bootstrap cutoff (percent)')
parser.add_argument('-select_best_supertree', nargs=3,
                    help='Select best supertree, provide path to Fasturec output file and name to the file, index of selected tree')
parser.add_argument('-visualize_newick', type=Path,
                    help='Visualize newick tree, provide path to file containing tree in newick format')


if __name__ == '__main__':
    ncbi_dir = './data/ncbi_genomes/'
    db_dir = './data/db'
    clusters_dir = './data/clusters'
    consensus_dir = './data/consensus'

    args = parser.parse_args()
    if args.cluster:
        accessions_to_names_dict, accession_to_protein_id_list_dict, protein_id_to_sequence_dict = cluster_sequences(args.cluster)

    if args.approach == 'one_to_one':
        protein_id_to_sequence_dict = get_protein_id_to_sequence(ncbi_dir)
        calculate_gene_trees_ml(0, protein_id_to_sequence_dict, args)
    if args.approach == 'allow_paralogs':
        protein_id_to_sequence_dict = get_protein_id_to_sequence(ncbi_dir)
        calculate_gene_trees_ml(0, protein_id_to_sequence_dict, args)
    if args.approach == 'bootstrap':
        protein_id_to_sequence_dict = get_protein_id_to_sequence(ncbi_dir)
        calculate_gene_trees_ml(100, protein_id_to_sequence_dict, args)

    if args.consensus_tree:
        calculate_constenus_tree(args.consensus_tree[0], args.consensus_tree[1])
    if args.supertree:
        calculate_supertree(args.supertree[0])
    if args.supertree_bootstrap:
        bootstrap_supertree(args.supertree_bootstrap[0], args.supertree_bootstrap[1])

    if args.select_best_supertree:
        select_best_supertree(args.select_best_supertree[0], args.select_best_supertree[1], int(args.select_best_supertree[2]))
    if args.visualize_newick:
        visualize_newick(args.visualize_newick)







