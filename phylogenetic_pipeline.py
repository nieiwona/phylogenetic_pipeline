from Bio import Entrez, SeqIO, Phylo, AlignIO
import gzip
import os
import shutil
import urllib.request as request
from contextlib import closing
import pandas as pd
import logging
import sys
from statistics import median
from shutil import copyfile
from collections import Counter
import time
import matplotlib.pyplot as plt
from io import StringIO
from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus, strict_consensus
import re
from itertools import permutations


ncbi_dir = './data/ncbi_genomes/'
db_dir = './data/db'
clusters_dir = './data/clusters'
consensus_dir = './data/consensus'


def setup_logging(arg):
    logging.basicConfig(filename='pipeline.log',
                        filemode='a',
                        format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                        datefmt='%H:%M:%S',
                        level=logging.INFO)
    logging.info(f'Start with argument {arg}')


# Stage 1
def parse_genomes_ids(path):
    path = os.path.relpath(path)
    accessions_to_names_dict = {}
    with open(path, 'r') as f:
        for line in f:
            words = line.strip().split()
            acc = words[0]
            name = words[1]
            accessions_to_names_dict[acc] = name
    logging.info(f"Path to accessions file {path}")
    logging.info(f"Accessions of genomes to be analysed and organisms names: {accessions_to_names_dict}")
    return accessions_to_names_dict


def make_empty_folder(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
        print(f"Directory {dir_name} created")
        logging.info(f"Directory {dir_name} created")


def download_file_from_ftp(remote_path, local_path):
    with closing(request.urlopen(remote_path)) as remote_file:
        with open(local_path, "wb") as local_file:
            shutil.copyfileobj(remote_file, local_file)


def unzip_file(input_file_path, output_file_path):
    with gzip.open(input_file_path, "rb") as input_file:
        with open(output_file_path, "wb") as output_file:
            shutil.copyfileobj(input_file, output_file)


def get_ncbi_genomes(base_path, accessions_to_names_dict):
    start_time = time.time()
    make_empty_folder(base_path)
    Entrez.email = "iwonaa.gozdziewska@gamil.com"
    for assembly in accessions_to_names_dict.keys():
        output_file_name = assembly + ".fna"
        if os.path.isfile(base_path + output_file_name):
            print(f"Assembly {assembly} is already downloaded")
            logging.info(f"Assembly {assembly} is already downloaded")
            continue
        try:
            search_handle = Entrez.esearch(db="assembly", term=assembly)
            assembly_id = Entrez.read(search_handle)["IdList"][0]
        except IndexError:
            print(f"Entrez.esearch did not find assembly {assembly}, ignoring")
            logging.info(f"Entrez.esearch did not find assembly {assembly}, ignoring")
            continue
        summary_handle = Entrez.esummary(db="assembly", id=assembly_id, report="full")
        assembly_info = Entrez.read(summary_handle)
        remote_folder_path = assembly_info["DocumentSummarySet"]["DocumentSummary"][0]["FtpPath_RefSeq"] # FtpPath_GeneBank
        print(remote_folder_path)
        remote_file_name = remote_folder_path.rsplit("/", 1)[-1] + "_protein.faa.gz"
        print(remote_file_name)
        remote_file_path = remote_folder_path + "/" + remote_file_name
        print(remote_file_path)
        print(f"Downloading: {remote_file_path}")
        logging.info(f"Downloading: {remote_file_path}")
        download_file_from_ftp(remote_file_path, base_path + remote_file_name)
        unzip_file(base_path + remote_file_name, base_path + output_file_name)
        os.remove(base_path + remote_file_name)
    logging.info(f'Downloading completed. Execution time {round((time.time() - start_time) / 60, 2)} min')


def get_accession_to_protein_id_list(ncbi_dir):
    accession_to_protein_id_list_dict = {}
    for fn in os.listdir(ncbi_dir):
        if fn.endswith('.fna'):
            accession = os.path.splitext(fn)[0]
            protein_id_lst = set()
            fn_path = f'{ncbi_dir}/{fn}'
            with open(fn_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        protein_id = line.split()[0]
                        protein_id = protein_id.split('>')[1]
                        protein_id_lst.add(protein_id)
            logging.info(f"Genome {accession} include {len(protein_id_lst)} proteins")
            accession_to_protein_id_list_dict[accession] = protein_id_lst
    return accession_to_protein_id_list_dict


def conacatenate_fasta_files(ncbi_dir):
    outfile_path = f'{ncbi_dir}/concatenated_fasta_files.faa'
    with open(outfile_path, 'w') as outfile:
        for filename in os.listdir(ncbi_dir):
            if filename.endswith('.fna'):
                filename_path = f'{ncbi_dir}/{filename}'
                with open(filename_path) as infile:
                    for line in infile:
                        outfile.write(line)
    logging.info(f"File {ncbi_dir}/concatenated_fasta_files.faa created")


def get_protein_id_to_sequence(ncbi_dir):
    outfile_path = f'./{ncbi_dir}concatenated_fasta_files.faa'

    protein_id_to_sequence_dict = {}
    for record in SeqIO.parse(outfile_path, "fasta"):
        protein_id_to_sequence_dict[record.id] = str(record.seq)

    return protein_id_to_sequence_dict


def clustering(ncbi_dir, db_dir, clusters_dir):
    outfile_path = f'./{ncbi_dir}concatenated_fasta_files.faa'

    make_empty_folder(db_dir)
    make_empty_folder(clusters_dir)
    start_time = time.time()
    print('Clustering...')
    logging.info('Clustering with parameters -min-seq-id 0.5 -c 0.8 --cov-mode 1')
    create_db = os.system(f'mmseqs createdb {outfile_path} {db_dir}/db')
    execute_clustering = os.system(f'mmseqs cluster {db_dir}/db {clusters_dir}/db_clu {clusters_dir}/tmp -c 0.8 --cov-mode 0 --min-seq-id 0.5') # default
    output_to_tsv = os.system(f'mmseqs createtsv {db_dir}/db {db_dir}/db {clusters_dir}/db_clu {clusters_dir}/DB_clu.tsv')
    logging.info(f'Clustering completed. Execution time {round((time.time() - start_time) / 60, 2)} min')
    print('Clustering completed')


def file_for_cluster(clusters_dir):
    logging.info('Writting clusters to .cluster files')
    db_clu_tsv_path = f'{clusters_dir}/DB_clu.tsv'
    db_clu_tsv = pd.read_csv(db_clu_tsv_path, sep='\t', header=None, names=['representative', 'elements'])
    logging.info(f"Number of clusters: {len(db_clu_tsv.representative.unique())}")
    logging.info(f"Number of proteins clustered: {len(db_clu_tsv.elements.unique())}")

    make_empty_folder(f'{clusters_dir}/cluster')
    previous = None
    i = 0

    for index, row in db_clu_tsv.iterrows():
        if row[0] == previous:
            fn = f'{i}.cluster'
            file_path = f'{clusters_dir}/cluster/{fn}'
            with open(file_path, 'a') as f:
                f.write(f'\n>{row[1]}')
        if row[0] != previous:
            i += 1
            fn = f'{i}.cluster'
            file_path = f'{clusters_dir}/cluster/{fn}'
            with open(file_path, 'w') as nf:
                nf.write(f'>{row[1]}') # row[0] will always be the first in cluster as row[1]

                previous = row[0]
    logging.info("Writting clusters to .cluster files completed")


def add_accession_to_protein_id_in_cluster(clusters_dir, accession_to_protein_id_list_dict, accessions_to_names_dict):
    cluster_files_path = f'{clusters_dir}/cluster'
    for filename in os.listdir(cluster_files_path):
        f = open(f'{cluster_files_path}/{filename}', 'r+')
        protein_ids_with_acccessions = []
        for line in f:
            line = line.rstrip()
            protein_id = line.split('>')[1]
            for acc in accession_to_protein_id_list_dict.keys():
                if protein_id in accession_to_protein_id_list_dict[acc]:
                    protein_ids_with_acccessions.append(f'{accessions_to_names_dict[acc]} {protein_id}')
        f.seek(0)
        for el in protein_ids_with_acccessions:
            f.write(f'>{el}\n')
        f.truncate()
        f.close()


def filter_clusters(clusters_dir):
    cluster_files_path = f'{clusters_dir}/cluster'
    allow_paralogs_path = f'{clusters_dir}/cluster_allow_paralogs'
    one_to_one_path = f'{clusters_dir}/cluster_one_to_one'
    bootstrap_path = f'{clusters_dir}/cluster_bootstrap'
    make_empty_folder(allow_paralogs_path)
    make_empty_folder(one_to_one_path)
    make_empty_folder(bootstrap_path)

    stats = []
    for rm_fn in os.listdir(cluster_files_path): # filter remaining
        if rm_fn.endswith('.cluster'):
            with open(f'{cluster_files_path}/{rm_fn}', 'r') as rf:
                accessions_lst = []
                for line in rf:
                    line = line.rstrip()
                    line_els = line.replace('>', '').split(' ')
                    org = line_els[0]
                    accessions_lst.append(org)
                acc_counts = Counter(accessions_lst)
                stats.append(len(accessions_lst))
                # Allow paralogs
                if len(set(accessions_lst)) == 10:
                    if all(i <= 2 for i in list(acc_counts.values())): # up to 2 paralogs
                        copyfile(f'{cluster_files_path}/{rm_fn}', f'{allow_paralogs_path}/{rm_fn}')
                        # One to one
                        if all(i == 1 for i in list(acc_counts.values())):
                            copyfile(f'{cluster_files_path}/{rm_fn}', f'{one_to_one_path}/{rm_fn}')
                            copyfile(f'{cluster_files_path}/{rm_fn}', f'{bootstrap_path}/{rm_fn}')


    stats_dict = Counter(stats)
    plt.bar(range(len(stats_dict)), list(stats_dict.values()), align='center')
    plt.xlabel('Size of clusters')
    plt.ylabel('Frequency')
    plt.title('Frequency of the size of clusters')
    plt.savefig('clusters_histogram.png')
    logging.info(f'Cluster size frequency: {stats_dict}')
    print(f'Mean size of clusters: {sum(stats) / len(stats)}')
    logging.info(f'Mean size of clusters: {sum(stats) / len(stats)}')
    print(f'Median size of clusters:  {median(stats)}')
    logging.info(f'Median size of clusters:  {median(stats)}')


def cluster_sequences(Path):
    ncbi_dir = './data/ncbi_genomes/'
    db_dir = './data/db'
    clusters_dir = './data/clusters'

    setup_logging('clustering')
    accessions_to_names_dict = parse_genomes_ids(Path)
    get_ncbi_genomes(ncbi_dir, accessions_to_names_dict)
    accession_to_protein_id_list_dict = get_accession_to_protein_id_list(ncbi_dir)
    conacatenate_fasta_files(ncbi_dir)
    protein_id_to_sequence_dict = get_protein_id_to_sequence(ncbi_dir)
    clustering(ncbi_dir, db_dir, clusters_dir)
    file_for_cluster(clusters_dir)
    add_accession_to_protein_id_in_cluster(clusters_dir, accession_to_protein_id_list_dict, accessions_to_names_dict)
    filter_clusters(clusters_dir)

    return accessions_to_names_dict, accession_to_protein_id_list_dict, protein_id_to_sequence_dict


# Stage 2
def cluster_to_fasta(cluster_dir, protein_id_to_sequence_dict):
    fasta_path = f'{cluster_dir}/fasta'
    make_empty_folder(fasta_path)

    for cluster_fn in os.listdir(cluster_dir):
        if cluster_fn.endswith('.cluster'):
            cluster_number = os.path.splitext(cluster_fn)[0]
            faa_name = f'{cluster_number}.faa'
            with open(f'{cluster_dir}/{cluster_fn}', 'r') as cluster_f:
                with open(f'{fasta_path}/{faa_name}', 'a') as fasta_f:
                    orgs = []
                    for line in cluster_f:
                        els = line.replace('>', '').rstrip().split(' ')
                        protein_id = els[1]
                        org = els[0]
                        protein_seq = protein_id_to_sequence_dict[protein_id]
                        if org not in orgs:
                            fasta_f.write(f'>{org} {protein_id}\n{protein_seq}\n')
                            orgs.append(org)
                        else:
                            fasta_f.write(f'>{org}_ {protein_id}\n{protein_seq}\n')


def calculate_multialignment_for_cluster(cluster_dir):
    logging.info('Calculating multialignment')
    start_time = time.time()
    clustalw2_path = './tools/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2'

    aln_path = f'{cluster_dir}/aln'
    make_empty_folder(aln_path)
    fasta_files_path = f'{cluster_dir}/fasta'

    for fasta_fn in os.listdir(fasta_files_path):
        fasta_fn_path = f'{fasta_files_path}/{fasta_fn}'

        cluster_number = os.path.splitext(fasta_fn)[0]
        aln_fn = f'{cluster_number}.aln'
        aln_fn_path = f'{aln_path}/{aln_fn}'
        print(aln_fn_path)
        mln_command = os.system(f'{clustalw2_path} -align -TYPE=PROTEIN -infile={fasta_fn_path} -outfile={aln_fn_path}')
    print(f'Multialignment completed. Execution time {round((time.time() - start_time) / 60, 2)} min')
    logging.info(f'Multialignment completed. Execution time {round((time.time() - start_time) / 60, 2)} min')


def aln_to_phylip(clusters_dir, ml_dir):
    aln_path = f'{clusters_dir}/aln'
    make_empty_folder(ml_dir)

    for aln_fn in os.listdir(aln_path):
        aln_fn_path = f'{aln_path}/{aln_fn}'
        with open(aln_fn_path, 'r') as file:
            filedata = file.read()
        filedata = filedata.replace('_', ' ')
        with open(aln_fn_path, 'w') as file:
            file.write(filedata)
        cluster_number = os.path.splitext(aln_fn)[0]
        cluster_fn = f'{cluster_number}.phylip'
        phylip_fn_path = f'{ml_dir}/{cluster_fn}'
        seqret_command = os.system(f'seqret -sequence aln::{aln_fn_path} -outseq phylip::{phylip_fn_path}')


def calculate_ml_trees(bootstrap, ml_dir):
    logging.info('Calculating ML trees with parameters phyml -i -d aa -m LG -b 0 -c 4 --run_id=lg')
    start_time = time.time()
    for phylip_fn in os.listdir(ml_dir):
        phylip_fn_path = f'{ml_dir}/{phylip_fn}'
        phyml_command = os.system(f'phyml -i {phylip_fn_path} -d aa -m LG -b {bootstrap} -c 4 --run_id=lg')
    logging.info(f'ML trees calculated. Execution time {round((time.time() - start_time) / 60, 2)} min')


def calculate_gene_trees_ml(bootstrap, protein_id_to_sequence_dict, args):
    setup_logging(args.approach)
    cluster_to_fasta(f'{clusters_dir}/cluster_{args.approach}', protein_id_to_sequence_dict)
    calculate_multialignment_for_cluster(f'{clusters_dir}/cluster_{args.approach}')
    aln_to_phylip(f'{clusters_dir}/cluster_{args.approach}', f'{clusters_dir}/cluster_{args.approach}/ml')
    calculate_ml_trees(bootstrap, f'{clusters_dir}/cluster_{args.approach}/ml')


# Stage 3
def calculate_constenus_tree(ml_dir, name):
    make_empty_folder(f'./data/consensus')
    newicks = {}
    for fn in os.listdir(ml_dir):
        if fn.endswith('.phylip_phyml_tree_lg.txt'):
            cluster = int(fn.split('.')[0])
            with open(f'{ml_dir}/{fn}') as newick:
                for line in newick:
                    line = line.rstrip()
                    newicks[cluster] = line


    def read_newick(treedata):
        handle = StringIO(treedata)
        return Phylo.read(handle, "newick")

    trees = [read_newick(newicks[key]) for key in newicks.keys()]
    majority_tree = majority_consensus(trees, 0.5)
    Phylo.write(majority_tree, f'./data/consensus/majority_consensus_{name}.newick', "newick")


def calculate_supertree(ml_trees_path):
    all_trees_path = f'{ml_trees_path}/all_ml_trees.txt'
    for fn in os.listdir(ml_trees_path):
        if fn.endswith('.phylip_phyml_tree_lg.txt'):
            with open(all_trees_path, 'a') as outfile:
                with open(f'{ml_trees_path}/{fn}', 'r') as infile:
                    for line in infile:
                        outfile.write(line)
    os.system(f'./tools/fasturec32bit/fasturec -G {all_trees_path} -Y -r10 -k1 -j10')


def select_best_supertree(supertrees_path, name, number):
    with open(supertrees_path, 'r') as f:
        number = number - 1
        best = f.readlines()
        best = best[number]
        score_and_tree = best.split()
        score = score_and_tree[0]
        tree = score_and_tree[1]
    with open(f'./best_supertree_{name}_score{score}.txt', 'w') as infile:
        infile.write(tree)


def visualize_newick(path):
    tree = Phylo.read(path, "newick")
    tree.ladderize()  # Flip branches so deeper clades are displayed at top
    Phylo.draw(tree)


def bootstrap_supertree(ml_trees_path, cutoff):
    cutoff = int(cutoff)
    all_trees_path = f'{ml_trees_path}/all_ml_trees.txt'
    for fn in os.listdir(ml_trees_path):
        if fn.endswith('.phylip_phyml_tree_lg.txt'):
            with open(all_trees_path, 'a') as outfile:
                with open(f'{ml_trees_path}/{fn}', 'r') as infile:
                    for line in infile:
                        bootstrap_support = [x for x in re.findall(r'\)(.*?)\:', str(line))]
                        bootstrap_support = list(map(int, bootstrap_support))
                        if all(i >= cutoff for i in bootstrap_support):
                            without_branch_support = re.sub(r'(?<=\)).+?(?=\:)', '', str(line))
                            outfile.write(without_branch_support)
    os.system(f'./tools/fasturec32bit/fasturec -G {all_trees_path} -Y -r10 -k1 -j10'