"""
Zacząć od publikacji z drzewem gatunków np.
Z blast interesuje nas e-value dla każdej pary sekwencji
MCL buduje na podstawie tego graf, przeprowadza oblicenia matematyczne i klastruje
zapisać gdzieś nazwa-gen, kiedy zastępuję

W zależności od narzędzia:
miękkie klastrowanie - jeden obiekt może trafić do jednego lub więcej klastrów
twarde klastrowanie - jeden obiekt, jeden klaster
Ile mamy klastró i kiedy stwierdzamy, ze jany klaster jest śmieciem, np drzewa o małym rozmiarze, kiedy wiele gatunków
Na pewno usuwamy małe klastry

Dla każdej grupy multialignment, w obrębie klastra

Najlepiej ML - PhyML
NJ - najgorsze drzewa, najmniej prawdopodobne, dziwactwa w postaci ujemnych wartości krawędzi, wtedy należy usunąć

Drzewo dla klastra - drzewa genów

6. Ze wszystkich drzew jedno drzewo - drzewo gatunków
Consensus lub superdrzewo - o tym jeszcze zostanie powiedziane

7. Porównać z drzewem NCBI, RF distance

Dodatkowe pkt.
Filtr: Usunąć drzewa o kiepskim wsparciu, żeby nie wprowadzały szumu i nie zniekształcały drzewa gatunków

https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001342
"""

# Program będzie napisany tak, żeby wywoływać z command line kolejne eatpy

import os
import pandas as pd
from Bio import Entrez, SeqIO, AlignIO

# Read in genome accessions and shortcuts for organisms names from file as dict
def parse_genomes_ids(path):
    path = os.path.relpath(path)
    accessions_to_names_dict = {}

    with open(path, 'r') as f:
        for line in f:
            words = line.strip().split()
            acc = words[0]
            name = words[1]
            accessions_to_names_dict[acc] = name
    return(accessions_to_names_dict)


def make_empty_folder(dir_name):
    try:
        os.makedirs(dir_name)
        print(f"Directory {dir_name} created")
    except FileExistsError:
        print(f"Directory {dir_name} already exists")
        # jeśli chcesz mogę usunąć i stworzyć kolejne
        # for fileName in os.listdir(dirName):
        #     os.remove(os.path.join(dirName, fileName))


# Download genomes in Gene Bank format by genome accessionss and save in dedicated directory
def get_ncbi_genomes(genome_accessions, ncbi_dir):
    make_empty_folder(ncbi_dir)

    for genome_id in genome_accessions.keys(): # iterate over genome accessions in dict
        Entrez.email = "i.gozdziewska@student.uw.edu.pl"
        filename = f'data/ncbi_genomes/{genome_id}.gb'
        print(f'Downloading: {filename}')
        handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")

        with open(filename, 'w') as f:
            f.write(handle.read())


# Convert gb format to fasta, keep protein_id as name of gene
def gb_to_fasta(ncbi_dir):
    for filename in os.listdir(ncbi_dir):
        seq_record = os.path.splitext(filename)[0]
        gb_file = os.path.join(ncbi_dir, filename)
        faa_filename = os.path.join(ncbi_dir, f'{seq_record}.faa')
        input_handle = open(gb_file, "r")
        output_handle = open(faa_filename, "w")

        for seq_record in SeqIO.parse(input_handle, "genbank"):
            print("Dealing with GenBank record {}".format(seq_record.id))
            for seq_feature in seq_record.features:
                if seq_feature.type == "CDS":
                    if 'translation' in seq_feature.qualifiers:
                        CDS_seq = seq_feature.qualifiers['translation'][0]
                        protid = seq_feature.qualifiers['protein_id'][0]
                        output_handle.write(">" + protid + "\n" + str(CDS_seq) + "\n")
        output_handle.close()
        input_handle.close()


# Change protein name TODO
def add_organism_name(accessions_to_names_dict, ncbi_dir):
    for filename in os.listdir(ncbi_dir):
        if filename.endswith(".faa"):
            ncbi_id = os.path.splitext(filename)[0] # without extension
            old_file = os.path.join(ncbi_dir, filename)
            with open(old_file, 'r+') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith('>'):
                        words = line.split()
                        protein_id = words[0].split('>')[1]
                        lines[i] = f'>{protein_id}_{accessions_to_names_dict[ncbi_id]}\n'
                f.seek(0)
                f.truncate()
                for line in lines:
                    f.write(line)


def conacatenate_fasta_files(ncbi_dir):
    outfile_path = f'./{ncbi_dir}concatenated_fasta_files.faa'
    with open(outfile_path, 'w') as outfile:
        for filename in os.listdir(ncbi_dir):
            if filename.endswith(".faa"):
                filename_path = f'{ncbi_dir}/{filename}'
                with open(filename_path) as infile:
                    for line in infile:
                        outfile.write(line)


def get_protein_id_to_sequence(ncbi_dir):
    outfile_path = f'./{ncbi_dir}concatenated_fasta_files.faa'

    protein_id_to_sequence_dict = {}
    for record in SeqIO.parse(outfile_path, "fasta"):
        protein_id_to_sequence_dict[record.id] = str(record.seq)

    return protein_id_to_sequence_dict


# TODO tuning parametrów klastrowania
def clustering(ncbi_dir, db_dir, clusters_dir):
    outfile_path = f'./{ncbi_dir}concatenated_fasta_files.faa'

    make_empty_folder(db_dir)
    make_empty_folder(clusters_dir)

    create_db = os.system(f'mmseqs createdb {outfile_path} {db_dir}/db')
    execute_clustering = os.system(f'mmseqs cluster {db_dir}/db {clusters_dir}/db_clu {clusters_dir}/tmp')
    output_to_tsv = os.system(f'mmseqs createtsv {db_dir}/db {db_dir}/db {clusters_dir}/db_clu {clusters_dir}/DB_clu.tsv')
    output_to_fasta = os.system(f'mmseqs result2flat {db_dir}/db {db_dir}/db {clusters_dir}/db_clu {clusters_dir}/DB_clu_seq.fasta')


def file_for_cluster(clusters_dir):
    db_clu_tsv_path = f'{clusters_dir}/DB_clu.tsv'
    db_clu_tsv = pd.read_csv(db_clu_tsv_path, sep='\t', header=None, names=['representative', 'elements'])
    print("Number of clusters: ", len(db_clu_tsv.representative.unique()))
    print("Number of proteins clustered: ", len(db_clu_tsv.elements.unique()))

    make_empty_folder(f'{clusters_dir}/clusters')
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
            with open(file_path, 'w') as f:
                f.write(f'>{row[1]}') # row[0] will always be the first in cluster as row[1]

        previous = row[0]


def filter_clusters(clusters_dir, cutoff=10):
    cluster_files_path = f'{clusters_dir}/cluster'

    for fn in os.listdir(cluster_files_path):
        if fn.endswith('.cluster'):
            num_lines = sum(1 for line in open(f'{cluster_files_path}/{fn}'))
            if num_lines < cutoff:
                os.remove(f'{cluster_files_path}/{fn}')

    organisms_list = []
    for rm_fn in os.listdir(cluster_files_path): # filter remaining
        if rm_fn.endswith('.cluster'):
            with open(f'{cluster_files_path}/{rm_fn}', 'r') as rf:
                for line in rf:
                    line = line.rstrip()
                    organism = line.split('_')[1]
                    organisms_list.append(organism)
                if len(set(organisms_list)) == 1:
                    os.remove(f'{cluster_files_path}/{rf}')



def cluster_to_fasta(clusters_dir, protein_id_to_sequence_dict):
    cluster_files_path = f'{clusters_dir}/cluster'

    fasta_path = f'{clusters_dir}/fasta'
    make_empty_folder(fasta_path)

    for cluster_fn in os.listdir(cluster_files_path):
        cluster_number = os.path.splitext(cluster_fn)[0]
        faa_name = f'{cluster_number}.faa'
        with open(f'{cluster_files_path}/{cluster_fn}', 'r') as cluster_f:
            with open(f'{fasta_path}/{faa_name}', 'a') as fasta_f:
                for line in cluster_f:
                    line = line.rstrip()
                    protein_id = line.split('>')[1]
                    # organism_name = protein_id.split('_')[1]
                    protein_seq = protein_id_to_sequence_dict[protein_id]
                    fasta_f.write(f'>{protein_id}\n{protein_seq}\n')


# TODO tuning parametrów aln
def calculate_multialignment_for_cluster(clusters_dir):
    fasta_files_path = f'{clusters_dir}/fasta'
    clustalw2_path = './tools/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2'

    aln_path = f'{clusters_dir}/aln'
    make_empty_folder(aln_path)

    for fasta_fn in os.listdir(fasta_files_path):
        fasta_fn_path = f'{fasta_files_path}/{fasta_fn}'

        cluster_number = os.path.splitext(fasta_fn)[0]
        aln_fn = f'{cluster_number}.aln'
        aln_fn_path = f'{aln_path}/{aln_fn}'
        print(aln_fn_path)
        mln_command = os.system(f'{clustalw2_path} -infile={fasta_fn_path} -outfile={aln_fn_path}')

# Used seqret to convert hem_sub_alpha1_aln.fasta to seqret_interleaved_out.phylip
# Used PhyML 3.0 to calculate ML tree, output: seqret_interleaved_out_phylip_phyml



if __name__ == '__main__':
    ncbi_dir = './data/ncbi_genomes'
    db_dir = './data/db'
    clusters_dir = './data/clusters'
    accessions_to_names_dict = parse_genomes_ids('genomes.txt')
    # get_ncbi_genomes(accessions_to_names_dict, ncbi_dir)
    # gb_to_fasta(ncbi_dir)
    # add_organism_name(accessions_to_names_dict, ncbi_dir)
    # conacatenate_fasta_files(ncbi_dir)
    # protein_id_to_sequence_dict = get_protein_id_to_sequence(ncbi_dir)
    #
    # clustering(ncbi_dir, db_dir, clusters_dir)
    # file_for_cluster(clusters_dir)
    # filter_clusters(clusters_dir)
    # cluster_to_fasta(clusters_dir, protein_id_to_sequence_dict)
    #
    # calculate_multialignment_for_cluster(clusters_dir)

# TODO change names to only organism
# TODO ML trees
# TODO consensus tree

# TODO trees without paralogs vs one-to_one



