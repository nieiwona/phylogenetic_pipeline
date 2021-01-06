"""
Zacząć od publikacji z drzewem gatunków np.
Z blast interesuje nas e-value dla każdej pary sekwencji
MCL buduje na podstawie tego graf, przeprowadza oblicenia matematyczne i klastruje
zapisać gdzieś nazwa-gen, kiedy zastępuję

W zależności od narzędzia:
miękkie klastrowanie - jeden obiekt może trafić do jednego lub więcej klastró
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

# Part I - downloading proteomes from NCBI

import os
from Bio import Entrez, SeqIO

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
def change_protein_names(accessions_to_names_dict, ncbi_dir):
    dicts_list = []
    for filename in os.listdir(ncbi_dir):
        if filename.endswith(".faa"):
            ncbi_id = os.path.splitext(filename)[0] # without extension
            old_file = os.path.join(ncbi_dir, filename)
            dict = {}
            with open(old_file, 'r+') as f:
                lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith('>'):
                        words = line.split()
                        protein_id = words[0].split('>')[1]
                        genome_protein_id = f'{ncbi_id}_{protein_id}'
                        dict[genome_protein_id] = i
                        lines[i] = f'>{i}_{accessions_to_names_dict[ncbi_id]}\n'
                f.seek(0)
                f.truncate()
                for line in lines:
                    f.write(line)

            dicts_list.append(dict)
    return dicts_list


if __name__ == '__main__':
    ncbi_dir = './data/ncbi_genomes'
    accessions_to_names_dict = parse_genomes_ids('genomes.txt')
    # get_ncbi_genomes(genome_accessions, ncbi_dir)
    # gb_to_fasta(ncbi_dir)
    list_of_genome_protein_ids_names_dicts = change_protein_names(accessions_to_names_dict, ncbi_dir)
