from Bio import Entrez
from Bio import SeqIO
import time

Entrez.email = 'youremail@example.com'

def mt_check(species_path):
    species_list = []
    with open(species_path, "r") as file:
        for line in file:
            columns = line.strip().split()
            if columns: 
                species_list.append(columns[0])

    species_list = list(dict.fromkeys(species_list))

    taxa_names = []

    for i in species_list:
        query = f'{i}[Organism] AND mitochondrion[filter] AND 8000:60000[Sequence Length]'

        handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)

        record = Entrez.read(handle)
        handle.close()

        sequence_ids = record['IdList']

        if len(sequence_ids) > 0:
            taxa_names.append(i)

        time.sleep(2)

    taxa_names = list(dict.fromkeys(taxa_names))

    with open("mt_species.txt", "a") as file:
        file.write("\n".join(taxa_names))

def pt_check(species_path):
 
    species_list = []
    with open(species_path, "r") as file:
        for line in file:
            columns = line.strip().split()
            if columns: 
                species_list.append(columns[0])

    species_list = list(dict.fromkeys(species_list))

    taxa_names = []

    for i in species_list:
        query = f'{i}[Organism] AND chloroplast[filter] AND 80000:600000[Sequence Length]'

        handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)  # You can adjust the retmax value

        record = Entrez.read(handle)
        handle.close()

        sequence_ids = record['IdList']

        if len(sequence_ids) > 0:
            taxa_names.append(i)

        time.sleep(2)

    taxa_names = list(dict.fromkeys(taxa_names))

    with open("pt_species.txt", "a") as file:
        file.write("\n".join(taxa_names))

