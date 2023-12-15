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
        query = f'{i}[Organism] AND mitochondrion[filter] AND 10000:25000[Sequence Length]'

        handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)

        record = Entrez.read(handle)
        handle.close()

        sequence_ids = record['IdList']

        for seq_id in sequence_ids:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")

            record = SeqIO.read(handle, "genbank")

            taxa = ""
            description = record.description
            description = description.split(" ")
            description = [x for x in description if x != "UNVERIFIED:" and "TPA" not in x]  
            taxa += description[0] + "_" + description[1]
            
            taxa_names.append(taxa)
            time.sleep(3)
    
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
        query = f'{i}[Organism] AND chloroplast[filter] AND 100000:280000[Sequence Length]'

        handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)  # You can adjust the retmax value

        record = Entrez.read(handle)
        handle.close()

        sequence_ids = record['IdList']

        for seq_id in sequence_ids:
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")

            record = SeqIO.read(handle, "genbank")

            taxa = ""
            description = record.description
            description = description.split(" ")
            description = [x for x in description if x != "UNVERIFIED:" and "TPA" not in x]  
            taxa += description[0] + "_" + description[1]
            taxa_names.append(taxa)
            time.sleep(3)

    taxa_names = list(dict.fromkeys(taxa_names))

    with open("pt_species.txt", "a") as file:
        file.write("\n".join(taxa_names))


