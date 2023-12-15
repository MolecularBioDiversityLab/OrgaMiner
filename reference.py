from Bio import Entrez
from Bio import SeqIO
import os
import time
Entrez.email = 'youremail@example.com'

def ref_finder(genus, logging):

    max_retries = 3
    retries = 0
    
    while retries < max_retries:
        try:
            if not os.path.exists(f"references/{genus}.fa"):
                query = f'{genus}[Organism] AND mitochondrion[filter] AND 10000:25000[Sequence Length]'

                handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)

                record = Entrez.read(handle)
                handle.close()

                sequence_ids = record['IdList']

                if len(sequence_ids) > 0:
                    sequence_ids = record['IdList'][0]

                    handle = Entrez.efetch(db="nucleotide", id=sequence_ids, rettype="fasta", retmode="text")

                    record = SeqIO.read(handle, "fasta")

                    fasta = ""
                    header = ">" + record.description + "\n"
                    sequence = str(record.seq)
                    fasta = header + sequence

                    with open(f"references/{genus}.fa", "a") as file:
                        file.write(fasta)

                    return True

                else:
                    return False
                
            else:
                return True
            
        except Exception as e:
            if logging:
                logging.error(f"An error occurred while searching for reference fasta : {e}")
            retries += 1
            if logging:
                logging.info(f"Retry attempt {retries}/{max_retries}")
            time.sleep(5)
        
