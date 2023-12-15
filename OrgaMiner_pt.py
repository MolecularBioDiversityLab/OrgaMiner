import argparse 
import subprocess
import os
import gzip
import pandas as pd 
from datetime import datetime
import logging

log_format = '%(asctime)s - %(levelname)s - %(message)s'
date = str(datetime.now().strftime("%H%M-%y%m%d"))
file = f"info-{date}.log"

usage = """OrgaMiner_pt.py [-h] [taxa_file] [--pt-check] [--download {fasterq-dump,aspera,curl}] [--DNA ['-F', '-t', '-s', '--max-reads', '-P'] | --RNA]"""

parser = argparse.ArgumentParser(description='Process DNA or RNA data for plastid genome', usage=usage, add_help=False, epilog="-"*50)
parser.add_argument("-h", "--help", action="store_true", help="Show this help message and exit")
parser.add_argument('taxa_file', nargs="?", help= 'Input text file containing taxa, SRA and reference information.')
parser.add_argument("--pt-check", action="store_true", help="Discard species whose plastid genome is already assembled and available on NCBI.")
parser.add_argument('--download', help = "Choose one of the given download options. (default = fasterq-dump)", choices=["fasterq-dump", "aspera", "curl"], default="fasterq-dump")
parser.add_argument("--skip-download", action ="store_true", help="Skip download process to use present fastq files.")
parser.add_argument("--skip-trim", action ="store_true", help="Skip trimming process.")
parser.add_argument("-R", "--remove", action="store_true", help = "Remove fastq files after each process to save memory.")

parser._optionals.title = "General options"
args, remaining_args = parser.parse_known_args()

requirement = not parser.parse_known_args()[0].help and args.taxa_file

seq_group_parser = parser.add_mutually_exclusive_group(required=requirement)
seq_group_parser.add_argument("--DNA", action="store_true", help='Run GetOrganelle script for DNA data')
seq_group_parser.add_argument("--RNA", action="store_true", help='Process RNA reads')
args, remaining_args = parser.parse_known_args()

dna_parser = argparse.ArgumentParser(prog="OrgaMiner.py [taxa_file] --DNA", add_help=False, epilog="-"*50)
dna_parser.add_argument("-F", type=str, required = True, help = "Target organelle genome type(s). mandatory", choices=["embplant_pt", "other_pt"])
dna_parser.add_argument('-t', type=int,  metavar="THREADS", help='Number of threads to use for GetOrganelle.py (default=1)', default=1) 
dna_parser.add_argument('-s', type=str, metavar="SEED_FILE", help='Input fasta format file as initial seed') # !!!
dna_parser.add_argument('--max-reads', type=int, metavar="MAX_READS", help=' Maximum number of reads to be used per file. (default=1.5E7)', default=1.5E7) 
dna_parser.add_argument('-P', type=int, metavar="PRE_GROUPED", help='Pre-grouping value. Default: 200000', default=200000)
dna_parser._optionals.title = "Options for DNA"


def help_message():
    parser.print_help()
    dna_parser.print_help()

if args.help:
    help_message() 
else:
    if not args.taxa_file:
        parser.error("You must provide a path to the input file.")
    with open(args.taxa_file, 'r') as f:
        logging.basicConfig(filename=file, format=log_format, level=logging.DEBUG)

        print("""
                     <:::::::::::::::::::::::::::::::::::::::::::::::::::>
                    =|   ____                  __  __ _                  |=
                    =|  / __ \                |  \/  (_)                 |=
                    =| | |  | |_ __ __ _  __ _| \  / |_ _ __   ___ _ __  |=
                    =| | |  | | '__/ _` |/ _` | |\/| | | '_ \ / _ \ '__| |= 
                    =| | |__| | | | (_| | (_| | |  | | | | | |  __/ |    |= 
                    =|  \____/|_|  \__, |\__,_|_|  |_|_|_| |_|\___|_|    |= 
                    =|              __/ |                                |= 
                    =|             |___/                                 |=
                     <:::::::::::::::::::::::::::::::::::::::::::::::::::>
""")

        expected_line_length = None
        is_sra_given = False

        sra_species = {}
        species_list = []
        sra_list = []

        if args.DNA:
            dna_args = dna_parser.parse_args(remaining_args)
            for line in f:
                line = line.strip()
                line = line.split()
                taxa = line[0]

                if expected_line_length is None:
                    expected_line_length = len(line)

                if len(line) != expected_line_length:
                    error_message = "Inconsistent number of column in the input file. The presence or absence of SRA information should be consistent across all lines."
                    logging.error(f"{error_message}")
                    print(f"{error_message}")
                    raise Exception(error_message)
                
                if args.skip_download:
                    if len(line) != 2:
                        error_message = "If --skip-download is enabled, you must provide both species names and fastq file paths in the input file."
                        logging.error(f"{error_message}")
                        print(f"{error_message}")
                        raise Exception(error_message)
                    else:
                        species_list.append(taxa)
                        sra = line[1]
                        sra_list.append(sra)
                        sra_species[sra] = taxa
                else:
                    if len(line) == 2: # taxa, sra
                        species_list.append(taxa)
                        is_sra_given = True
                        sra = line[1]
                        sra_list.append(sra)
                        sra_species[sra] = taxa
                    else:
                        logging.info(f"Starting SRA Accession Search: Initiating ESearch for {taxa}")
                        print(f"Starting SRA Accession Search: Initiating ESearch for {taxa}")
                        cmd = f"esearch -db sra -query '{taxa}[Organism] AND GENOMIC[Source] AND WGS[Strategy]' | efetch -format runinfo"
                        output = subprocess.check_output(cmd, shell=True, text=True)
                        logging.info(f"SRA Accession Search Completed: ESearch Finished for {taxa}")
                        print(f"SRA Accession Search Completed: ESearch Finished for {taxa}")
                        with open("SRA_meta.txt", 'a') as sra_file:
                            sra_file.write(output)

        else: #RNA
            for line in f:
                line = line.strip()
                line = line.split()
                taxa = line[0] 

                if expected_line_length is None:
                    expected_line_length = len(line)

                if len(line) != expected_line_length:
                    error_message = "Inconsistent number of column in the input file. The presence or absence of SRA/Reference information should be consistent across all lines."
                    logging.error(f"{error_message}")
                    print(f"{error_message}")
                    raise Exception(error_message)
                
                if args.skip_download:
                    species_list.append(taxa)
                    if len(line) == 1:
                        error_message = "If --skip-download is enabled, you must provide both species names and fastq file paths in the input file."
                        logging.error(f"{error_message}")
                        print(f"{error_message}")
                        raise Exception(error_message)
                    elif len(line) == 2: # taxa, sra
                        sra = line[1]
                        sra_list.append(sra)
                        sra_species[sra] = taxa
                   
                else:
                    if len(line) == 2: # taxa, sra
                        is_sra_given = True
                        sra = line[1]
                        sra_list.append(sra)
                        sra_species[sra] = taxa
                    else: # taxa
                        logging.info(f"Starting SRA Accession Search: Initiating ESearch for {taxa}")
                        print(f"Starting SRA Accession Search: Initiating ESearch for {taxa}")
                        cmd = f"esearch -db sra -query '{taxa}[Organism] AND TRANSCRIPTOMIC[Source] AND RNA-Seq[Strategy]' | efetch -format runinfo"
                        output = subprocess.check_output(cmd, shell=True, text=True)
                        logging.info(f"SRA Accession Search Completed: ESearch Finished for {taxa}")
                        print(f"SRA Accession Search Completed: ESearch Finished for {taxa}")
                        with open("SRA_meta.txt", 'a') as sra_file:
                            sra_file.write(output)

################# READ SELECTING FROM SRA_meta.txt
    if is_sra_given:
        with open('sra_species.txt', 'w') as file:
            for key, value in sra_species.items():
                file.write(f'{key} {value}\n')

        with open('sra.txt', 'w') as f:
            f.write('\n'.join(sra_list))

        with open('species.txt', 'w') as f:
            f.write('\n'.join(list(dict.fromkeys(species_list))))

        os.makedirs("result")

    else:
        df = pd.read_csv('SRA_meta.txt', header=None)
        df.columns = df.iloc[0]
        df = df.drop_duplicates()
        df = df[df['ScientificName'].str.contains(" ")] # discards only genus entries
        df = df[~df['ScientificName'].str.contains(r' x |subsp\.|sp\.|cf\.|[0-9]', case=False, na=False)] # discards entries includes these expressions
        if args.DNA:
            df = df[(df["LibrarySelection"] == "PCR") | (df["LibrarySelection"] == "HYBRID") | (df["LibrarySelection"] == "RANDOM") | (df["LibrarySelection"] == "RANDOM PCR") | (df["LibrarySelection"] == "Restriction Digest")]
        else:
            df = df[~df['LibrarySelection'].str.contains(r'Oligo', case=False, na=False)] # discards Oligo_dt RNA-Seq reads
        df['bases'] = pd.to_numeric(df['bases'])
        df['size_MB'] = pd.to_numeric(df['size_MB'])
        df['avgLength'] = pd.to_numeric(df['avgLength'])
        df = df[df["avgLength"] < 1000]
        df = df[~((df['Platform'].str.contains('NANOPORE|PACBIO')) & (df['avgLength'] == 0))]
        df_sorted = df.sort_values(by=['ScientificName', 'bases'], ascending=[True, False])
        df2 = df_sorted.groupby('ScientificName')


        def check_size(group):
            """
            If there are reads higher than 15 Gb, checks if the lower ones reach at least 4 Gb and if does uses them, if does not takes the
            lowest value higher than 15Gb. If there is not any read higher than 15 Gb, takes all of them.
            """
            group = group[group["size_MB"] != 0]

            df_1 = group[group["size_MB"] > 15000]
            df_2 = group[group["size_MB"] < 15000]
            sum = df_2["size_MB"].sum()

            if df_1["size_MB"].min() > 0:
                if sum > 4000:
                    df = df_2
                else:
                    df = df_1[df_1["size_MB"] == df_1["size_MB"].min()]
            else :
                df = group

            return df

        results = []
        for name, group in df2:
            current_sum = 0
            rows_to_keep = []

            group = check_size(group)

            for index, row in group.iterrows():
                current_sum += row['size_MB']
                if current_sum < 5000:
                    rows_to_keep.append(index)
                else:
                    rows_to_keep.append(index)
                    break

            results.append(group.loc[rows_to_keep])


        df_final = pd.concat(results)
        
        os.system("mv SRA_meta.txt SRA_meta_initial.txt")
        df_final.to_csv("SRA_meta.txt", index=False)

        with open("SRA_meta.txt", "r") as file:
            for line in file:
                line = line.replace(",", "\t")

                columns = line.split("\t")
                column1 = columns[0] # Run
                column29 = columns[28] # ScientificName

                column1 = column1.replace(" ", "_")
                column29 = column29.replace(" ", "_")

                if "ScientificName" not in line :
                    sra_species[column1] = column29 

        sra_list = list(sra_species.keys())
        species_list = list(dict.fromkeys(sra_species.values()))

        with open('sra_species.txt', 'w') as file:
            for key, value in sra_species.items():
                file.write(f'{key} {value}\n')

        with open('sra.txt', 'w') as f:
            f.write('\n'.join(sra_list))

        with open('species.txt', 'w') as f:
            f.write('\n'.join(species_list))

################# MITOGENOME CHECK 

        os.makedirs("result")
        os.system("mv SRA_meta_initial.txt result/")

        if args.pt_check:
            logging.info("Plastid Genome Check Begins: Identifying Taxa with Complete Plastid Genomes")
            print("Plastid Genome Check Begins: Identifying Taxa with Complete Plastid Genomes")
            import check
            check.pt_check("species.txt")
            logging.info("Plastid Genome  Check Completed: Taxa with Complete Plastid Genomes Identified >> result/pt_species.txt")
            print("Plastid Genome Check Completed: Taxa with Complete Plastid Genomes Identified >> result/pt_species.txt")
            species_pt = open("pt_species.txt", "r").read().split("\n")
            sra_species2 = {}
            for sra,sp in sra_species.items():
                if sp not in species_pt:
                    sra_species2[sra] = sp

            sra_species = sra_species2 # remove species that are present in pt_species.txt
            del sra_species2

            df = pd.read_csv('SRA_meta.txt')
            df['ScientificName_t'] = df['ScientificName'].str.replace(' ', '_')
            df = df[df['ScientificName_t'].isin(species_pt)] 
            df = df.drop(columns=['ScientificName_t'])
            os.system("rm SRA_meta.txt")
            df.to_csv("SRA_meta.txt", index=False)

            os.system("mv pt_species.txt result/")
            os.system("mv species.txt all_species.txt && mv all_species.txt result/")
            os.system("mv sra.txt all_sra.txt && mv all_sra.txt result/") 
            # result/all* files contain both species/reads with and without plastid genome

            sra_list = list(sra_species.keys())
            species_list = list(dict.fromkeys(sra_species.values()))

            os.remove("sra_species.txt")

            with open('sra_species.txt', 'w') as file:
                for key, value in sra_species.items():
                    file.write(f'{key} {value}\n')

            with open('species.txt', 'w') as f:
                f.write('\n'.join(species_list))

            with open('sra.txt', 'w') as f:
                f.write('\n'.join(sra_list))

    species_sras = {} # {"spe1" : ["sra1", "sra2", "sra3"]}

    for sra, species in sra_species.items():

        if species not in species_sras:
            species_sras[species] = []
        species_sras[species].append(sra)


    for SPECIES, SRAS in species_sras.items(): # ONE BY ONE

        logging.info(f"Starting Process for {SPECIES}")
        print(f"Starting Process for {SPECIES}")

        os.makedirs(f"{SPECIES}")

        os.system(f"mkdir result/{SPECIES}") 

        if not args.skip_download:
            sras = " ".join(SRAS)
            
            if args.download == "aspera":
                command = f"kingfisher get -r {sras} -m ena-ascp ena-ftp"
            elif args.download == "curl":
                command = f"kingfisher get -r {sras} -m ena-ftp"
            else:
                command = f"kingfisher get -r {sras} -m prefetch"

            os.chdir(f"{SPECIES}")

            logging.info(f"{SPECIES} - Starting Download for SRA Data")
            print(f"{SPECIES} - Starting Download for SRA Data")

            os.system(command)

            logging.info(f"{SPECIES} - SRA Data Download Completed")
            print(f"{SPECIES} - SRA Data Download Completed")

            os.chdir("..")

        else:
            for sra in SRAS:
                os.system(f"mv {sra}* {SPECIES}")

################# FASTQC STAT

        logging.info(f"{SPECIES} - Generating Sequence Statistics ...")
        print(f"{SPECIES} - Generating Sequence Statistics ...")
        os.system(f"seqkit stat -a {SPECIES}/* >> result/{SPECIES}/stats.txt")
        logging.info(f"{SPECIES} - Sequence Statistics Written")
        print(f"{SPECIES} - Sequence Statistics Written")

        logging.info(f"{SPECIES} - Running FastQC ...")
        print(f"{SPECIES} - Running FastQC ...")
        os.system(f"fastqc {SPECIES}/*")
        logging.info(f"{SPECIES} - FastQC Completed ")
        print(f"{SPECIES} - FastQC Completed ")

        ################# QUALITY REPORT

        logging.info(f"{SPECIES} - MultiQC Analysis Begins")
        print(f"{SPECIES} - MultiQC Analysis Begins")
        os.system(f"multiqc {SPECIES} -o {SPECIES}")
        os.system(f"mv {SPECIES}/multiqc* result/{SPECIES}") 
        logging.info(f"{SPECIES} - MultiQC Analysis Completed")
        print(f"{SPECIES} - MultiQC Analysis Completed")

        os.system(f"mv {SPECIES}/*fastqc* result/{SPECIES}")
        os.system(f"mv {SPECIES}/*multiqc* result/{SPECIES}")

################# ADAPTER TRIMMING

        if not args.skip_trim:

            for root, dirs, files in os.walk(SPECIES):
                files = sorted(files)
                for file in files:

                    if "_1" in file and "fastqc" not in file:
                        reverse = file.replace("_1", "_2")
                        if reverse in files:
                            trim_command =  f"trim_galore --gzip --paired -o {SPECIES} {SPECIES}/{file} {SPECIES}/{reverse}"

                            logging.info(f"{SPECIES} - Adapter Trimming Begins for {file} {reverse}")
                            print(f"{SPECIES} - Adapter Trimming Begins for {file} {reverse}")
                            subprocess.call(trim_command, shell=True)
                            logging.info(f"{SPECIES} - Adapter Trimming Completed for {file} {reverse}")
                            print(f"{SPECIES} - Adapter Trimming Completed for {file} {reverse}")

                            if args.remove: 
                                os.system(f"rm {SPECIES}/{file} {SPECIES}/{reverse}")
                            else:
                                os.system(f"mv {SPECIES}/{file} {SPECIES}/{reverse} result/{SPECIES}")

                        else:
                            trim_command = f"trim_galore --gzip -o {SPECIES} {SPECIES}/{file}"

                            logging.info(f"{SPECIES} - Adapter Trimming Begins for {file}")
                            print(f"{SPECIES} - Adapter Trimming Begins for {file}")
                            subprocess.call(trim_command, shell=True)
                            logging.info(f"{SPECIES} - Adapter Trimming Completed for {file}")
                            print(f"{SPECIES} - Adapter Trimming Completed for {file}")

                            if args.remove: 
                                os.system(f"rm {SPECIES}/{file}")
                            else:
                                os.system(f"mv {SPECIES}/{file} result/{SPECIES}")

                    elif "_1" not in file and "_2" not in file and "fastqc" not in file:
                        trim_command = f"trim_galore --gzip -o {SPECIES} {SPECIES}/{file}"

                        logging.info(f"{SPECIES} - Adapter Trimming Begins for {file}")
                        print(f"{SPECIES} - Adapter Trimming Begins for {file}")
                        subprocess.call(trim_command, shell=True)
                        logging.info(f"{SPECIES} - Adapter Trimming Completed for {file}")
                        print(f"{SPECIES} - Adapter Trimming Completed for {file}")

                        if args.remove: 
                            os.system(f"rm {SPECIES}/{file}")
                        else:
                            os.system(f"mv {SPECIES}/{file} result/{SPECIES}")

            os.system(f"mv {SPECIES}/*txt result/{SPECIES}")


    ################# SEQKIT STAT FOR TRIMMED READS

            logging.info(f"{SPECIES} - Generating Sequence Statistics for Trimmed Reads ...")
            print(f"{SPECIES} - Generating Sequence Statistics for Trimmed Reads ...")
            os.system(f"seqkit stat -a {SPECIES}/*gz >> result/{SPECIES}/stats.txt")
            logging.info(f"{SPECIES} - Sequence Statistics for Trimmed Reads Completed")        
            print(f"{SPECIES} - Sequence Statistics for Trimmed Reads Completed")        

################# COMBINING READS

        current_dir = os.getcwd()
        gz_file_paths_1 = {}
        gz_file_paths_2 = {}
        gz_file_paths_single = {}

        for root, dirs, files in os.walk(f"{SPECIES}"):
            files = sorted(files)
            for file in files:
                if "_1" not in file and "_2" not in file:
                    folder_path = os.path.relpath(root, current_dir)  # Get the relative path of the folder
                    gz_file_paths_single.setdefault(folder_path, []).append(os.path.join(root, file))

                if "_1" in file :
                    reverse = file.replace("_1", "_2")
                    if reverse in files:
                        folder_path = os.path.relpath(root, current_dir)  # Get the relative path of the folder
                        gz_file_paths_1.setdefault(folder_path, []).append(os.path.join(root, file))
                    else:
                        folder_path = os.path.relpath(root, current_dir)  # Get the relative path of the folder
                        gz_file_paths_single.setdefault(folder_path, []).append(os.path.join(root, file))

                if '_2' in file :
                    folder_path = os.path.relpath(root, current_dir)  # Get the relative path of the folder
                    gz_file_paths_2.setdefault(folder_path, []).append(os.path.join(root, file))


        for folder_path, gz_file_list in gz_file_paths_single.items():
            leng = len(gz_file_list)
            if leng > 0:
                if leng == 1:
                    logging.info(f"{SPECIES} - Only One Single-end Read found, Passing Combining")
                    print(f"{SPECIES} - Only One Single-end Read found, Passing Combining")

                    if gz_file_list[0].endswith("gz") :
                        os.system(f"mv {gz_file_list[0]} {SPECIES}/combined_single.fastq.gz")
                    else:
                        os.system(f"mv {gz_file_list[0]} {SPECIES}/combined_single.fastq")

                else:
                    logging.info(f"{SPECIES} - {leng} Single-end Reads Found, Combining Begins")
                    print(f"{SPECIES} - {leng} Single-end Reads Found, Combining Begins")

                    if gz_file_list[0].endswith("gz"):    
                        output_filename = os.path.join(current_dir, folder_path, 'combined_single.fastq.gz')
                        with gzip.open(output_filename, 'wb') as output_file:
                            for gz_file_path in gz_file_list:
                                with gzip.open(gz_file_path, 'rb') as input_file:
                                    output_file.write(input_file.read())

                    else:
                        output_filename = os.path.join(current_dir, folder_path, 'combined_single.fastq.gz')
                        with gzip.open(output_filename, 'wb') as output_file:
                            for gz_file_path in gz_file_list:
                                with open(gz_file_path, 'rb') as input_file:
                                    output_file.write(input_file.read())

                    logging.info(f"{SPECIES} - Single-end Read Combining Is Completed")
                    print(f"{SPECIES} - Single-end Read Combining Is Completed")

                    if args.remove:
                        os.system(f"rm {' '.join(gz_file_list)}")
                    else:
                        os.system(f"mv {' '.join(gz_file_list)} result/{SPECIES}")
                    
            else:
                logging.info(f"{SPECIES} - Not Found Any Single-end Read")
                print(f"{SPECIES} - Not Found Any Single-end Read")


        for folder_path, gz_file_list in gz_file_paths_1.items():
            leng = len(gz_file_list)
            if leng > 0:
                if leng == 1:
                    logging.info(f"{SPECIES} - Only One Forward Read Found, Passing Combining")
                    print(f"{SPECIES} - Only One Forward Read Found, Passing Combining")

                    if gz_file_list[0].endswith("gz") :
                        os.system(f"mv {gz_file_list[0]} {SPECIES}/combined_1.fastq.gz")
                    else:
                        os.system(f"mv {gz_file_list[0]} {SPECIES}/combined_1.fastq")

                else:
                    logging.info(f"{SPECIES} - {leng} Forward Reads Found, Combining Begins")
                    print(f"{SPECIES} - {leng} Forward Reads Found, Combining Begins")

                    if gz_file_list[0].endswith("gz"):    
                        output_filename = os.path.join(current_dir, folder_path, 'combined_1.fastq.gz')
                        with gzip.open(output_filename, 'wb') as output_file:
                            for gz_file_path in gz_file_list:
                                with gzip.open(gz_file_path, 'rb') as input_file:
                                    output_file.write(input_file.read())

                    else:
                        output_filename = os.path.join(current_dir, folder_path, 'combined_1.fastq.gz')
                        with gzip.open(output_filename, 'wb') as output_file:
                            for gz_file_path in gz_file_list:
                                with open(gz_file_path, 'rb') as input_file:
                                    output_file.write(input_file.read())

                    logging.info(f"{SPECIES} - Forward Read Combining Is Done")
                    print(f"{SPECIES} - Forward Read Combining Is Done")

                    if args.remove:
                        os.system(f"rm {' '.join(gz_file_list)}")
                    else:
                        os.system(f"mv {' '.join(gz_file_list)} result/{SPECIES}")
                    
            else:
                logging.info(f"{SPECIES} - Not Found Any Forward Read")
                print(f"{SPECIES} - Not Found Any Forward Read")


        for folder_path, gz_file_list in gz_file_paths_2.items():
            leng = len(gz_file_list)
            if leng > 0:
                if leng == 1:
                    logging.info(f"{SPECIES} - Only One Reverse Read Found, Passing Combining")
                    print(f"{SPECIES} - Only One Reverse Read Found, Passing Combining")

                    if gz_file_list[0].endswith("gz") :
                        os.system(f"mv {gz_file_list[0]} {SPECIES}/combined_2.fastq.gz")
                    else:
                        os.system(f"mv {gz_file_list[0]} {SPECIES}/combined_2.fastq")

                else:
                    logging.info(f"{SPECIES} - {leng} Reverse Reads Found, Combining Begins")
                    print(f"{SPECIES} - {leng} Reverse Reads Found, Combining Begins")

                    if gz_file_list[0].endswith("gz"):    
                        output_filename = os.path.join(current_dir, folder_path, 'combined_2.fastq.gz')
                        with gzip.open(output_filename, 'wb') as output_file:
                            for gz_file_path in gz_file_list:
                                with gzip.open(gz_file_path, 'rb') as input_file:
                                    output_file.write(input_file.read())

                    else:
                        output_filename = os.path.join(current_dir, folder_path, 'combined_2.fastq.gz')
                        with gzip.open(output_filename, 'wb') as output_file:
                            for gz_file_path in gz_file_list:
                                with open(gz_file_path, 'rb') as input_file:
                                    output_file.write(input_file.read())

                    logging.info(f"{SPECIES} - Reverse Read Combining Is Done")
                    print(f"{SPECIES} - Reverse Read Combining Is Done")

                    if args.remove:
                        os.system(f"rm {' '.join(gz_file_list)}")
                    else:
                        os.system(f"mv {' '.join(gz_file_list)} result/{SPECIES}")
                    
            else:
                logging.info(f"{SPECIES} - Not Found Any Reverse Read")
                print(f"{SPECIES} - Not Found Any Reverse Read")


################# ASSEMBLY

        forward = " "
        reverse = " "
        single = " "

        for root, dirs, files in os.walk(SPECIES):
            if "combined_1.fastq.gz" in files or "combined_1.fastq" in files:
                forward = f"-1 combined_1.fastq*"
            if "combined_2.fastq.gz" in files or "combined_2.fastq" in files:
                reverse = f"-2 combined_2.fastq*"
            if "combined_single.fastq.gz" in files or "combined_single.fastq" in files:
                if args.DNA:
                    single = f"-u combined_single.fastq*"
                if args.RNA:
                    single = f"-S combined_single.fastq*"


        if args.DNA:
            suffix = ".fasta"
            logging.info(f"{SPECIES} - DNA Assembly Begins")
            print(f"{SPECIES} - DNA Assembly Begins")
            
            if dna_args.s:
                cmd = f"get_organelle_from_reads.py {forward} {reverse} {single} -F {dna_args.F} -o assembly -t {dna_args.t} --max-reads {dna_args.max_reads} -P {dna_args.P} -s {args.s}"
            else:
                cmd = f"get_organelle_from_reads.py {forward} {reverse} {single} -F {dna_args.F} -o assembly -t {dna_args.t} --max-reads {dna_args.max_reads} -P {dna_args.P}"
            
            os.chdir(f"{SPECIES}")
            subprocess.call(cmd, shell=True)

            logging.info(f"{SPECIES} - DNA Assembly Is Completed")
            print(f"{SPECIES} - DNA Assembly Is Completed")

            os.chdir("..")
            os.system(f"mv {SPECIES}/assembly/*fasta {SPECIES}")

        if args.remove:
            os.system(f"rm {SPECIES}/combined*")

        os.system(f"mv {SPECIES}/* result/{SPECIES}")
        os.system(f"rmdir {SPECIES}")