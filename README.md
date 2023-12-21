# Introduction

**OrgaMiner** is a pipeline designed for automated assembly of mitochondrial and plastid genomes of any taxa from the reads available in the NCBI SRA database. This pipeline comprises three essential stages: *"Downloading_SRA"*, *"Trimming_and_read_quality_assessment"*, *"Assembly"*, and *"Annotation"*. In the *"Downloading_SRA"* phase, desired SRAs (Sequence Read Archives) from the NCBI database are retrieved, employing a range of available options, which will be detailed later. Moving on to the *"Trimming_and_read_quality_assessment"* step, its primary objective is to enhance the quality of the data in the Fastq format. This is achieved by applying Trim Galore wrapper for trimming and FastQC for generating read quality reports for the SRAs. Then, quality reports are summarized using MultiQC. The *"Assembly"* stage is dedicated to the assembly of organelle genomes from the input DNA or RNA reads. Finally, the assembled mitochondrial genomes are annotated with MitoZ.

# Quick start

Dependencies:

- [python](https://www.python.org/) (v3.8 or higher)
- [kingfisher-download](https://github.com/wwood/kingfisher-download) (v.0.3.1)
- [Biopython](https://biopython.org/) (v1.78)
- [Pandas](https://pandas.pydata.org/) (v1.5.2)
- [cutadapt](https://github.com/marcelm/cutadapt) (v4.4)
- [entrez-direct](https://github.com/Klortho/edirect)
- [sra-tools](https://github.com/ncbi/sra-tools) (v3.0 or higher)
- [Seqkit](https://github.com/shenwei356/seqkit)
- [FastQC](https://github.com/s-andrews/FastQC)
- [MultiQC](https://github.com/ewels/MultiQC) (v1.11)
- [Trim Galore](https://github.com/FelixKrueger/TrimGalore) (v0.6.10)
- [GetOrganelle](https://github.com/Kinggerm/GetOrganelle) (v1.7.7.0) 
- [MITGARD](https://github.com/pedronachtigall/MITGARD)
- [Samtools](https://www.htslib.org/) (v1.9 or higher)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (v2.3.5.1 or higher)
- [Minimap2](https://github.com/lh3/minimap2) (v2.17 or higher)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) (v2.8.5 or higher)
- [SPAdes](https://cab.spbu.ru/software/spades/) (v3.13.1 or higher)
- [MitoZ](https://github.com/linzhi2013/MitoZ) (v3.6) 

**IMPORTANT NOTE**: Downloading and compiling the MITGARD repository is required for working with RNA data. Check [MITGARD github page](https://github.com/pedronachtigall/MITGARD) for further guidance.

**IMPORTANT NOTE**: When installing GetOrganelle, the preferred databases must be installed by using `get_organelle_config.py --add {database}` command. Check [GetOrganelle github page](https://github.com/Kinggerm/GetOrganelle) for further guidance.

**IMPORTANT NOTE**: Before running the script, please replace 'youremail@example.com' in the 5th line of check.py and reference.py with your email address for the Entrez.email parameter.

**WARNING**: Using `-–remove` option will cause all fastq files to be deleted permanently. Use with caution.

**NOTE**: Using mamba for installation of the dependencies is highly recommended.

For Aspera Installation (optional but recommended for downloading fastq files in .gz format):

`wget https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/0adrj/0/ibm-aspera-connect_4.1.3.93_linux.tar.gz`

`tar zxvf ibm-aspera-connect_4.1.3.93_linux.tar.gz`

`bash ibm-aspera-connect_4.1.3.93_linux.sh`

# Steps to follow

1. Install the dependencies and compile if needed (see above)
2. Download OrgaMiner repository
3. Generate a metadata file to define taxa, and include the SRA accession if selected by the user
4. Run OrgaMiner script with proper options (see below for further explanation)

# Usage

`OrgaMiner.py [-h]
 [taxa_file] [--mt-check] [--download {fasterq-dump,aspera,curl}] [--annotate {Chordata,Arthropoda,Echinodermata,Annelida-segmented-worms,Bryozoa,Mollusca,Nematoda,Nemertea-ribbon-worms,Porifera-sponges}] [--DNA ['-F', '-t', '-s', '--max-reads', '-P'] | --RNA ['-a', '-L', '-r', '-C', '-g', '-M', '']]`

Process DNA or RNA data for mitochondrial genome

**Positional arguments:**

    taxa_file             Input text file containing taxa, SRA and reference information.

**General options:**

    -h, --help          Show this help message and exit
    --mt-check          Discard species whose mitochondrial genome is already assembled and available on NCBI.
    --download {fasterq-dump,aspera,curl}
                        Choose one of the given download options. (default = fasterq-dump)
    --skip-download     Skip download process to use present fastq files.
    --skip-trim         Skip trimming process.
    -R, --remove        Remove fastq files after each process to save memory.
    --annotate {Chordata,Arthropoda,Echinodermata,Annelida-segmented-worms,Bryozoa,Mollusca,Nematoda,Nemertea-ribbon-worms,Porifera-sponges}
                        Annotate output assembly using MITOZ. Clade must be specified.


    --DNA   Run GetOrganelle script for DNA data

        Options for DNA / GetOrganelle:
	
        -F {embplant_mt,animal_mt}
                              Target organelle genome type(s). mandatory

        -t THREADS            Number of threads to use for GetOrganelle.py (default=1)

        -s SEED_FILE          Input fasta format file as initial seed

        --max-reads MAX_READS
                              Maximum number of reads to be used per file. (default=1.5E7)

        -P PRE_GROUPED        Pre-grouping value. Default: 200000


    --RNA   Run MITGARD.py script for RNA data

        Options for RNA / MITGARD:

        -a ASSEMBLER, --assembler ASSEMBLER
                            sets the assembler tools to be used to generate the mitochondrial contigs.
                            trinity,spades,mitoz. default=trinity,spades

        -L {N,R}, --low_coverage {N,R}
                    to decide what to do with low coverage regions. N to assign N's in the low coverage regions. Use R to assign the Reference nucleotides in the low coverage regions. (default=N)

        -r {True,False}, --rearrangement {True,False}
                                to turn on/off an additional step to check for rearrangements using the mitochondrial contigs generated. To turn on set True.

        -C CLADE, --clade CLADE
                                set the taxonomic clade to be used in the '--clade' parameter of MitoZ assembler

        -g GENETIC_CODE, --genetic_code GENETIC_CODE
                                set the genetic code to be used in the '--genetic_code' parameter of MitoZ assembler [Default=2]

        -M MEMORY, --memory MEMORY
                                Max memory usage to be passed to Trinity assembler [default=4G], use the same format as stated by Trinity assembler

        -c CPU, --cpu CPU       Number of threads to be used in each step of MITGARD.



`OrgaMiner_pt.py [-h] [taxa_file] [--pt-check] [--download {fasterq-dump,aspera,curl}] [--DNA ['-F', '-t', '-s', '--max-reads', '-P'] | --RNA]`

Process DNA or RNA data for plastid genome

**Positional arguments:**

    taxa_file             Input text file containing taxa, SRA and reference information.

**General options:**

    -h, --help            Show this help message and exit
    --pt-check            Discard species whose plastid genome is already assembled and available on
                        NCBI.
    --download {fasterq-dump,aspera,curl}
                          Choose one of the given download options. (default = fasterq-dump)
    --skip-download       Skip download process to use present fastq files.
    --skip-trim           Skip trimming process.
    -R, --remove          Remove fastq files after each process to save memory.

    --DNA                 Run GetOrganelle script for DNA data

        -F {embplant_pt,other_pt}
                            Target organelle genome type(s). mandatory
        -t THREADS            Number of threads to use for GetOrganelle.py (default=1)
        -s SEED_FILE          Input fasta format file as initial seed
        --max-reads MAX_READS
                            Maximum number of reads to be used per file. (default=1.5E7)
        -P PRE_GROUPED        Pre-grouping value. Default: 200000

    --RNA                 Process RNA reads

**Important Note**: OrgaMiner_pt.py script does not handle assembly step for RNA data. 

Example usage for OrgaMiner.py script with DNA data

`python OrgaMiner.py taxa_metafile.txt –-mt-check --download aspera –-annotate Mollusca –-remove --DNA -F animal_mt -t 16`
 
Example usage for OrgaMiner.py script with RNA data

`python OrgaMiner.py taxa_metafile.txt –-skip-download –-annotate Chordata –-skip-trim --RNA -a trinity,spades,mitoz –-clade Chordata`

Example usage for OrgaMiner_pt.py script with DNA data

`python OrgaMiner_pt.py taxa_metafile.txt --pt-check --DNA -F embplant_pt`

Example usage for OrgaMiner_pt.py script with RNA data

`python OrgaMiner_pt.py taxa_metafile.txt --download curl --RNA`


Example Metadata Files
---
**Warning**: When specifying species names in the input file, please ensure there are no spaces, and instead, use underscores (_) to separate words in the name. For example, 'Acanthopleura_granulata'

- Usage with Selected Taxa Only

If only the scientific name of taxa is to be given as input, the input meta file can be created like this:

    taxa1
	taxa2
	taxa3

- Usage with Selected SRA Accessions

If the SRA accessions to be used in assembly step already chosen by user, the input file can be created like this:

    species1	sra1
    species1	sra2
    species2	sra3

- Usage with Existing Fastq files

To process downloaded fastq files, ensure they are all located in the current working directory. Additionally, provide a list of file names in the input file, each corresponding to its respective taxa like this,

    species1	sra1_1.fastq.gz
	species1	sra1_2.fastq.gz
	species2	sra2.fastq.gz

or just specify the prefixes of those files without extensions.

    species1	sra1
	species2	sra2


To run the script with this input file, `–-skip-download` option must be given.

**Warning**: Be cautious when utilizing the `--remove` option with this usage, as it results in the deletion of all fastq files including the files you specify in the input file. Please be mindful of this consequence before proceeding with the command.

**Warning**: Pair-end read files must include "_1" or "_2" in their names.

- Usage with Specifying Reference Fasta Files (for --RNA)

If the species, corresponding read files and their reference fasta files to be used in assembly step are chosen by user before the process, these could be given in input file which the first column is taxa, the second is read files and the third is reference fasta files.

    species1	sra1.fastq	ref1.fasta
	species1	sra2.fastq	ref1.fasta
	species2	sra3.fastq	ref2.fasta

Then the script will use these fasta files instead of running reference.py to find and download a reference fasta.

Output Files and Directories
---

species.txt ->  List of all species came from esearch step (or the species specified by user)

sra.txt -> List of all SRA accessions used in the pipeline (or SRA accessions specified bu user)

sra_species.txt -> SRA accessions and corresponding species.

info-xxxx.log -> Log file of the process (xxxx stands for the date and time of the run)

SRA_meta.txt -> esearch results after a series of filters applied (results without filter is saved in result/SRA_meta_initial.txt)

result/ -> Stores all files and directories, excluding the ones mentioned above.

result/mt_species.txt or result/pt_species.txt -> list of species to discard whose organelle genome is present in NCBI database. (for --mt-check or --pt-check option)

result/references/ -> All the reference fasta files used in assembly (for --RNA option) 

All other output files specific to species is stored in its own folder in result directory.
