# UCT_metagenomics

This pipeline was developed for the metagenomic classification Illumina paired-end reads through the uct_meta.py script. 
Targeted read evaluation can also be performed using the target_read_retrieval.py script for the specific further 
evaluation of reads classified by kraken2.

## Dependencies

This pipeline is run in a conda environment which can be created using the file uct_meta_env.yml.

## Usage - uct_meta.py

The uct_meta.py script is run using python. The script takes as input a directory with paired-end reads. 
The filename convention is as follows: 

`<sample_name>_R<1/2>.fastq.gz`

The number of chars in the sample name should be consistent (including negative control)

Basic usage is as follows:

`uct_meta.py -e <string_length_identifying_reads> -i <input_folder> -o <output_folder> -t <threads> -m <trimmomatic_path> -1 <path_to_kraken2_db1> -2 <path_to_kraken2_db2> -d <path_to_fasta_file_for_db_creation> -k <path_to_kneaddata_db> -g <string_present_in_all_unprocessed_reads> '-f' <string_present_in_all_forward_reads> '-r' <string_present_in_all_reverse_reads> '-a' <metadata_file>`

The script takes the following arguments:

- `--identifier_length, -e` - (Required) The string length of the input reads' filename that identifies the pair of 
reads as a unit. For example: [sample_1_R1.fastq.gz, sample_1_R2.fastq.gz] the identifier length would be 8.

- `--metadata, -a` - (Required) The filepath to the metadata file.

- `--input_folder, -i` - (Default=CWD) The path to the folder containing the paired-end reads to be processed.

- `--output_folder, -o` - (Default=CWD) The path to the folder where the output files will be written.

- `--threads, -t` - (Default=1) The number of threads to be used for the analysis.

- `--trimmomatic, -m` - (Default=bin) The path to the trimmomatic jar file.

- `--kraken2_db1, -1` - (Default=Path specified in data_and_variable.py file) The path to the first kraken2 
database to use for the analysis. The analysis allows for usage of up to 2 databases.

- `--kraken2_db2, -2` - (Default=Path specified in data_and_variable.py file) The path to the second kraken2 
database to use for the analysis.

- `--diamond_db, -d` - (Default=Path specified in data_and_variable.py file) The path to the fasta file to use 
to construct a DIAMOND database.

- `--kneaddata_db, -k` - (Default=Path specified in data_and_variable.py file) The path to the kneaddata database 
to use for human read subtraction and contaminant removal.

- `--gen_id', '-g` - (Default='fastq.gz') A string identifying all unprocessed reads.

- `--r1_id, -f` - (Default='_R1.fastq.gz') A string identifying all forward reads.

- `---r2_id, -r` - (Default='_R2.fastq.gz') A string identifying all reverse reads.

The metadata file (-a) is csv with the headers:

`sample,run`

The sample name should correspond to the substring generated when the read filenames are split by the identifier 
length (-e). Run is an integer. 

## Usage - target_read_retrieval.py

This script should be run once the uct_meta.py script has been run and potentially significant pathogens have been 
selected for further follow-up. Prior to running target_read_retrieval:

- Generate a csv file for targets to follow up with the following headers

`sample,organism,taxid,reference,reference_prefix,kraken_database`

The sample name should correspond to the substring generated when the read filenames are split by the identifier. 
For example, the sample name for [sample_1_R1.fastq.gz, sample_1_R2.fastq.gz] would be <sample_1>. The organism is the 
organism name as it appears in the Kraken 2 report. The taxid is the NCBI Taxonomy ID for the reads to be retrieved as 
it appears in the Kraken 2 report. The reference is the name of the reference genome as it appears in the specified 
references folder. The reference_prefix is the prefix to be used for the output files. The kraken_database is the 
database that was used to detect the targeted pathogen (either 1 or 2).

- Reference genomes for the targets for further evaluation should be retrieved and placed in a folder which is passed 
to the script and specified in the reference csv specified above.


