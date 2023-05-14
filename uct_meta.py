import argparse
import logging
import os
import services
import data_and_variables as _dv
import knead_data_run
import subprocess
import shutil
import kraken2_command_generation as kraken2
import extract_from_krakenfile as extract
import generate_meta_contigs_spades as spades
import interleave_paired_reads as interleave
import diamond_assign as diamond

# The input file for the pipeline
parser = argparse.ArgumentParser(description='This pipeline analyses metagenomic data. Created by GJK Marais.')
parser.add_argument('--input_folder', required=True, help='The file path that contains Illumina paired-end reads')
parser.add_argument('--output_folder', required=True, help='Result folders will be output to this folder')
parser.add_argument('--threads', default=1, help='The number of threads to use, default is 1')
parser.add_argument('--identifier_length', type=int, required=True,
                    help='The string length that identifies the sequencing read')
parser.add_argument('--trimmomatic', default=None, help='trimmomatic path')
parser.add_argument('--kraken2_db1', default=None, help='kraken2 db1 path')
parser.add_argument('--kraken2_db2', default=None, help='kraken2 db2 path')
parser.add_argument('--diamond_db', default=None, help='diamond db path')
parser.add_argument('--kneaddata_db', default=None, help='knead_data path')
args = parser.parse_args()

# Set variables from input
input_folder = args.input_folder
output_folder = args.output_folder
id_length = args.identifier_length
threads = args.threads
trimmomatic_path = args.trimmomatic
kraken2_db1_path = args.kraken2_db1
kraken2_db2_path = args.kraken2_db2
diamond_db_path = args.diamond_db
kneaddata_db_path = args.kneaddata_db

# Logging
logger = logging.getLogger(__name__)
handler = logging.FileHandler(os.path.join(output_folder, "UCT_metagenomics.log"))
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.debug("Metagenomics pipeline started")

# Assign database paths
if trimmomatic_path is None:
    trimmomatic_path = _dv.trimmomatic_path
if kraken2_db1_path is None:
    kraken2_db1_path = _dv.kraken2_db1_path
if kraken2_db2_path is None:
    kraken2_db2_path = _dv.kraken2_db2_path
if diamond_db_path is None:
    diamond_db_path = _dv.diamond_db_path
if kneaddata_db_path is None:
    kneaddata_db_path = _dv.kneaddata_db_path

# Generate a list with all valid filenames
file_names_all = services.filename_list_generate('001.fastq.gz', input_folder)

# Separate filenames by id and into a list of pair tuples
paired_sorted_dict = services.paired_read_list_generate(id_length,
                                                        'R1_001.fastq.gz',
                                                        'R2_001.fastq.gz',
                                                        file_names_all,
                                                        input_folder)

# Generate output directories
knead_directory = os.path.join(output_folder, 'knead_data_results')
kraken_source_directory = os.path.join(output_folder, 'kraken2_source_files')
kraken2_directory_db1 = os.path.join(output_folder, 'kraken2_db1_results')
kraken2_directory_db2 = os.path.join(output_folder, 'kraken2_db2_results')
extracted_reads_db1_directory = os.path.join(output_folder, 'extracted_reads_db1')
extracted_reads_db2_directory = os.path.join(output_folder, 'extracted_reads_db2')
diamond_output = os.path.join(output_folder, 'diamond_output')

extract_db_list = [extracted_reads_db1_directory, extracted_reads_db2_directory]

directories_to_create = [knead_directory,
                         kraken2_directory_db1,
                         kraken2_directory_db2,
                         kraken_source_directory,
                         extracted_reads_db1_directory,
                         extracted_reads_db2_directory,
                         diamond_output]

for directory in directories_to_create:
    if not os.path.exists(directory):
        os.makedirs(directory)
        logger.debug(f'Output directory created: {directory}')

# Generate kneaddata commands
kneaddata_commands = knead_data_run.generate_kneaddata_commands(paired_sorted_dict,
                                                                knead_directory,
                                                                threads,
                                                                trimmomatic_path,
                                                                kneaddata_db_path)

# Run kneaddata commands
for knead_command in kneaddata_commands:
    try:
        subprocess.run(knead_command, check=True)
        logging.debug(f'Kneaddata command run successful: {knead_command}')
    except subprocess.CalledProcessError as e:
        logging.debug(f'Kneaddata command failed: {knead_command}')

# Generate kraken2 source directory
kraken2_source_files = {}
folders_knead = os.listdir(knead_directory)
directory_paths_knead = [(os.path.join(knead_directory, folder), folder)
                         for folder in folders_knead
                         if os.path.isdir(os.path.join(knead_directory, folder))]
for directory in directory_paths_knead:
    logging.debug(f'Kneaddata read output retrieval: {directory}')
    filename_1 = f'{directory[1]}_L001_R1_001_kneaddata_paired_1.fastq'
    filename_2 = f'{directory[1]}_L001_R1_001_kneaddata_paired_2.fastq'
    file_1 = os.path.join(directory[0], filename_1)
    file_2 = os.path.join(directory[0], filename_2)
    copy_path_1 = os.path.join(kraken_source_directory, filename_1)
    copy_path_2 = os.path.join(kraken_source_directory, filename_2)
    if os.path.isfile(file_1) and os.path.isfile(file_2):
        logging.debug(f'Kraken2 source files retrieved from: {directory}')
        shutil.copy(file_1, kraken_source_directory)
        shutil.copy(file_2, kraken_source_directory)
    kraken2_source_files[directory[1]] = (copy_path_1, copy_path_2)

# Generate kraken2 commands
kraken2_commands_db1 = kraken2.generate_kraken2_commands(kraken2_source_files,
                                                         kraken2_directory_db1,
                                                         threads,
                                                         kraken2_db1_path)

kraken2_commands_db2 = kraken2.generate_kraken2_commands(kraken2_source_files,
                                                         kraken2_directory_db2,
                                                         threads,
                                                         kraken2_db2_path)

kraken_commands_list = [kraken2_commands_db1, kraken2_commands_db2]

# Run kraken2 commands
for db, kraken_commands in enumerate(kraken_commands_list):
    logging.debug(f'Kraken2 db{db+1} run')
    for command in kraken_commands:
        try:
            subprocess.run(command, check=True)
            logging.debug(f'Kraken2 command db{db+1} successful: {command}')
        except subprocess.CalledProcessError as e:
            logging.debug(f'Kraken2 command db{db+1} failed: {command}')

# Create krakenfile dictionary to be used to extract relevant krakenfiles
krakenfiles = {}
krakenresults_folders = [kraken2_directory_db1, kraken2_directory_db2]
for db, krakenfolder in enumerate(krakenresults_folders):
    krakenfile_dict = {}
    folders_results = os.listdir(krakenfolder)
    paths_results = [(os.path.join(krakenfolder, folder), folder)
                     for folder in folders_results
                     if os.path.isdir(os.path.join(krakenfolder, folder))]
    for directory in paths_results:
        logging.debug(f'Krakenfile retrieval for db{db+1}: {directory}')
        krakenfile_name = f'{directory[1]}_krakenfile'
        file_path_krakenfile = os.path.join(directory[0], krakenfile_name)
        query_file_path_1 = os.path.join(kraken_source_directory, f'{directory[1]}_L001_R1_001_kneaddata_paired_1.fasta')
        query_file_path_2 = os.path.join(kraken_source_directory, f'{directory[1]}_L001_R1_001_kneaddata_paired_2.fasta')
        if os.path.isfile(file_path_krakenfile):
            logging.debug(f'krakenfile added to krakenfiles dictionary: {directory[1]}')
        krakenfile_dict[directory[1]] = (file_path_krakenfile, query_file_path_1, query_file_path_2)
    krakenfiles[db] = krakenfile_dict

# Retrieve unassigned reads
for db, krakenfiles_dict in krakenfiles.items():
    logging.debug(f'Extracting unassigned reads for db{db+1}')
    for sample, krakenfile in krakenfiles_dict.items():
        logging.debug(f'Extracting unassigned reads for: {sample}')
        output_directory = os.path.join(extract_db_list[db], sample)
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
            logging.debug(f'Output directory created: {output_directory}')
        extract_reads_command = extract.extract_reads_command(krakenfile[0], 'Unassigned_reads', '0', krakenfile[1],
                                                              krakenfile[2], output_directory)
        try:
            subprocess.run(extract_reads_command, check=True)
            logging.debug(f'Extracting unassigned reads successful: {extract_reads_command}')
        except subprocess.CalledProcessError as e:
            logging.debug(f'Extracting unassigned reads failed: {extract_reads_command}')

# Generate contigs and interleaved file for unassigned reads
# Create output file structure
unassigned_dict = {}
for db, reads_folder in enumerate(extract_db_list):
    unassigned_local = {}
    folders_extract = os.listdir(reads_folder)
    paths_extract = [(os.path.join(reads_folder, folder), folder)
                     for folder in folders_extract
                     if os.path.isdir(os.path.join(reads_folder, folder))]
    for directory in paths_extract:
        logging.debug(f'Generating contigs for: {directory[1]}')
        output_directory = os.path.join(directory[0], 'Unassigned_reads')
        unassigned_read_1 = os.path.join(output_directory, 'Unassigned_reads_1.fastq')
        unassigned_read_2 = os.path.join(output_directory, 'Unassigned_reads_2.fastq')
        contigs_folder = os.path.join(output_directory, 'contigs')
        if not os.path.exists(contigs_folder):
            os.makedirs(contigs_folder)
            logging.debug(f'Contigs directory created: {contigs_folder}')

        # Spades command
        contigs_command = spades.meta_contigs(unassigned_read_1, unassigned_read_2, contigs_folder)
        try:
            subprocess.run(contigs_command, check=True)
            logging.debug(f'Generating contigs successful: {contigs_command}')
        except subprocess.CalledProcessError as e:
            logging.debug(f'Generating contigs failed: {contigs_command}')

        # create interleaved files
        interleave_command = interleave.interleave_paired_reads(unassigned_read_1, unassigned_read_2,
                                                                os.path.join(output_directory, 'interleaved.fasta'))

        try:
            subprocess.run(interleave_command, check=True)
            logging.debug(f'Interleaving reads successful: {interleave_command}')
        except subprocess.CalledProcessError as e:
            logging.debug(f'Interleaving reads failed: {interleave_command}')

        contigs_file = os.path.join(contigs_folder, 'contigs.fasta')
        interleaved_file = os.path.join(output_directory, 'interleaved.fasta')
        if os.path.isfile(contigs_file) and os.path.isfile(interleaved_file):
            unassigned_local[directory[1]] = (contigs_file, interleaved_file)
        elif os.path.isfile(interleaved_file):
            unassigned_local[directory[1]] = (None, interleaved_file)
        else:
            unassigned_local[directory[1]] = (None, None)
    unassigned_dict[db] = unassigned_local

# Create diamond database
create_db_command = diamond.diamond_build_db(diamond_db_path, output_folder, 'diamond_db')
created_diamond_db = create_db_command[1]
try:
    subprocess.run(create_db_command[0], check=True)
    logging.debug(f'Create diamond database: {create_db_command[0]}')
except subprocess.CalledProcessError as e:
    logging.debug(f'Create diamond database failed: {create_db_command[0]}')

# Run diamond on unassigned reads
for db, unassigned_dict in unassigned_dict.items():
    db_folder = os.path.join(diamond_output, f'kraken_db{db+1}')
    if not os.path.exists(db_folder):
        os.makedirs(db_folder)
    for sample, unassigned_files in unassigned_dict.items():
        logging.debug(f'Running diamond on unassigned reads for: {sample}')
        output_directory = os.path.join(db_folder, sample)
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        if unassigned_files[0] is not None:
            diamond_command_contigs = diamond.diamond_classify(created_diamond_db,
                                                               unassigned_files[0],
                                                               output_directory,
                                                               f'{sample}_contigs')
            try:
                subprocess.run(diamond_command_contigs, check=True)
                logging.debug(f'Diamond command successful: {diamond_command_contigs}')
            except subprocess.CalledProcessError as e:
                logging.debug(f'Diamond command failed: {diamond_command_contigs}')
        if unassigned_files[1] is not None:
            diamond_command_interleaved = diamond.diamond_classify(created_diamond_db,
                                                                   unassigned_files[1],
                                                                   output_directory,
                                                                   f'{sample}_interleaved')
            try:
                subprocess.run(diamond_command_interleaved, check=True)
                logging.debug(f'Diamond command successful: {diamond_command_interleaved}')
            except subprocess.CalledProcessError as e:
                logging.debug(f'Diamond command failed: {diamond_command_interleaved}')

