import argparse
import logging
import os
import services
import data_and_variables as _dv
import knead_data_run
import subprocess
import shutil
import kraken2_command_generation as kraken2

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
handler = logging.FileHandler(os.path.join(output_folder, "Run.log"))
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
paired_sorted_dict = services.paired_read_list_generate(id_length, 'R1_001.fastq.gz',
                                                        'R2_001.fastq.gz', file_names_all)

# Generate output directories
knead_directory = os.path.join(output_folder, 'knead_data_results')
kraken_source_directory = os.path.join(output_folder, 'kraken2_source_files')
kraken2_directory_db1 = os.path.join(output_folder, 'kraken2_db1_results')
kraken2_directory_db2 = os.path.join(output_folder, 'kraken2_db2_results')

directories_to_create = [knead_directory, kraken2_directory_db1, kraken2_directory_db2, kraken_source_directory]

for directory in directories_to_create:
    if not os.path.exists(directory):
        os.makedirs(directory)
        logger.debug(f'Output directory created: {directory}')

# Generate kneaddata commands
kneaddata_commands = knead_data_run.generate_kneaddata_commands(paired_sorted_dict, knead_directory, threads,
                                                                trimmomatic_path, kneaddata_db_path)

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
    # TODO: get correct path for knead output reads
    filename_1 = f'{directory[1]}_paired_1.fasta'
    filename_2 = f'{directory[1]}_paired_2.fasta'
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


krakenfiles = {}
krakenresults_folders = [kraken2_directory_db1, kraken2_directory_db2]
for db, krakenfolder in enumerate(krakenresults_folders):
    folders_results = os.listdir(krakenfolder)
    paths_results = [(os.path.join(krakenfolder, folder), folder)
                     for folder in folders_results
                     if os.path.isdir(os.path.join(krakenfolder, folder))]
    for directory in paths_results:
        logging.debug(f'Kraken2 output read output retrieval: {directory}')
        # TODO: get correct path for knead output reads
        filename_1 = f'{directory[1]}_paired_1.fasta'
        filename_2 = f'{directory[1]}_paired_2.fasta'
        file_1 = os.path.join(directory[0], filename_1)
        file_2 = os.path.join(directory[0], filename_2)
        copy_path_1 = os.path.join(kraken_source_directory, filename_1)
        copy_path_2 = os.path.join(kraken_source_directory, filename_2)
        if os.path.isfile(file_1) and os.path.isfile(file_2):
            logging.debug(f'Kraken2 source files retrieved from: {directory}')
            shutil.copy(file_1, kraken_source_directory)
            shutil.copy(file_2, kraken_source_directory)
        kraken2_source_files[directory[1]] = (copy_path_1, copy_path_2)