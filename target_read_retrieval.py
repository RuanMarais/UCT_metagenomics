import argparse
import logging
import os
import services
import csv
import extract_from_krakenfile
import subprocess

# The input file for the pipeline
parser = argparse.ArgumentParser(description='This pipeline analyses metagenomic data. Created by GJK Marais.')
parser.add_argument('--input_folder_knead', required=True, help='The file path that contains the kneaddata paired '
                                                                'output files used for kraken analysis')
parser.add_argument('--input_folder_kraken2', required=True, help='Result folders for kraken2 results')
parser.add_argument('--output_folder', required=True, help='Result folders will be output to this folder')
parser.add_argument('--identifier_length', type=int, required=True,
                    help='The string length that identifies the sequencing read')
parser.add_argument('--reference_file', required=True, help='The csv file that refers to the taxid and reference to '
                                                            'use for each identifier')
args = parser.parse_args()

# Logging
logger = logging.getLogger(__name__)
handler = logging.FileHandler("Run.log")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.debug("Targeted read retrieval started")

# Set variables from input
input_folder_knead = args.input_folder_knead
output_folder = args.output_folder
id_length = args.identifier_length
reference_file = args.reference_file
kraken2_dir = args.input_folder_kraken2

# Generate a list with all valid filenames
file_names_all = services.filename_list_generate('001.fastq.gz', input_folder_knead)

# Generate a data dict which has the identifier as the key and the paired reads as a tuple
paired_sorted_dict = services.paired_read_list_generate(id_length, 'R1_001.fastq.gz',
                                                        'R2_001.fastq.gz', file_names_all)

# Generate krakenfile dict
kraken2_krakenfiles = {}
folders_kraken = os.listdir(kraken2_dir)
directory_paths_kraken = [(os.path.join(kraken2_dir, folder), folder)
                          for folder in folders_kraken
                          if os.path.isdir(os.path.join(kraken2_dir, folder))]
for directory in directory_paths_kraken:
    logging.debug(f'Kraken2 krakenfile retrieval: {directory}')
    # TODO: get correct path for kraken2 output reads
    filename_1 = f'{directory[1]}'
    file_1 = os.path.join(directory[0], filename_1)
    if os.path.isfile(file_1):
        logging.debug(f'Kraken2 krakenfile assigned to dictionary: {directory[1]}')
        kraken2_krakenfiles[directory[1]] = file_1


# Create an empty dictionary to store the data
data_dict_reference_file = {}

with open(reference_file, mode='r') as csvfile:
    logger.debug("Data retrieval from reference file started")
    # Read the CSV file using DictReader, which returns an iterator producing dictionaries
    try:
        csv_reader = csv.DictReader(csvfile)

        # Iterate through each row in the CSV file
        for row in csv_reader:
            # Get the values for 'identifier', 'taxid', 'reference' and 'reference_prefix'
            identifier = row['identifier']
            organism = row['organism']
            taxid = row['taxid']
            reference = row['reference']
            reference_prefix = row['reference_prefix']

            # Create a tuple with 'taxid' and 'reference' and store it as the value for 'identifier' in the dictionary
            data_dict_reference_file[identifier] = (organism, taxid, reference, reference_prefix)
    except:
        logger.debug("Data retrieval from reference file failed")

# Retrieve detected reads
extract_kraken_reads_commands = extract_from_krakenfile.extract_targeted_reads(kraken2_krakenfiles,
                                                                               data_dict_reference_file,
                                                                               paired_sorted_dict,
                                                                               output_folder)

# Path dict for alignment
alignment_dict = {}
for extraction_data in extract_kraken_reads_commands:
    alignment_dict[extraction_data[0]] = extraction_data[2]
    try:
        subprocess.run(extraction_data[1], check=True)
        logging.debug(f'krakenfile extraction command: {extraction_data[1]}')
    except subprocess.CalledProcessError as e:
        logging.debug(f'krakenfile extraction command failed: {extraction_data[1]}')


for command_dict in subprocess_commands_dict_list:
        try:
            subprocess.run(command_dict['BWA'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Candidate reference index for BWA: {command_dict['BWA']} with error code: {e.returncode}")
        try:
            subprocess.run(command_dict['toBAM'], check=True, stdout=open(command_dict['BAM'], 'wb'))
        except subprocess.CalledProcessError as e:
            error_log.append(f"Conversion to BAM: {command_dict['toBAM']} with error code: {e.returncode}")
        try:
            subprocess.run(command_dict['sortBAM'], check=True)
        except subprocess.CalledProcessError as e:
            error_log.append(f"Sorting BAM: {command_dict['sortBAM']} with error code: {e.returncode}")


