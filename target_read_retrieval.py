import argparse
import logging
import os
import services
import csv
import extract_from_krakenfile
import bowtie2_commands as bowtie2
import subprocess
import pysam
from collections import defaultdict
import numpy as np
import reference_fragmentation as frag

# The input file for the pipeline
parser = argparse.ArgumentParser(description='This pipeline analyses metagenomic data. Created by GJK Marais.')
parser.add_argument('--input_folder_knead', required=True, help='The file path that contains the kneaddata paired '
                                                                'output files used for kraken analysis')
parser.add_argument('--input_folder_kraken2_1', required=True, help='Result folder for kraken2 results from database 1')
parser.add_argument('--input_folder_kraken2_2', required=True, help='Result folder for kraken2 results from database 2')
parser.add_argument('--output_folder', required=True, help='Result folders will be output to this folder')
parser.add_argument('--identifier_length', type=int, required=True,
                    help='The string length that identifies the sequencing read')
parser.add_argument('--threads', type=int, default=1,
                    help='Number of threads to use for bowtie2 alignment')
parser.add_argument('--reference_file', required=True, help='The csv file that refers to the taxid and reference to '
                                                            'use for each identifier')
parser.add_argument('--references_directory', required=True, help='Directory for ncbi references')
args = parser.parse_args()

# Logging
logger = logging.getLogger(__name__)
handler = logging.FileHandler("Read retrieval.log")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
handler.setLevel(logging.DEBUG)
logger.setLevel(logging.DEBUG)
logger.debug("Targeted read retrieval started")

# Set variables from input
input_folder_knead = args.input_folder_knead
main_output = args.output_folder
output_folder = os.path.join(main_output, 'targeted_reads_retrieval')
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
id_length = args.identifier_length
reference_file = args.reference_file
kraken2_dir_1 = args.input_folder_kraken2_1
kraken2_dir_2 = args.input_folder_kraken2_2
references_directory = args.references_directory
threads = args.threads

# Generate a list with all valid filenames
file_names_all = services.filename_list_generate('.fastq', input_folder_knead)

# Generate a data dict which has the identifier as the key and the paired reads as a tuple
paired_sorted_dict = services.paired_read_list_generate(id_length,
                                                        'L001_R1_001_kneaddata_paired_1.fastq',
                                                        'L001_R1_001_kneaddata_paired_2.fastq',
                                                        file_names_all,
                                                        input_folder_knead)

# Generate krakenfile dict
kraken2_krakenfiles_1 = {}
folders_kraken = os.listdir(kraken2_dir_1)
directory_paths_kraken = [(os.path.join(kraken2_dir_1, folder), folder)
                          for folder in folders_kraken
                          if os.path.isdir(os.path.join(kraken2_dir_1, folder))]
for directory in directory_paths_kraken:
    logger.debug(f'Kraken2 krakenfile retrieval: {directory}')
    filename_1 = f'{directory[1]}_krakenfile'
    file_1 = os.path.join(directory[0], filename_1)
    if os.path.isfile(file_1):
        logger.debug(f'Kraken2 krakenfile assigned to dictionary: {directory[1]}')
        kraken2_krakenfiles_1[directory[1]] = file_1

# Generate krakenfile dict
kraken2_krakenfiles_2 = {}
folders_kraken = os.listdir(kraken2_dir_2)
directory_paths_kraken = [(os.path.join(kraken2_dir_2, folder), folder)
                          for folder in folders_kraken
                          if os.path.isdir(os.path.join(kraken2_dir_2, folder))]
for directory in directory_paths_kraken:
    logger.debug(f'Kraken2 krakenfile retrieval: {directory}')
    filename_1 = f'{directory[1]}_krakenfile'
    file_1 = os.path.join(directory[0], filename_1)
    if os.path.isfile(file_1):
        logger.debug(f'Kraken2 krakenfile assigned to dictionary: {directory[1]}')
        kraken2_krakenfiles_2[directory[1]] = file_1

# Create an empty dictionary to store the data
data_dict_reference_file = defaultdict(list)

with open(reference_file, mode='r', encoding='utf-8-sig') as csvfile:
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
            kraken_database = row['kraken_database']

            # Create a tuple with 'taxid' and 'reference' and store it as the value for 'identifier' in the dictionary
            data_dict_reference_file[identifier].append([organism, taxid, reference, reference_prefix, kraken_database])
    except:
        logger.debug("Data retrieval from reference file failed")

# Generate organism reference dict
organism_reference_dict = {}
for key, values in data_dict_reference_file.items():
    for value in values:
        organism_reference_dict[value[0]] = (value[2], value[3])


# Retrieve detected reads
extract_kraken_reads_commands = extract_from_krakenfile.extract_targeted_reads(kraken2_krakenfiles_1,
                                                                               kraken2_krakenfiles_2,
                                                                               data_dict_reference_file,
                                                                               paired_sorted_dict,
                                                                               output_folder,
                                                                               references_directory)

# Path dict for alignment
alignment_dict = defaultdict(list)
for extraction_data in extract_kraken_reads_commands:
    # alignment dictionary: key = identifier, value = [organism, output_paths, reference, out_dir]
    alignment_dict[extraction_data[0]].append([extraction_data[3],
                                               extraction_data[2],
                                               extraction_data[4],
                                               extraction_data[5]])
    try:
        subprocess.run(extraction_data[1], check=True)
        logger.debug(f'krakenfile extraction command: {extraction_data[1]}')
    except subprocess.CalledProcessError as e:
        logger.debug(f'krakenfile extraction command failed: {extraction_data[1]}')


# index references
for organism, reference_data in organism_reference_dict.items():
    command = bowtie2.index_reference(reference_data[0], reference_data[1], references_directory)
    try:
        subprocess.run(command, check=True)
        logger.debug(f'Indexing command: {command}')
    except subprocess.CalledProcessError as e:
        logger.debug(f'Indexing command failed: {command}')

    # Generate reference fragmented reads and kraken analysis
    command_reference_reads = frag.generate_pirs_commands()


# Generate alignment commands dictionary
command_dict_defaultdict = defaultdict(list)
for sample, alignment_data_items in alignment_dict.items():
    for alignment_data in alignment_data_items:
        alignment_commands = bowtie2.generate_bowtie2_command(alignment_data[0],
                                                              alignment_data[2],
                                                              alignment_data[1][0],
                                                              alignment_data[1][1],
                                                              alignment_data[3],
                                                              threads)
        command_dict_defaultdict[sample].append(alignment_commands)

align_results_folder = os.path.join(output_folder, 'align_results')
if not os.path.exists(align_results_folder):
    os.makedirs(align_results_folder)
for sample, command_dict_list in command_dict_defaultdict.items():
    sample_results_folder = os.path.join(align_results_folder, sample)
    if not os.path.exists(sample_results_folder):
        os.makedirs(sample_results_folder)
    # Run bowtie2 commands
    for command_dict in command_dict_list:
        organism = command_dict['organism']
        align_command = command_dict['align']
        view_command = command_dict['view']
        sort_command = command_dict['sort']
        index_command = command_dict['index']
        sorted_bamfile = command_dict['bamfile']
        out_file_name = f'{organism}_alignment_positions.txt'
        out_file = os.path.join(sample_results_folder, out_file_name)

        try:
            subprocess.run(align_command, check=True)
        except subprocess.CalledProcessError as e:
            logger.debug(f'Align command failed: {align_command}')

        try:
            p1 = subprocess.Popen(view_command, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(sort_command, stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
            p2.communicate()
        except subprocess.CalledProcessError as e:
            logger.debug(f'Bam processing failed: {view_command}, {sort_command}')

        try:
            subprocess.run(index_command, check=True)
        except subprocess.CalledProcessError as e:
            logger.debug(f'Bam index failed: {index_command}')

        bamfile = pysam.AlignmentFile(sorted_bamfile, 'rb')

        with open(out_file, 'w') as f:
            for read in bamfile:
                avg_phred_score = np.mean(read.query_qualities) if read.query_qualities is not None else 'N/A'
                f.write(f'Read: {read.query_name}, '
                        f'Start: {read.reference_start}, '
                        f'End: {read.reference_end}, '
                        f'Length: {read.query_length},'
                        f'Mapping Quality: {read.mapping_quality}, '
                        f'Average Phred Score: {avg_phred_score}\n')
        bamfile.close()
