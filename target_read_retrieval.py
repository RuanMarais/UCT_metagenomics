"""
******************************************************************************************************************
target_read_retrieval.py
This function retrieves and analyses reads for specific organisms from a kraken2 based metagenomic analysis.
The details of the analysis is passed to the function in a csv file
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""

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
import shutil
import meta_analysis as meta
import data_and_variables as _dv

# The input file for the pipeline
parser = argparse.ArgumentParser(description='This pipeline analyses metagenomic data. Created by GJK Marais.')

# Required arguments
parser.add_argument('--identifier_length', '-e', type=int, required=True,
                    help='The string length that identifies the sequencing read')
parser.add_argument('--reference_file', '-r', required=True,
                    help='The csv file that refers to the taxid and reference to use for each sample')
parser.add_argument('--references_directory', '-d', required=True,
                    help='Directory for ncbi references')

# Input/Output folders
parser.add_argument('--uct_meta_output', '-k', Required=True,
                    help='The file path that uct_meta.py outputs to')
parser.add_argument('--output_folder', '-o', default=None,
                    help='Result folders will be output to this folder')

# System parameters
parser.add_argument('--threads', '-t', type=int, default=1,
                    help='Number of threads to use for bowtie2 alignment')
parser.add_argument('--pirs_env', '-p', default=None,
                    help='conda environment for pirs')
parser.add_argument('--trimmomatic', '-m', default=None,
                    help='trimmomatic path')

# Database paths
parser.add_argument('--kraken2_db1', '-i', default=None,
                    help='kraken2 db1 path')
parser.add_argument('--kraken2_db2', '-I', default=None,
                    help='kraken2 db2 path')
parser.add_argument('--kneaddata_db', '-K', default=None,
                    help='knead_data path')
args = parser.parse_args()

# Set variables from input
id_length = args.identifier_length
file_length = id_length + 3
reference_file = args.reference_file
references_directory = args.references_directory
threads = args.threads
main_output = args.output_folder
if main_output is None:
    main_output = os.getcwd()
uct_meta_output = args.uct_meta_output

output_folder = os.path.join(main_output, 'targeted_reads_retrieval')
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
input_folder_knead = os.path.join(uct_meta_output, 'kraken2_source_files')
kraken2_dir_1 = os.path.join(uct_meta_output, 'kraken2_db1_results')
kraken2_dir_2 = os.path.join(uct_meta_output, 'kraken2_db2_results')
pirs_env = args.pirs_env
if pirs_env is None:
    pirs_env = _dv.pirs_env_name
trimmomatic_path = args.trimmomatic
if trimmomatic_path is None:
    trimmomatic_path = _dv.trimmomatic_path
kraken2_db1 = args.kraken2_db1
if kraken2_db1 is None:
    kraken2_db1 = _dv.kraken2_db1_path
kraken2_db2 = args.kraken2_db2
if kraken2_db2 is None:
    kraken2_db2 = _dv.kraken2_db2_path
kneaddata_db = args.kneaddata_db
if kneaddata_db is None:
    kneaddata_db = _dv.kneaddata_db_path

# Alignment output for reference evaluation
reference_alignment_output = os.path.join(main_output, 'reference_alignment_output')
if not os.path.exists(reference_alignment_output):
    os.makedirs(reference_alignment_output)

# Logging file generation
logger = logging.getLogger(__name__)
log_file = os.path.join(main_output, 'Read_retrieval.log')
handler = logging.FileHandler(log_file)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
handler.setLevel(logging.DEBUG)
logger.setLevel(logging.DEBUG)
logger.info("Targeted read retrieval started")

# Generate a list with all valid filenames
file_names_all = services.filename_list_generate('.fastq', input_folder_knead)

# Generate a data dict which has the identifier as the key and the paired reads as a tuple
paired_sorted_dict = services.paired_read_list_generate(id_length,
                                                        file_length,
                                                        '_kneaddata_paired_1.fastq',
                                                        '_kneaddata_paired_2.fastq',
                                                        file_names_all,
                                                        input_folder_knead,
                                                        logger)

# Generate krakenfile dict for kraken2 database 1
kraken2_krakenfiles_1 = {}
folders_kraken = os.listdir(kraken2_dir_1)
directory_paths_kraken = [(os.path.join(kraken2_dir_1, folder), folder)
                          for folder in folders_kraken
                          if os.path.isdir(os.path.join(kraken2_dir_1, folder))]
for directory in directory_paths_kraken:
    logger.info(f'Kraken2 krakenfile retrieval: {directory}')
    filename_1 = f'{directory[1]}_krakenfile'
    file_1 = os.path.join(directory[0], filename_1)
    if os.path.isfile(file_1):
        logger.info(f'Kraken2 krakenfile assigned to dictionary: {directory[1]}')
        kraken2_krakenfiles_1[directory[1]] = file_1

# Generate krakenfile dict for kraken2 database 1
kraken2_krakenfiles_2 = {}
folders_kraken = os.listdir(kraken2_dir_2)
directory_paths_kraken = [(os.path.join(kraken2_dir_2, folder), folder)
                          for folder in folders_kraken
                          if os.path.isdir(os.path.join(kraken2_dir_2, folder))]
for directory in directory_paths_kraken:
    logger.info(f'Kraken2 krakenfile retrieval: {directory}')
    filename_1 = f'{directory[1]}_krakenfile'
    file_1 = os.path.join(directory[0], filename_1)
    if os.path.isfile(file_1):
        logger.info(f'Kraken2 krakenfile assigned to dictionary: {directory[1]}')
        kraken2_krakenfiles_2[directory[1]] = file_1

# Create an empty dictionary to store the data from the provided reference_file
data_dict_reference_file = defaultdict(list)
reference_dict = {}
reference_prefix_dict = {}

with open(reference_file, mode='r', encoding='utf-8-sig') as csvfile:
    logger.info("Data retrieval from reference file started")
    # Read the CSV file using DictReader, which returns an iterator producing dictionaries
    try:
        csv_reader = csv.DictReader(csvfile)

        # Iterate through each row in the CSV file
        for row in csv_reader:
            # Get the values for 'identifier', 'taxid', 'reference', 'reference_prefix' and 'kraken_database'
            # for each row. This matches the reference_file passed to the function.
            identifier = f"{row['sample']}_R1"
            organism = row['organism']
            taxid = row['taxid']
            reference = row['reference']
            reference_prefix = row['reference_prefix']
            kraken_database = row['kraken_database']

            # Create tuples in dictionaries for further analysis.
            data_dict_reference_file[identifier].append([organism, taxid, reference, reference_prefix, kraken_database])
            reference_dict[reference] = (organism, reference_prefix, kraken_database)
            reference_prefix_dict[reference_prefix] = reference
    except KeyError:
        logger.error("Data retrieval from reference file failed: Ensure that the reference "
                     "file has the correct headers")
    except UnicodeDecodeError:
        logger.error("Data retrieval from reference file failed: Ensure that the reference "
                     "file is in the correct format: UTF-8 encoded csv")
    except FileNotFoundError:
        logger.error("Data retrieval from reference file failed: Ensure that the reference "
                     "file is in the correct location")
    except csv.Error:
        logger.error("Data retrieval from reference file failed: Ensure that the reference "
                     "file is in the correct format: UTF-8 encoded csv")



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
                                                                               references_directory,
                                                                               logger)

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
        logger.info(f'krakenfile extraction command successful')
    except subprocess.CalledProcessError:
        logger.debug(f'krakenfile extraction command failed: {extraction_data[1]}')

# Generate fragmented read folder
# This command is used to generate synthetic reads from the reference genomes to evaluate the percentage that is
# reference specific. The coverage is set at 10X with a read length of 100 and fragment length of 180
pirs_commands = frag.generate_pirs_commands(references_directory,
                                            pirs_env,
                                            _dv.coverage,
                                            _dv.read_length,
                                            _dv.fragment_length,
                                            threads)

# Run the pirs commands to generate raw synthetic reads
for command in pirs_commands:
    try:
        subprocess.run(command, check=True)
        logger.info(f'pIRS synthetic read generation')
    except subprocess.CalledProcessError as e:
        logger.error(f'pIRS synthetic read generation failed: {command}')

# Pirs data organisation object
fragmented_reads_dict = {}

# Check that the pirs output directory exists
pirs_dir = os.path.join(references_directory, 'pirs_output')
if os.path.isdir(pirs_dir):
    logger.info(f'pIRS output directory found: {pirs_dir}')
    files_pirs = os.listdir(pirs_dir)
    # iterate over the files in the pirs output directory which should contain the folders for the generated
    # raw synthetic reads and confirm they are directories and not equal to synthetic_reads which indicates repeat
    # analysis of the same data
    directory_paths_pirs = [(os.path.join(pirs_dir, ref), ref)
                            for ref in files_pirs
                            if os.path.isdir(os.path.join(pirs_dir, ref)) and ref != 'synthetic_reads']

    # Create a directory for the renamed synthetic reads
    source_dir_kraken = os.path.join(pirs_dir, 'synthetic_reads')
    if not os.path.exists(source_dir_kraken):
        os.makedirs(source_dir_kraken)

    # Check that raw synthetic reads were generated
    if len(directory_paths_pirs) > 0:

        # Iterate over the folders containing raw synthetic reads
        for directory in directory_paths_pirs:

            # Retrieve the information about the reference genome to populate the dict that will be used to run kraken2
            reference_org_name = reference_dict[directory[1]][1]
            kraken_db = reference_dict[directory[1]][2]
            organism_name = reference_dict[directory[1]][0]
            org_file = os.path.join(source_dir_kraken, reference_org_name)

            # Create an organism-specific directory for the renamed synthetic reads
            if not os.path.exists(org_file):
                os.makedirs(org_file)

            # Filepaths for the synthetic reads generated by Pirs
            read_1_target = os.path.join(directory[0],
                                         f'{directory[1]}_{_dv.read_length}_{_dv.fragment_length}_1.fq.gz')
            read_2_target = os.path.join(directory[0],
                                         f'{directory[1]}_{_dv.read_length}_{_dv.fragment_length}_2.fq.gz')

            # Check that the synthetic reads exist
            if os.path.isfile(read_1_target) and os.path.isfile(read_2_target):

                # New filenames
                read_1_out = os.path.join(org_file, f'{reference_org_name}_1.fastq.gz')
                read_2_out = os.path.join(org_file, f'{reference_org_name}_2.fastq.gz')

                # Move the synthetic reads to the organism-specific directory to limit data usage and facilitate ease
                # of access for manual and automated review and downstream analysis
                shutil.move(read_1_target, read_1_out)
                shutil.move(read_2_target, read_2_out)

                # Populate the fragmented reads dict which will be used to run kraken2
                fragmented_reads_dict[reference_org_name] = (org_file, organism_name, kraken_db)
                logger.info(f'Pirs synthetic reads moved: {directory[1]}')
            else:
                logger.error(f'Pirs raw synthetic reads not found: {directory[1]}')
    else:
        logger.error(f'Pirs output directory empty: {pirs_dir}')
else:
    logger.error(f'Pirs output directory not found: {pirs_dir}')

# Create a directory for the kraken2 output
alignment_dict_pirs = {}

# Iterate over the synthetic data dictionary to run kraken2
for organism_key, fragmented_reads_data in fragmented_reads_dict.items():

    # Retrieve data needed to run kraken2
    input_folder = fragmented_reads_data[0]
    organism_name = fragmented_reads_data[1]
    kraken_db = fragmented_reads_data[2]
    id_length = len(organism_key)
    file_length = id_length + 1
    reference = os.path.join(references_directory, reference_prefix_dict[organism_key])

    # file identifiers
    gen_id = 'fastq.gz'
    r1_id = '_1.fastq.gz'
    r2_id = '_2.fastq.gz'

    # Check which kraken2 database to use
    kraken2_db1_in = None
    kraken2_db2_in = None
    if int(kraken_db) == 1:
        kraken2_db1_in = kraken2_db1
    elif int(kraken_db) == 2:
        kraken2_db2_in = kraken2_db2

    # Run kraken2
    merged_reads = meta.meta_analysis(input_folder,
                                      id_length,
                                      input_folder,
                                      threads,
                                      trimmomatic_path,
                                      kneaddata_db,
                                      logger,
                                      gen_id,
                                      r1_id,
                                      r2_id,
                                      file_length,
                                      kraken2_db1_in,
                                      kraken2_db2_in,
                                      diamond_db_path=None,
                                      targeted=organism_key,
                                      target_name=organism_name,
                                      level=_dv.species)
    alignment_dict_pirs[organism_key] = (merged_reads, reference)

# index references
for organism, reference_data in organism_reference_dict.items():
    command = bowtie2.index_reference(reference_data[0], reference_data[1], references_directory)
    try:
        subprocess.run(command, check=True)
        logger.info(f'Indexing command: {organism}')
    except subprocess.CalledProcessError:
        logger.error(f'Indexing command failed: {command}')

# Generate alignment commands dictionary for reference synthetic read alignment
command_dict_defaultdict_ref = defaultdict(list)
for ref_name, alignment_data_items in alignment_dict_pirs.items():
    ref_dir = os.path.join(reference_alignment_output, ref_name)
    if not os.path.exists(ref_dir):
        os.makedirs(ref_dir)
    for alignment_data in alignment_data_items:
        alignment_commands = bowtie2.generate_bowtie2_command(ref_name,
                                                              alignment_data_items[1],
                                                              alignment_data_items[0][0],
                                                              alignment_data_items[0][1],
                                                              ref_dir,
                                                              threads)
        command_dict_defaultdict_ref[ref_name].append(alignment_commands)

# Run alignment
reference_folder_output = os.path.join(output_folder, 'reference_evaluation')
if not os.path.exists(reference_folder_output):
    os.makedirs(reference_folder_output)
for ref_name, command_dict in command_dict_defaultdict_ref.items():
    align_output = bowtie2.run_alignment(command_dict, logger)
    output_file = os.path.join(reference_folder_output, f'{ref_name}')
    bowtie2.coverage_command(align_output[1], output_file)

# Generate alignment commands dictionary for targeted retrieval
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
        align_output = bowtie2.run_alignment(command_dict, logger)
        organism = align_output[0]
        sorted_bamfile = align_output[1]
        out_file_name = f'{organism}_alignment_positions.txt'
        out_file = os.path.join(sample_results_folder, out_file_name)
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
