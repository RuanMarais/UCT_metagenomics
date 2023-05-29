"""
******************************************************************************************************************
uct_meta.py
This is the main script for the UCT_metagenomics pipeline.
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""

import argparse
import logging
import os
import data_and_variables as _dv
import meta_analysis as meta
import parse_kraken2_report as parse
from collections import defaultdict
import zscore_analyse as zscore

# The input file for the pipeline
parser = argparse.ArgumentParser(description='This pipeline analyses metagenomic data. Created by GJK Marais.')

# Required parameters
parser.add_argument('--identifier_length', '-e', type=int, required=True,
                    help='The string length that identifies the sequencing read. Unique to each pair of reads')
parser.add_argument('--filename_length', '-l', type=int, required=True,
                    help='The string length that identifies the sequencing read. Up to the filetype suffix')

# Input/output
parser.add_argument('--input_folder', '-i', default=None,
                    help='The file path that contains Illumina paired-end reads')
parser.add_argument('--output_folder', '-o', default=None,
                    help='Result folders will be output to this folder')

# System
parser.add_argument('--threads', '-t', default=1,
                    help='The number of threads to use, default is 1')
parser.add_argument('--trimmomatic', '-m', default=None,
                    help='trimmomatic path')

# Database paths
parser.add_argument('--kraken2_db1', '-1', default=None,
                    help='kraken2 db1 path')
parser.add_argument('--kraken2_db2', '-2', default=None,
                    help='kraken2 db2 path')
parser.add_argument('--diamond_db', '-d', default=None,
                    help='diamond db path')
parser.add_argument('--kneaddata_db', '-k', default=None,
                    help='knead_data path')

# Input filetype data
parser.add_argument('--gen_id', '-g', default='fastq.gz',
                    help='string used to identify all raw sequencing reads')
parser.add_argument('--r1_id', '-f', default='_R1.fastq.gz',
                    help='string used to identify read 1')
parser.add_argument('--r2_id', '-r', default='_R2.fastq.gz',
                    help='string used to identify read 2')

# metadata file
parser.add_argument('--metadata', '-a', required=True,
                    help='metadata file')
args = parser.parse_args()

# Set variables from input
input_folder = args.input_folder
if input_folder is None:
    input_folder = os.getcwd()
output_folder = args.output_folder
if output_folder is None:
    output_folder = os.getcwd()
id_length = args.identifier_length
threads = args.threads
trimmomatic_path = args.trimmomatic
kraken2_db1_path = args.kraken2_db1
kraken2_db2_path = args.kraken2_db2
diamond_db_path = args.diamond_db
kneaddata_db_path = args.kneaddata_db
gen_id = args.gen_id
r1_id = args.r1_id
r2_id = args.r2_id
file_length = args.filename_length
metadata = args.metadata

# Logging
logger = logging.getLogger(__name__)
handler = logging.FileHandler(os.path.join(output_folder, "UCT_metagenomics.log"))
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
logger.info("Metagenomics pipeline started")

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

# Run the metagenomic analysis pipeline
meta_output = meta.meta_analysis(input_folder,
                                 id_length,
                                 output_folder,
                                 threads,
                                 trimmomatic_path,
                                 kneaddata_db_path,
                                 logger,
                                 gen_id,
                                 r1_id,
                                 r2_id,
                                 file_length,
                                 kraken2_db1_path,
                                 kraken2_db2_path,
                                 diamond_db_path)

sample_dict = {}
# Import metadata
metadata_dict = zscore.importer(metadata)
run_list = []
for val in metadata_dict:
    run = val['run']
    if run not in run_list:
        run_list.append(run)
    sample = val['sample']
    sample_dict[sample] = {}
    genus_dict = {}
    species_dict = {}
    result_dict = {}
    sample_dict[sample]['metadata'] = val
    try:
        sample_dict[sample]['metadata']['nonhuman_reads'] = meta_output[1][sample][1]
        sample_dict[sample]['metadata']['total_reads'] = meta_output[1][sample][0]
    except KeyError:
        logger.info(f'No nonhuman reads data for {sample}')
    sample_dict[sample]['genus_dict'] = genus_dict
    sample_dict[sample]['species_dict'] = species_dict
    sample_dict[sample]['results'] = result_dict

# all species and genus
species_dict = defaultdict(list)
genus_dict = defaultdict(list)

for output_directory in meta_output[0]:
    # Get the sample data from the kraken output directory
    folders_kraken = os.listdir(output_directory)
    directory_paths_kraken = [(os.path.join(output_directory, folder), folder)
                              for folder in folders_kraken
                              if os.path.isdir(os.path.join(output_directory, folder))]
    for sample_directory in directory_paths_kraken:
        kraken_report = os.path.join(sample_directory[0], f'{sample_directory[1]}.txt')
        if os.path.isfile(kraken_report):
            parsed_kraken_data = parse.generate_kraken2_report_csv(kraken_report, sample_directory[0], logger)
            if parsed_kraken_data is not None:
                if parsed_kraken_data:
                    for data in parsed_kraken_data:
                        if int(data['Number of reads covered']) > 0 and data['Rank code'] == _dv.species:
                            species_dict[data['Scientific name']].append((sample_directory[1],
                                                                          int(data['Number of reads covered'])))
                        if int(data['Number of reads covered']) > 0 and data['Rank code'] == _dv.genus:
                            genus_dict[data['Scientific name']].append((sample_directory[1],
                                                                        int(data['Number of reads covered'])))
                else:
                    logger.warning(f'Kraken2 report empty {sample_directory[1]}')
            else:
                logger.warning(f'Kraken2 report parse failed for {sample_directory[1]}')
        else:
            logger.warning(f'No kraken2 report found for {sample_directory[1]}')

# Organise species and genus data
zscore.process_read_data(species_dict, sample_dict, logger)
zscore.process_read_data(genus_dict, sample_dict, logger)

# z-score analysis
zscore_output = os.path.join(output_folder, 'zscore_results')
if not os.path.isdir(zscore_output):
    os.mkdir(zscore_output)
zscore.z_score_analysis(sample_dict, run_list, zscore_output)
