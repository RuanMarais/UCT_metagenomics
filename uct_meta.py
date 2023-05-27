import argparse
import logging
import os
import data_and_variables as _dv
import meta_analysis as meta

# The input file for the pipeline
parser = argparse.ArgumentParser(description='This pipeline analyses metagenomic data. Created by GJK Marais.')

# Required parameters
parser.add_argument('--identifier_length', type=int, required=True,
                    help='The string length that identifies the sequencing read. Unique to each pair of reads')
parser.add_argument('--filename_length', type=int, required=True,
                    help='The string length that identifies the sequencing read. Up to the filetype suffix')

# Input/output
parser.add_argument('--input_folder', default=None,
                    help='The file path that contains Illumina paired-end reads')
parser.add_argument('--output_folder', default=None,
                    help='Result folders will be output to this folder')

# System
parser.add_argument('--threads', default=1,
                    help='The number of threads to use, default is 1')
parser.add_argument('--trimmomatic', default=None,
                    help='trimmomatic path')

# Database paths
parser.add_argument('--kraken2_db1', default=None,
                    help='kraken2 db1 path')
parser.add_argument('--kraken2_db2', default=None,
                    help='kraken2 db2 path')
parser.add_argument('--diamond_db', default=None,
                    help='diamond db path')
parser.add_argument('--kneaddata_db', default=None,
                    help='knead_data path')

# Input filetype data
parser.add_argument('--gen_id', default='fastq.gz', help='string used to identify all raw sequencing reads')
parser.add_argument('--r1_id', default='R1_001.fastq.gz', help='string used to identify read 1')
parser.add_argument('--r2_id', default='R2_001.fastq.gz', help='string used to identify read 2')

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

meta.meta_analysis(input_folder,
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
