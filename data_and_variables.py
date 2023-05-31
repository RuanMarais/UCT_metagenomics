"""
******************************************************************************************************************
data_and_variable.py
Data and variables for the UCT_metagenomics pipeline
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""
# This file is used to store variable for the pipeline to streamline its use when running multiple
# times with similar variables

# z-score analysis variables
z_score_threshold = 2
minimum_reads = 2
minimum_reads_genus = 2
flag_species = ['Human immunodeficiency virus 1']
flag_genus = ['Mycobacterium']
exclude_species = ['Homo sapiens']
exclude_genus = ['Homo']

# Database locations
kraken2_db1_path = None
kraken2_db2_path = None
kneaddata_db_path = None
diamond_db_path = None

# Package paths
trimmomatic_path = None

# kraken2 Rank codes
species = 'S'
genus = 'G'

# Pirs parameters
pirs_env_name = None
coverage = 10
read_length = 100
fragment_length = 180

# kneaddata read count retrieval strings
total_reads = 'READ COUNT: raw pair1 : Initial number of reads'
final_reads = 'READ COUNT: final pair1 : Total reads after merging results from multiple databases'
kneaddata_filesize_trf_threshold_mb = 20