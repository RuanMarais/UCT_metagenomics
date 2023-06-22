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

# The species to always output for pathologist review (e.g. HIV, TB) and the species to exclude (e.g. human)
flag_species = ['Human immunodeficiency virus 1', 'Mycobacterium tuberculosis']
flag_genus = ['Mycobacterium']
exclude_species = ['Homo sapiens']
exclude_genus = ['Homo']

# Database paths
kraken2_db1_path = None
kraken2_db2_path = None
kneaddata_db_path = None
diamond_db_path = None

# Package paths
trimmomatic_path = None

# kraken2 Rank codes from the kraken 2 reports
species = 'S'
genus = 'G'

# pirs parameters
pirs_env_name = None
coverage = 10
read_length = 100
fragment_length = 180

# kneaddata read count retrieval strings and filesize threshold for running trf
total_reads = 'READ COUNT: raw pair1 : Initial number of reads'
final_reads = 'READ COUNT: final pair1 : Total reads after merging results from multiple databases'
kneaddata_filesize_trf_threshold_mb = 20