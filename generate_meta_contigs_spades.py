"""
******************************************************************************************************************
generate_meta_contigs_spades.py
This function generates contigs using SPAdes from metagenomic data
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""


def meta_contigs(input_read_1_path,
                 input_read_2_path,
                 output_dir):
    """
    This function generates contigs from 2 paired-end metagenomic reads using SPAdes

    :param input_read_1_path: filepath to the forward paired-end read
    :param input_read_2_path: filepath to the reverse paired-end read
    :param output_dir: The directory filepath to save the contigs
    :return: The SPAdes command to run using subprocess
    """
    return ['spades.py', '--meta', '-1', input_read_1_path, '-2', input_read_2_path, '-o',
            output_dir]
