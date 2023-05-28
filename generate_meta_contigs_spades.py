"""
******************************************************************************************************************
generate_meta_contigs_spades.py
This script generates contigs using SPAdes from metagenomic data
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""


def meta_contigs(input_read_1_path, input_read_2_path, output_dir):
    """

    :param input_read_1_path:
    :param input_read_2_path:
    :param output_dir:
    :return:
    """
    return ['spades.py', '--meta', '-1', input_read_1_path, '-2', input_read_2_path, '-o',
            output_dir]
