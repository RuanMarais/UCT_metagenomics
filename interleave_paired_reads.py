"""
******************************************************************************************************************
interleave_paired_reads.py
This function interleaves paired-end reads
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""


def interleave_paired_reads(input_read_1_path,
                            input_read_2_path,
                            output_path):
    """
    This function interleaves paired-end reads using PEAR

    :param input_read_1_path: The filepath to the forward paired-end read
    :param input_read_2_path: The filepath to the reverse paired-end read
    :param output_path: The filepath to save the interleaved reads
    :return: The PEAR command to run using subprocess
    """
    return ['pear', '-f', input_read_1_path, '-r', input_read_2_path, '-o', output_path]