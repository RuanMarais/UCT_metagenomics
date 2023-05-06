"""
******************************************************************************************************************
This this script interleaves paired-end reads
******************************************************************************************************************
"""


def interleave_paired_reads(input_read_1_path, input_read_2_path, output_path):
    return ['pear', '-f', input_read_1_path, '-r', input_read_2_path, '-o', output_path]