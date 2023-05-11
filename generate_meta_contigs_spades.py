"""
******************************************************************************************************************
This this script generates contigs using SPAdes from metagenomic data
******************************************************************************************************************
"""


def meta_contigs(input_read_1_path, input_read_2_path, output_dir):
    return ['spades.py', '--meta', '-1', input_read_1_path, '-2', input_read_2_path, '-o',
            output_dir]
