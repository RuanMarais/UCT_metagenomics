"""
******************************************************************************************************************
kraken2_command_generation.py
Kraken2 command generation functions
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""

import os


def kraken2_command(input_read_1,
                    input_read_2,
                    kraken2_db,
                    threads,
                    output_directory,
                    report_prefix):
    """
    Generate kraken2 command to run using subprocess

    :param input_read_1: The filepath to the forward paired-end read
    :param input_read_2: The filepath to the reverse paired-end read
    :param kraken2_db: The filepath to the kraken2 database
    :param threads: The number of threads to use
    :param output_directory: The directory to output the kraken2 report and krakenfile
    :param report_prefix: The prefix to use for the kraken2 report and krakenfile
    :return: The kraken2 command to run using subprocess
    """
    report = os.path.join(output_directory, f'{report_prefix}.txt')
    krakenfile = os.path.join(output_directory, f'{report_prefix}_krakenfile')
    return ['kraken2', '--db', kraken2_db, '--threads', str(threads), '--report', report, '--paired', input_read_1,
            input_read_2, '--output', krakenfile]


def generate_kraken2_commands(reads_dict,
                              output_directory,
                              threads,
                              kraken2_db_path):
    """
    This function generates the list of kraken2 commands to run using subprocess

    :param reads_dict: The dictionary with each key referring to a tuple representing the forward and
    reverse paired-end reads
    :param output_directory: The filepath to the output directory
    :param threads: The number of threads to use
    :param kraken2_db_path: The filepath to the kraken2 database
    :return: The list of kraken2 commands to run using subprocess
    """
    commands = []
    for key, reads_tuple in reads_dict.items():
        participant_output = os.path.join(output_directory, key)
        if not os.path.exists(participant_output):
            os.makedirs(participant_output)
        command = kraken2_command(reads_tuple[0], reads_tuple[1], kraken2_db_path, threads, participant_output, key)
        commands.append(command)
    return commands
