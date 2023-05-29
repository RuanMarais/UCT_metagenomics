"""
******************************************************************************************************************
knead_data_run.py
These functions generate kneaddata commands for read quality control and human read depletion
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""
import os


def kneaddata_command_generate(input_read_1,
                               input_read_2,
                               reference_db,
                               output_directory,
                               threads,
                               trimmomatic_path=None):
    """
    This function generates a kneaddata command to run using subprocess

    :param input_read_1: The filepath to the forward paired-end read
    :param input_read_2: The filepath to the reverse paired-end read
    :param reference_db: The filepath to the reference database
    :param output_directory: The directory filepath to save the kneaddata output
    :param threads: The number of threads to use
    :param trimmomatic_path: The filepath to the trimmomatic jar file
    :return: The kneaddata command to run using subprocess
    """
    if trimmomatic_path is None:
        if 'NC' in input_read_1:
            return ['kneaddata', '--input', input_read_1, '--input', input_read_2, '--reference-db', reference_db,
                    '--output', output_directory, '--threads', str(threads), '--bypass-trf']
        else:
            return ['kneaddata', '--input', input_read_1, '--input', input_read_2, '--reference-db', reference_db,
                    '--output', output_directory, '--threads', str(threads)]
    else:
        if 'NC' in input_read_1:
            return ['kneaddata', '--input', input_read_1, '--input', input_read_2, '--reference-db', reference_db,
                    '--output', output_directory, '--threads', str(threads), '--trimmomatic', trimmomatic_path,
                    '--bypass-trf']
        else:
            return ['kneaddata', '--input', input_read_1, '--input', input_read_2, '--reference-db', reference_db,
                    '--output', output_directory, '--threads', str(threads), '--trimmomatic', trimmomatic_path]


def generate_kneaddata_commands(raw_reads_dict,
                                output_directory,
                                threads,
                                trimmomatic_path,
                                kneaddata_db_path):
    """
    This function generates a list of kneaddata commands to run using subprocess

    :param raw_reads_dict: The dictionary of raw reads. Each item is a tuple of the forward and reverse read
    :param output_directory: The directory filepath to save the kneaddata output
    :param threads: The number of threads to use
    :param trimmomatic_path: The filepath to the trimmomatic jar file
    :param kneaddata_db_path: The filepath to the kneaddata database
    :return:
    """
    commands = []
    for key, reads_tuple in raw_reads_dict.items():
        participant_output = os.path.join(output_directory, key)
        if not os.path.exists(participant_output):
            os.makedirs(participant_output)
        command = kneaddata_command_generate(reads_tuple[0], reads_tuple[1], kneaddata_db_path,
                                             participant_output, threads, trimmomatic_path)
        commands.append(command)
    return commands
