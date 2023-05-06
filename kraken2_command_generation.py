"""
******************************************************************************************************************
Kraken2 command generation
******************************************************************************************************************
"""

import os


def kraken2_command(input_read_1, input_read_2, kraken2_db, threads, output_directory, report_prefix):
    return ['kraken2', '--db', kraken2_db, '--threads', threads, '--report',
            os.path.join(output_directory, f'{report_prefix}.txt'), '--paired', input_read_1, input_read_2]


def generate_kraken2_commands(reads_dict, output_directory, threads, kraken2_db_path):
    commands = []
    for key, reads_tuple in reads_dict.items():
        participant_output = os.path.join(output_directory, key)
        if not os.path.exists(participant_output):
            os.makedirs(participant_output)
        command = kraken2_command(reads_tuple[0], reads_tuple[1], kraken2_db_path, threads, participant_output, key)
        commands.append(command)
    return commands
