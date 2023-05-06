"""
******************************************************************************************************************
This this script uses generates kneaddata commands for read QC
******************************************************************************************************************
"""
import os


def kneaddata_command_generate(input_read_1, input_read_2, reference_db, output_directory, threads, trimmomatic_path):
    return ['kneaddata', '--input', input_read_1, '--input', input_read_2, '--reference-db', reference_db,
            '--output', output_directory, '--threads', threads, '--trimmomatic', trimmomatic_path]


def generate_kneaddata_commands(raw_reads_dict, output_directory, threads, trimmomatic_path, kneaddata_db_path):
    commands = []
    for key, reads_tuple in raw_reads_dict.items():
        participant_output = os.path.join(output_directory, key)
        if not os.path.exists(participant_output):
            os.makedirs(participant_output)
        command = kneaddata_command_generate(reads_tuple[0], reads_tuple[1], kneaddata_db_path,
                                             participant_output, threads, trimmomatic_path)
        commands.append(command)
    return commands
