"""
******************************************************************************************************************
reference_fragmentation.py
Generation of simulated reads that is run through the pipeline
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""
import os


def pirs_command(reference_path,
                 output_path,
                 coverage,
                 read_length,
                 fragment_length,
                 threads,
                 pirs_env=None):
    """
    This function generates the PIRS command for subprocess to produce simulated reads using PIRS

    :param reference_path: The filepath to the reference genome
    :param output_path: The filepath to save the simulated reads
    :param coverage: The coverage to simulate
    :param read_length: The read length to simulate
    :param fragment_length: The fragment length to simulate
    :param threads: The number of threads to use
    :param pirs_env: The name of the conda environment to run PIRS in
    :return: The PIRS command to run using subprocess
    """
    if pirs_env is None:
        return ['pirs', 'simulate', reference_path, '-x', str(coverage), '-l', str(read_length), '-m',
                str(fragment_length), '-o', output_path, '-t', str(threads), '-z']
    else:
        return ['conda', 'run', '-n', pirs_env, 'pirs', 'simulate', reference_path, '-x', str(coverage), '-l',
                str(read_length), '-m', str(fragment_length), '-o', output_path, '-t', str(threads), '-z']


def generate_pirs_commands(source,
                           pirs_env_name,
                           coverage,
                           read_length,
                           fragment_length,
                           threads):
    """
    This function generates the list PIRS commands and the output directories for subprocess to produce simulated
    reads using PIRS
    :param source: The filepath to the directory containing the reference genomes to simulate reads from
    :param pirs_env_name: The name of the conda environment to run PIRS in
    :param coverage: The coverage to simulate
    :param read_length: The read length to simulate
    :param fragment_length: The fragment length to simulate
    :param threads: The number of threads to use
    :return: The list of command lists to run using subprocess
    """
    commands = []
    # Get the reference genome filepaths
    files = os.listdir(source)
    directory_paths_knead = [(os.path.join(source, ref), ref)
                             for ref in files
                             if os.path.isfile(os.path.join(source, ref))]

    # Create the output directory for simulated reads
    output_dir = os.path.join(source, 'pirs_output')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate over the reference genomes and generate the PIRS command
    for directory_path, ref in directory_paths_knead:
        # Create the output directory for the simulated reads
        output_path = os.path.join(output_dir, ref)
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_final = os.path.join(output_path, ref)
        # Generate the PIRS command
        command = pirs_command(directory_path,
                               output_final,
                               coverage,
                               read_length,
                               fragment_length,
                               threads,
                               pirs_env_name)
        commands.append(command)
    return commands
