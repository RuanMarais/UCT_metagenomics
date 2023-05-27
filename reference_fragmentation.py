"""
******************************************************************************************************************
Generation of simulated reads that is run through the pipeline
******************************************************************************************************************
"""
import os


def pirs_command(reference_path, output_path, coverage, read_length, fragment_length, threads, pirs_env=None):
    if pirs_env is None:
        return ['pirs', 'simulate', reference_path, '-x', str(coverage), '-l', str(read_length), '-m',
                str(fragment_length), '-o', output_path, '-t', str(threads), '-z']
    else:
        return ['conda', 'run', '-n', pirs_env, 'pirs', 'simulate', reference_path, '-x', str(coverage), '-l',
                str(read_length), '-m', str(fragment_length), '-o', output_path, '-t', str(threads), '-z']


def generate_pirs_commands(source, pirs_env_name, coverage, read_length, fragment_length, threads):
    commands = []
    files = os.listdir(source)
    directory_paths_knead = [(os.path.join(source, ref), ref)
                             for ref in files
                             if os.path.isfile(os.path.join(source, ref))]
    output_dir = os.path.join(source, 'pirs_output')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for directory_path, ref in directory_paths_knead:
        output_path = os.path.join(output_dir, ref)
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        output_final = os.path.join(output_path, ref)
        command = pirs_command(directory_path,
                               output_final,
                               coverage,
                               read_length,
                               fragment_length,
                               threads,
                               pirs_env_name)
        commands.append(command)
    return commands
