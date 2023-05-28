"""
******************************************************************************************************************
bowtie2_commands.py
General helper functions for the pipeline related to genome alignment with Bowtie2
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""
import os
import subprocess


def index_reference(reference_path,
                    reference_prefix,
                    location):
    """
    This function indexes the reference genome using bowtie2-build

    :param reference_path: Path to the reference genome being indexed
    :param reference_prefix: The reference short name used for files and variable names
    :param location: The filepath for the output
    :return: The command to be run using subprocess
    """
    ref_loc = os.path.join(location, reference_path)
    ref_out = os.path.join(location, reference_prefix)
    return ['bowtie2-build', ref_loc, ref_out]


def generate_bowtie2_command(organism,
                             reference,
                             read_1,
                             read_2,
                             output_dir,
                             threads):
    """
    This function generates the commands bowtie2 alignment. The commands are returned as a dictionary

    :param organism: The name of the organism being aligned for filename purposes
    :param reference: The filepath to the reference used for the alignment
    :param read_1: The filepath to the forward read file
    :param read_2: The filepath to the reverse read file
    :param output_dir: The output directory filepath for alignment files
    :param threads: The number of threads to use for the alignment
    :return: A dictionary with the keys: organism, out, align, view, sort, index, bamfile that can be used to
    sequentially run alignment functions using subprocess
    """
    # Align reads with Bowtie2
    sam_file = os.path.join(output_dir, f'{organism}.sam')
    cmd_align = ['bowtie2', '-x', reference, '-1', read_1, '-2', read_2, '-S', sam_file, '-p', str(threads),
                 '--very-sensitive']

    # Convert SAM to BAM, sort and index with samtools
    bam_file = os.path.join(output_dir, f'{organism}.bam')
    cmd_view = ['samtools', 'view', '-S', '-b', sam_file]
    cmd_sort = ['samtools', 'sort', '-o', bam_file]
    cmd_index = ['samtools', 'index', bam_file]

    return {'organism': organism,
            'out': output_dir,
            'align': cmd_align,
            'view': cmd_view,
            'sort': cmd_sort,
            'index': cmd_index,
            'bamfile': bam_file}


def run_alignment(command_dict,
                  logger):
    """
    This function runs the alignment commands sequentially using subprocess.

    :param command_dict: A dictionary with commands produced using generate_bowtie2_command function. The keys are:
    organism, out, align, view, sort, index, bamfile
    :param logger: A logging object that is used for logfile generation initialised in the primary pipeline
    :return: [0] The organism name, [1] The filepath to the sorted BAM file
    """
    # Retrieve commands from command dict
    organism = command_dict['organism']
    align_command = command_dict['align']
    view_command = command_dict['view']
    sort_command = command_dict['sort']
    index_command = command_dict['index']
    sorted_bamfile = command_dict['bamfile']

    # Run alignment
    try:
        subprocess.run(align_command, check=True)
        logger.info(f'Aligned using bowtie2: {organism}')
    except subprocess.CalledProcessError:
        logger.error(f'Align command failed: {align_command}')

    # Convert SAM to BAM and sort
    try:
        p1 = subprocess.Popen(view_command, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(sort_command, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p2.communicate()
        logger.info(f'Bam processing completed: {organism}')
    except subprocess.CalledProcessError:
        logger.error(f'Bam processing failed: {view_command}, {sort_command}')

    # Index bamfile
    try:
        subprocess.run(index_command, check=True)
        logger.info(f'Bam index completed: {organism}')
    except subprocess.CalledProcessError:
        logger.error(f'Bam index failed: {index_command}')

    return organism, sorted_bamfile


def coverage_command(bamfile,
                     output_file):
    """
    This function generates the BAM file evaluation command

    :param bamfile: BAM file to be evaluated
    :param output_file: The output path for the evaluation file
    :return: The command to be run using subprocess
    """
    return ['samtools', 'coverage', bamfile, '-o', output_file]
