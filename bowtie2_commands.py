"""
******************************************************************************************************************
General helper functions for the pipeline
******************************************************************************************************************
"""
import os
import subprocess


def index_reference(reference_path, reference_prefix, location):
    ref_loc = os.path.join(location, reference_path)
    ref_out = os.path.join(location, reference_prefix)
    return ['bowtie2-build', ref_loc, ref_out]


def generate_bowtie2_command(organism, reference, read_1, read_2, output_dir, threads):
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


def run_alignment(command_dict, logger):
    # Retrieve commands from command dict
    organism = command_dict['organism']
    align_command = command_dict['align']
    view_command = command_dict['view']
    sort_command = command_dict['sort']
    index_command = command_dict['index']
    sorted_bamfile = command_dict['bamfile']

    try:
        subprocess.run(align_command, check=True)
    except subprocess.CalledProcessError:
        logger.error(f'Align command failed: {align_command}')

    try:
        p1 = subprocess.Popen(view_command, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(sort_command, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        p2.communicate()
    except subprocess.CalledProcessError:
        logger.error(f'Bam processing failed: {view_command}, {sort_command}')

    try:
        subprocess.run(index_command, check=True)
    except subprocess.CalledProcessError:
        logger.error(f'Bam index failed: {index_command}')

    return organism, sorted_bamfile
