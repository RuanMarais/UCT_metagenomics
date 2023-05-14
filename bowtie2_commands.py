"""
******************************************************************************************************************
General helper functions for the pipeline
******************************************************************************************************************
"""
import os


def index_reference(reference_path, reference_prefix):
    return ['bowtie2-build', reference_path, reference_prefix]


def generate_bowtie2_command(organism, reference, read_1, read_2, output_dir, threads):
    # Align reads with Bowtie2
    sam_file = os.path.join(output_dir, f'{organism}.sam')
    cmd_align = ['bowtie2', '-x', reference, '-1', read_1, '-2', read_2, '-S', sam_file, '-p', str(threads)]

    # Convert SAM to BAM, sort and index with samtools
    bam_file = os.path.join(output_dir, f'{organism}.bam')
    cmd_view = ['samtools', 'view', '-S', '-b', sam_file]
    cmd_sort = ['samtools', 'sort', '-o', bam_file]
    cmd_index = ['samtools', 'index', bam_file]

    return {'out': output_dir,
            'align': cmd_align,
            'view': cmd_view,
            'sort': cmd_sort,
            'index': cmd_index,
            'bamfile': bam_file}