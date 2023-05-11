"""
******************************************************************************************************************
General helper functions for the pipeline
******************************************************************************************************************
"""


def index_reference(reference_path, reference_prefix):
    return ['bowtie2-build', reference_path, reference_prefix]


def generate_bowtie2_commands(alignment_files_dict, reference_dict, output_dir, threads):
    pass
    # bowtie2 - x reference_index - 1 reads_1.fastq - 2 reads_2.fastq - S output.sam - p 4
    #
    # convert_to_bam = ['samtools', 'view', '-b', '-S', output_sam]
    # output_bam_file = os.path.join(reference_output_directory, f'{assembly[0]}.bam')
    # sort_bam = ['samtools', 'sort', output_bam_file, '-o', sorted_output_bam_file]