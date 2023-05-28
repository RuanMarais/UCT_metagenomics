"""
******************************************************************************************************************
extract_from_krakenfile.py
This script extracts reads from a krakenfile matching a specific taxid using extract_kraken_reads.py
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""
import os


def extract_reads_command(krakenfile_path,
                          org,
                          taxid,
                          input_read_1_path,
                          input_read_2_path,
                          output_dir):
    """
    This function generates the command to run using subprocess to extract reads after kraken2
    classification based on taxid

    :param krakenfile_path: The filepath to the krakenfile
    :param org: The name of the organism being extracted for filename purposes
    :param taxid: The NCBI taxonomy ID of the organism
    :param input_read_1_path: The filepath to the forward read that was classified using kraken2
    :param input_read_2_path: The filepath to the reverse read that was classified using kraken2
    :param output_dir: The directory filepath to save the extracted reads
    :return: [0] The command to run using subprocess, [1] The output filepaths for the forward and reverse
    extracted reads, [2] The output directory filepath
    """
    org_output = os.path.join(output_dir, org)
    if not os.path.exists(org_output):
        os.makedirs(org_output)
    output_paths = (os.path.join(org_output, org + '_1.fastq'), os.path.join(org_output, org + '_2.fastq'))
    return ['extract_kraken_reads.py', '-k', krakenfile_path, '-t', taxid, '--fastq-output', '-s', input_read_1_path,
            '-s2', input_read_2_path, '-o', os.path.join(org_output, org + '_1.fastq'), '-o2',
            os.path.join(org_output, org + '_2.fastq')], output_paths, org_output


def extract_targeted_reads(krakenfile_dict_1,
                           krakenfile_dict_2,
                           target_dict,
                           raw_reads_dict,
                           output_dir,
                           reference_folder,
                           logger):
    """
    This function passes the variables from a target_dict object to extract_reads_command to create all
    commands necessary to run kraken2 taxid based read extraction

    :param krakenfile_dict_1: The dictionary with keys corresponding to the key in the target_dict referring to
    krakenfile filepaths. This is for krakenfiles generated with kraken2 database 1.
    :param krakenfile_dict_2: The dictionary with keys corresponding to the key in the target_dict referring to
    krakenfile filepaths. This is for krakenfiles generated with kraken2 database 2.
    :param target_dict:
    :param raw_reads_dict: The dictionary with keys corresponding to the key in the target_dict referring to
    the raw reads that were used for kraken2 classification. Reads are saved as a tuple (read_1 path, read_2 path)
    :param output_dir: The output directory filepath for the reads to be extracted to
    :param reference_folder: The directory where reference genomes are saved
    :param logger: A logging object for logfile creation
    """
    output_data = []
    for id_item, reference_data_items in target_dict.items():
        id_output = os.path.join(output_dir, id_item)
        if not os.path.exists(id_output):
            os.makedirs(id_output)
        for reference_data in reference_data_items:
            try:
                # Retrieve taxid from data file object
                if reference_data[1] != '':
                    taxid = str(reference_data[1])
                else:
                    logger.error(f'Error: {id_item} taxid is empty')
                    taxid = None

                # Retrieve kraken2 database used and the corresponding filepath to the generated krakenfile
                # for the classified reads
                if reference_data[4] != '':
                    kraken_dict = int(reference_data[4])
                    if kraken_dict == 1:
                        if os.path.isfile(krakenfile_dict_1[id_item]):
                            krakenfile = krakenfile_dict_1[id_item]
                        else:
                            logger.error(f'Error: {id_item} krakenfile does not exist')
                            krakenfile = None
                    elif kraken_dict == 2:
                        if os.path.isfile(krakenfile_dict_2[id_item]):
                            krakenfile = krakenfile_dict_2[id_item]
                        else:
                            logger.error(f'Error: {id_item} krakenfile does not exist')
                            krakenfile = None
                    else:
                        logger.error('Error: kraken2 database value is not 1 or 2')
                        krakenfile = None
                else:
                    logger.error(f'Error: {id_item} kraken2 database value is empty')
                    krakenfile = None

                # Retrieve the organism name
                if reference_data[0] != '':
                    organism = reference_data[0]
                else:
                    logger.error(f'Error: {id_item} organism name is empty')
                    organism = None

                # Retrieve the raw read filepaths that were classified using kraken2
                if os.path.isfile(raw_reads_dict[id_item][0]):
                    raw_read_1 = raw_reads_dict[id_item][0]
                else:
                    logger.error(f'Error: {id_item} raw read 1 file does not exist')
                    raw_read_1 = None
                if os.path.isfile(raw_reads_dict[id_item][1]):
                    raw_read_2 = raw_reads_dict[id_item][1]
                else:
                    logger.error(f'Error: {id_item} raw read 2 file does not exist')
                    raw_read_2 = None

                # Retrieve the reference filepath for the target
                if reference_data[3] != '':
                    reference_name = reference_data[3]
                    if os.path.isfile(os.path.join(reference_folder, reference_name)):
                        reference = os.path.join(reference_folder, reference_name)
                    else:
                        logger.error(f'Error: {id_item} reference file does not exist')
                        reference = None
                else:
                    logger.error(f'Error: {id_item} reference prefix is empty')
                    reference = None

                # Check if required variables are
                if krakenfile is not None and taxid is not None and raw_read_1 is not None and raw_read_2 is not None \
                        and reference is not None and organism is not None:

                    # Generate the read extraction commands and metadata for the generated variables
                    extract_data = extract_reads_command(krakenfile, organism, taxid, raw_read_1, raw_read_2, id_output)
                    # Retrieve the command to run with subprocess
                    command = extract_data[0]
                    # The filepaths that reads will be extracted to
                    output_paths = extract_data[1]
                    # The directory that the extracted reads will be saved to
                    out_dir = extract_data[2]
                    # The output array item
                    output = [id_item, command, output_paths, organism, reference, out_dir]
                    output_data.append(output)
            except IndexError:
                logger.error(f'Error: {id_item} reference_data is not the correct length (index out of range for '
                             f'object created from the specified file that contains the read extraction data)')
    return output_data


