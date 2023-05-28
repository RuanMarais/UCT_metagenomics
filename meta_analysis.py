"""
******************************************************************************************************************
meta_analysis.py
This is the primary metagenomic analysis file that organises read data, runs the read primary processing
for QC and performs read classification.
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""

import os
import services
import knead_data_run
import subprocess
import shutil
import kraken2_command_generation as kraken2
import extract_from_krakenfile as extract
import generate_meta_contigs_spades as spades
import interleave_paired_reads as interleave
import diamond_assign as diamond
import parse_kraken2_report as parse


def extract_reads(taxid,
                  filename,
                  krakenfiles,
                  logger,
                  extract_db_list=None,
                  output=None):
    """

    :param taxid: The NCBI Taxonomy ID to extract reads for
    :param filename: The filename to use for extracted reads
    :param krakenfiles: The dictionary of krakenfiles to extract reads from. The outer level has keys corresponding to
    kraken2 databases. The inner level has keys corresponding to the sample ID.
    :param logger: The logging object used to create the logfile
    :param extract_db_list: The dictionary corresponding to the kraken2 databases used for extraction. The keys refer
    to the kraken2 database and the item refers to the output folder for kraken2 results.
    :param output: The output directory to write the extracted reads to. If None then the directory specified in the
    extract_db_list dictionary is used.
    """
    for db, krakenfiles_dict in krakenfiles.items():
        logger.info(f'Extracting {filename} reads for db{db}')
        for sample, krakenfile in krakenfiles_dict.items():
            logger.info(f'Extracting {filename} reads for: {sample}')
            if output is None:
                output_directory = os.path.join(extract_db_list[db-1], sample)
            else:
                output_directory = os.path.join(output, sample)
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
                logger.info(f'Output directory created: {output_directory}')
            extract_reads_command = extract.extract_reads_command(krakenfile[0],
                                                                  filename,
                                                                  taxid,
                                                                  krakenfile[1],
                                                                  krakenfile[2],
                                                                  output_directory)
            try:
                subprocess.run(extract_reads_command[0], check=True)
                logger.info(f'Extracting {filename} reads successful')
            except subprocess.CalledProcessError:
                logger.error(f'Extracting {filename} reads failed: {extract_reads_command}')


def meta_analysis(input_folder,
                  id_length,
                  output_folder,
                  threads,
                  trimmomatic_path,
                  kneaddata_db_path,
                  logger,
                  general_filepath_id_string,
                  r1_file_id,
                  r2_file_id,
                  file_length,
                  kraken2_db1_path=None,
                  kraken2_db2_path=None,
                  diamond_db_path=None,
                  targeted=None,
                  target_name=None,
                  level=None):
    """
    This function performs the primary analysis of the metagenomic reads. It is used by both the uct_meta.py file and
    the target_read_retrieval.py file. For the former it takes a directory of raw data, performs QC and classification
    and used DIAMOND to further analyse unassigned reads. For the latter it takes a single pair of synthetic reads
    generated from a reference, performs metagenomic analysis and retrieves the reads that match the reference
    specifically. Only the latter function returns reads. These correspond to the reads that match the reference as
    merged forward and reverse read files.

    :param input_folder: The folder containing the raw read files
    :param id_length: The string length of in the filename that is unique to the sample's pair of reads
    :param output_folder: The folder to write the output to
    :param threads: The number of threads to use
    :param trimmomatic_path: The path to the trimmomatic jar file
    :param kneaddata_db_path: The path to the kneaddata database
    :param logger: The logging object used to create the logfile
    :param general_filepath_id_string: The string in the filename that is common to all reads
    :param r1_file_id: The string in the filename that is unique to all forward reads
    :param r2_file_id: The string in the filename that is unique to all reverse reads
    :param file_length: The string length of the filename excluding the file extension
    :param kraken2_db1_path: The path to the first kraken2 database
    :param kraken2_db2_path: The path to the second kraken2 database
    :param diamond_db_path: The path to the fasta file to use to create the DIAMOND database
    :param targeted: The string specifying the 'reference_prefix' from the reference_file from
    the target_read_retrieval function
    :param target_name: The string specifying the 'target_name' from the reference_file
    :param level: The taxonomic level to use for classification: species, genus
    :return: Merged reads for the specified target
    """
    # Generate a list with all valid filenames
    file_names_all = services.filename_list_generate(general_filepath_id_string, input_folder)

    # Separate filenames by id and into a list of pair tuples
    paired_sorted_dict = services.paired_read_list_generate(id_length,
                                                            file_length,
                                                            r1_file_id,
                                                            r2_file_id,
                                                            file_names_all,
                                                            input_folder,
                                                            logger)

    # Generate output directories
    knead_directory = os.path.join(output_folder, 'knead_data_results')
    kraken_source_directory = os.path.join(output_folder, 'kraken2_source_files')
    directories_to_create = [knead_directory,
                             kraken_source_directory]
    extract_db_list = [None, None]
    output_directories = []
    if kraken2_db1_path is not None:
        kraken2_directory_db1 = os.path.join(output_folder, 'kraken2_db1_results')
        output_directories.append(kraken2_directory_db1)
        extracted_reads_db1_directory = os.path.join(output_folder, 'extracted_reads_db1')
        directories_to_create.extend([kraken2_directory_db1, extracted_reads_db1_directory])
        extract_db_list[0] = extracted_reads_db1_directory
    if kraken2_db2_path is not None:
        kraken2_directory_db2 = os.path.join(output_folder, 'kraken2_db2_results')
        output_directories.append(kraken2_directory_db2)
        extracted_reads_db2_directory = os.path.join(output_folder, 'extracted_reads_db2')
        directories_to_create.extend([kraken2_directory_db2, extracted_reads_db2_directory])
        extract_db_list[1] = extracted_reads_db2_directory
    if diamond_db_path is not None:
        diamond_output = os.path.join(output_folder, 'diamond_output')
        directories_to_create.append(diamond_output)

    for directory in directories_to_create:
        if not os.path.exists(directory):
            os.makedirs(directory)
            logger.info(f'Output directory created: {directory}')

    # Generate kneaddata commands
    kneaddata_commands = knead_data_run.generate_kneaddata_commands(paired_sorted_dict,
                                                                    knead_directory,
                                                                    str(threads),
                                                                    trimmomatic_path,
                                                                    kneaddata_db_path)

    # Run kneaddata commands
    for knead_command in kneaddata_commands:
        try:
            subprocess.run(knead_command, check=True)
            logger.info(f'Kneaddata command run successful')
        except subprocess.CalledProcessError as e:
            logger.error(f'Kneaddata command failed: {knead_command}')

    # Generate kraken2 source directory
    kraken2_source_files = {}
    folders_knead = os.listdir(knead_directory)
    directory_paths_knead = [(os.path.join(knead_directory, folder), folder)
                             for folder in folders_knead
                             if os.path.isdir(os.path.join(knead_directory, folder))]
    for directory in directory_paths_knead:
        logger.info(f'Kneaddata read output retrieval: {directory}')
        filename_1 = f'{directory[1]}_kneaddata_paired_1.fastq'
        filename_2 = f'{directory[1]}_kneaddata_paired_2.fastq'
        file_1 = os.path.join(directory[0], filename_1)
        file_2 = os.path.join(directory[0], filename_2)
        copy_path_1 = os.path.join(kraken_source_directory, filename_1)
        copy_path_2 = os.path.join(kraken_source_directory, filename_2)
        if os.path.isfile(file_1) and os.path.isfile(file_2):
            logger.info(f'Kraken2 source files retrieved from: {directory}')
            shutil.copy(file_1, kraken_source_directory)
            shutil.copy(file_2, kraken_source_directory)
        kraken2_source_files[directory[1]] = (copy_path_1, copy_path_2)

    kraken_commands_list = []
    krakenresults_folders = []

    # Generate kraken2 commands
    # Generate commands for kraken2 database 1
    if kraken2_db1_path is not None:
        kraken2_commands_db1 = kraken2.generate_kraken2_commands(kraken2_source_files,
                                                                 kraken2_directory_db1,
                                                                 str(threads),
                                                                 kraken2_db1_path)
        kraken_commands_list.append((1, kraken2_commands_db1))
        krakenresults_folders.append((1, kraken2_directory_db1))

    # Generate commands for kraken2 database 2
    if kraken2_db2_path is not None:
        kraken2_commands_db2 = kraken2.generate_kraken2_commands(kraken2_source_files,
                                                                 kraken2_directory_db2,
                                                                 str(threads),
                                                                 kraken2_db2_path)
        kraken_commands_list.append((2, kraken2_commands_db2))
        krakenresults_folders.append((2, kraken2_directory_db2))

    # Run kraken2 commands
    for kraken_commands in kraken_commands_list:
        logger.info(f'Kraken2 db{kraken_commands[0]} run')
        for command in kraken_commands[1]:
            try:
                subprocess.run(command, check=True)
                logger.info(f'Kraken2 command db{kraken_commands[0]} successful')
            except subprocess.CalledProcessError as e:
                logger.error(f'Kraken2 command db{kraken_commands[0]} failed: {command}')

    if diamond_db_path is not None or targeted is not None:
        # Create krakenfile dictionary to be used to extract relevant krakenfiles
        krakenfiles = {}

        # Iterate through the kraken2 results folders
        for krakenfolder in krakenresults_folders:
            # Local dictionary to store krakenfiles
            krakenfile_dict = {}
            folders_results = os.listdir(krakenfolder[1])
            paths_results = [(os.path.join(krakenfolder[1], folder), folder)
                             for folder in folders_results
                             if os.path.isdir(os.path.join(krakenfolder[1], folder))]
            for directory in paths_results:
                logger.info(f'Krakenfile retrieval for db{krakenfolder[0]}: {directory}')
                logger.info(f'Kraken report retrieval for db{krakenfolder[0]}: {directory}')

                # Assign krakenfile and report names and filepaths
                krakenfile_name = f'{directory[1]}_krakenfile'
                report = f'{directory[1]}.txt'
                file_path_krakenfile = os.path.join(directory[0], krakenfile_name)
                file_path_report = os.path.join(directory[0], report)

                # Query filepaths for kraken2
                query_file_path_1 = os.path.join(kraken_source_directory,
                                                 f'{directory[1]}_kneaddata_paired_1.fastq')
                query_file_path_2 = os.path.join(kraken_source_directory,
                                                 f'{directory[1]}_kneaddata_paired_2.fastq')
                # Check if files exist
                if os.path.isfile(file_path_krakenfile) and os.path.isfile(file_path_report):
                    logger.info(f'kraken data added to krakenfiles dictionary: {directory[1]}')
                    krakenfile_dict[directory[1]] = (file_path_krakenfile,
                                                     query_file_path_1,
                                                     query_file_path_2,
                                                     file_path_report,
                                                     directory[0])
                else:
                    logger.error(f'kneaddata files necessary for kraken2 analysis not found: {directory[1]}')
            krakenfiles[krakenfolder[0]] = krakenfile_dict

        # This section of the analysis is for general metagenomic analysis of raw data using kraken2 from kneaddata
        # generated files. Unassigned files will also undergo further analysis using diamond blastx
        if diamond_db_path is not None and targeted is None:
            extract_reads('0', 'Unassigned_reads', krakenfiles, logger, extract_db_list=extract_db_list)
            # Generate contigs and interleaved file for unassigned reads
            # Create output file structure
            unassigned_dict = {}
            for db, reads_folder in enumerate(extract_db_list):
                unassigned_local = {}
                folders_extract = os.listdir(reads_folder)
                paths_extract = [(os.path.join(reads_folder, folder), folder)
                                 for folder in folders_extract
                                 if os.path.isdir(os.path.join(reads_folder, folder))]
                # Iterate through the unassigned reads folders
                for directory in paths_extract:
                    logger.info(f'Generating contigs for: {directory[1]}')
                    # Retrieve unassigned reads for further analysis

                    if os.path.isdir(os.path.join(directory[0], 'Unassigned_reads')):
                        output_directory = os.path.join(directory[0], 'Unassigned_reads')
                        if os.path.isfile(os.path.join(output_directory, 'Unassigned_reads_1.fastq')) and \
                                os.path.isfile(os.path.join(output_directory, 'Unassigned_reads_2.fastq')):
                            unassigned_read_1 = os.path.join(output_directory, 'Unassigned_reads_1.fastq')
                            unassigned_read_2 = os.path.join(output_directory, 'Unassigned_reads_2.fastq')

                            contigs_folder = os.path.join(output_directory, 'contigs')
                            if not os.path.exists(contigs_folder):
                                os.makedirs(contigs_folder)
                                logger.info(f'Contigs directory created: {contigs_folder}')

                            # Spades command to generate contigs
                            contigs_command = spades.meta_contigs(unassigned_read_1, unassigned_read_2, contigs_folder)
                            try:
                                subprocess.run(contigs_command, check=True)
                                logger.info(f'Generating contigs successful')
                            except subprocess.CalledProcessError as e:
                                logger.error(f'Generating contigs failed: {contigs_command}')

                            # create interleaved files
                            interleave_command = interleave.interleave_paired_reads(unassigned_read_1, unassigned_read_2,
                                                                                    os.path.join(output_directory,
                                                                                                 'interleaved.fasta.assembled.fastq'))

                            try:
                                subprocess.run(interleave_command, check=True)
                                logger.info(f'Interleaving reads successful')
                            except subprocess.CalledProcessError as e:
                                logger.error(f'Interleaving reads failed: {interleave_command}')

                            contigs_file = os.path.join(contigs_folder, 'contigs.fasta')
                            interleaved_file = os.path.join(output_directory, 'interleaved.fasta.assembled.fastq')
                            if os.path.isfile(contigs_file) and os.path.isfile(interleaved_file):
                                unassigned_local[directory[1]] = (contigs_file, interleaved_file)
                            elif os.path.isfile(interleaved_file):
                                unassigned_local[directory[1]] = (None, interleaved_file)
                            else:
                                unassigned_local[directory[1]] = (None, None)
                            unassigned_dict[db] = unassigned_local
                        else:
                            logger.error(f'Unassigned reads not found: {directory[1]}')
                    else:
                        logger.error(f'Unassigned reads folder not found: {directory[1]}')


            # Create diamond database
            create_db_command = diamond.diamond_build_db(diamond_db_path, output_folder, 'diamond_db')
            created_diamond_db = create_db_command[1]
            try:
                subprocess.run(create_db_command[0], check=True)
                logger.info(f'Create diamond database successful')
            except subprocess.CalledProcessError as e:
                logger.error(f'Create diamond database failed: {create_db_command[0]}')

            # Run diamond on unassigned reads
            for db, unassigned_dict in unassigned_dict.items():
                db_folder = os.path.join(diamond_output, f'kraken_db{db + 1}')
                if not os.path.exists(db_folder):
                    os.makedirs(db_folder)
                for sample, unassigned_files in unassigned_dict.items():
                    logger.info(f'Running diamond on unassigned reads for: {sample}')
                    output_directory = os.path.join(db_folder, sample)
                    if not os.path.exists(output_directory):
                        os.makedirs(output_directory)
                    if unassigned_files[0] is not None:
                        diamond_command_contigs = diamond.diamond_classify(created_diamond_db,
                                                                           unassigned_files[0],
                                                                           output_directory,
                                                                           f'{sample}_contigs')
                        try:
                            subprocess.run(diamond_command_contigs, check=True)
                            logger.info(f'Diamond command successful')
                        except subprocess.CalledProcessError as e:
                            logger.error(f'Diamond command failed: {diamond_command_contigs}')
                    if unassigned_files[1] is not None:
                        diamond_command_interleaved = diamond.diamond_classify(created_diamond_db,
                                                                               unassigned_files[1],
                                                                               output_directory,
                                                                               f'{sample}_interleaved')
                        try:
                            subprocess.run(diamond_command_interleaved, check=True)
                            logger.info(f'Diamond command successful')
                        except subprocess.CalledProcessError:
                            logger.error(f'Diamond command failed: {diamond_command_interleaved}')
            return output_directories
        elif diamond_db_path is None and targeted is not None:
            retrieved_files_r1 = []
            retrieved_files_r2 = []
            taxid_list = []
            parse_data = []
            read_1 = None
            read_2 = None
            krakenfile_retrieve = None
            if kraken2_db1_path is not None:
                save_path = extracted_reads_db1_directory
            elif kraken2_db2_path is not None:
                save_path = extracted_reads_db2_directory

            for db, kraken_reports in krakenfiles.items():
                key_item = f'{targeted}_1'
                if key_item in kraken_reports:
                    report_path = kraken_reports[key_item][3]
                    parse_data.extend(parse.generate_kraken2_report_csv(report_path, kraken_reports[key_item][4], logger))
                    read_1 = kraken_reports[key_item][1]
                    read_2 = kraken_reports[key_item][2]
                    krakenfile_retrieve = kraken_reports[key_item][0]

            for classification_item in parse_data:
                try:
                    if level in classification_item['Rank code']:
                        if int(classification_item['Number of reads assigned directly']) > 0:
                            if target_name in classification_item['Scientific name']:
                                taxid_list.append(classification_item['NCBI Taxonomy ID'])
                except KeyError:
                    logger.error(f'KeyError: {classification_item}')

            for taxid in taxid_list:
                logger.info(f'Extracting reads for: {taxid}')
                extract_reads_command = extract.extract_reads_command(krakenfile_retrieve,
                                                                      f'{target_name}_{taxid}',
                                                                      taxid,
                                                                      read_1,
                                                                      read_2,
                                                                      save_path)
                retrieved_files_r1.append(extract_reads_command[1][0])
                retrieved_files_r2.append(extract_reads_command[1][1])
                try:
                    subprocess.run(extract_reads_command[0], check=True)
                    logger.info(f'Extracting {taxid} reads successful')
                except subprocess.CalledProcessError as e:
                    logger.error(f'Extracting {taxid} reads failed: {extract_reads_command}')

            read_merge_1 = os.path.join(save_path, f'{target_name}_R1.fastq')
            read_merge_2 = os.path.join(save_path, f'{target_name}_R2.fastq')

            try:
                with open(read_merge_1, 'w') as outfile:
                    for fname in retrieved_files_r1:
                        with open(fname) as infile:
                            for line in infile:
                                outfile.write(line)
                with open(read_merge_2, 'w') as outfile:
                    for fname in retrieved_files_r2:
                        with open(fname) as infile:
                            for line in infile:
                                outfile.write(line)
                logger.info(f'Merged extracted kraken2 mapped reads for {target_name}')
                return read_merge_1, read_merge_2
            except IOError:
                logger.error(f'ERROR: Could not merge extracted reads for {target_name}')

    else:
        logger.error('No metagenomics use case specified. Ensure a kraken2 database is specified within the '
                     'reference_file when calling target_read_retrieval or a kraken2 database and diamond database is '
                     'specified when calling uct_meta')

