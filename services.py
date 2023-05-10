"""
******************************************************************************************************************
General helper functions for the pipeline
******************************************************************************************************************
"""

from collections import defaultdict
import os
import logging

# Logging
logger = logging.getLogger(__name__)
handler = logging.FileHandler("Run.log")
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)


def filename_list_generate(file_identifier, folder_path):
    filenames_output = []
    for file in os.listdir(folder_path):
        if os.path.isfile(os.path.join(folder_path, file)):
            if file_identifier in file:
                filenames_output.append(file)
    return filenames_output


def paired_read_list_generate(id_length, file_id_r1, file_id_r2, filename_list, file_folder):
    # Separate filenames by id
    sample_sorted_dict = defaultdict(list)
    for file in filename_list:
        unique_filename = file[:id_length]
        sample_sorted_dict[unique_filename].append(file)

    # Separate pairs
    paired_sorted_dict = {}
    for key, files in sample_sorted_dict.items():
        logger.debug(f'Read key: {key}')
        read_1 = None
        read_2 = None
        for file in files:
            if file_id_r1 in file:
                read_1 = file
                logger.debug(f'Read 1 found: {key}')
            elif file_id_r2 in file:
                read_2 = file
                logger.debug(f'Read 2 found: {key}')
        if read_1 is not None and read_2 is not None:
            read_out_1 = os.path.join(file_folder, read_1)
            read_out_2 = os.path.join(file_folder, read_2)
            output = (read_out_1, read_out_2)
            logger.debug(f'Reads assigned to key: {key}')
            paired_sorted_dict[key] = output

    return paired_sorted_dict
