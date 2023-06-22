"""
******************************************************************************************************************
zscore_analyse.py
This script uses the primary kraken2 classification output to generate a z-score analysis of the data
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""

import unicodecsv
import xlsxwriter
from collections import defaultdict
import numpy as np
from data_and_variables import z_score_threshold, minimum_reads, minimum_reads_genus, flag_species, flag_genus, \
    exclude_species, exclude_genus
import os


def process_read_data(tax_dict,
                      sample_dict,
                      sample_dict_field,
                      logger):
    """
    This function processes the read data from the kraken2 report and adds it to the sample dictionary

    :param tax_dict: The dictionary of species or genus data derived from a kraken2 report
    :param sample_dict: The sample data dictionary initialised from the metadata file
    :param logger: The logging object to generate the logfile
    """
    for tax, sample_data in tax_dict.items():
        for sample, sample_info in sample_dict.items():
            try:
                nonhuman_reads = float(sample_dict[sample]['metadata']['nonhuman_reads'])
                located = False
                for data_item in sample_data:
                    if data_item[0] == f'{sample}_R1':
                        read_data = float(data_item[1])
                        located = True
                if located:
                    ratio = read_data / nonhuman_reads
                    sample_dict[sample][sample_dict_field][tax] = [read_data, ratio]
                else:
                    sample_dict[sample][sample_dict_field][tax] = [0, 0]
            except KeyError:
                logger.error(f'Error accessing metadata for {sample}')


def importer(file, logger):
    """
    Helper function to import a csv file as a list of dictionaries

    :param logger: A logging object to generate the logfile
    :param file: csv file to be imported
    :return: list of dictionaries
    """
    output = []
    try:
        with open(file, 'rb') as f:
            reader = unicodecsv.DictReader(f)
            output.extend(list(reader))
            logger.info(f'Imported {file}')
        return output
    except IOError:
        logger.error(f'Error: File {file} not found')
        return None


def generate_export_list(nested_dict,
                         keys,
                         headings):
    """
    Helper function to generate the list of lists to be exported to excel
    Each top-level key in the nested dictionary corresponds to a separate sheet in the output Excel file.
    The headings parameter provides the column headers for each sheet.

    :param nested_dict: dict, A dictionary of dictionaries. Each top-level key in the nested_dict corresponds to a
    separate sheet in the output Excel file. The second level keys correspond to row identifiers,
    and their corresponding values are lists that are the remaining cell values in that row.
    :param keys: list, The keys for the top level of the nested dictionary that we want to include in the output.
    :param headings: list, A list of column header names for the Excel sheets. This should be a list of strings.

    :return: list, A list of lists. Each inner list represents a sheet in the Excel file, with the first sublist
    representing the column headers and the remaining sublists representing the rows in that sheet.

    """
    output = []
    for key in keys:
        output_add = [headings]
        for org, val in nested_dict[key].items():
            add_val = [org]
            add_val.extend(val)
            output_add.append(add_val)
        output.append(output_add)
    return output


def write_to_excel(nested_dict,
                   file_name,
                   keys,
                   headings,
                   output_folder,
                   logger):
    """
    Writes a nested dictionary to an Excel workbook. Each key in the outer dictionary
    corresponds to a sheet in the workbook. The inner dictionaries' key-value pairs
    correspond to rows within the respective sheet.

    :param nested_dict: dict, The nested dictionary to write to Excel. The outer dictionary's keys will be used as
    worksheet names, and the inner dictionaries will be written to their respective worksheets.
    :param file_name: str, The name of the Excel file to be created. Must include the '.xlsx' extension.
    :param keys: list, List of keys from the nested dictionary that define the order of worksheets within the created
    Excel workbook.
    :param headings: list, The headings to be used for each column in the worksheets. The order of the headings should
    correspond to the order of the values in the inner dictionaries of the nested dictionary.
    :param output_folder: str, The directory where the Excel file will be created. Must be a valid system path.
    :param logger: The logging object to generate the logfile
    :returns: None The function writes the contents of the nested_dict to an Excel file and does not return a value.
    """
    try:
        output_file = os.path.join(output_folder, file_name)
        workbook = xlsxwriter.Workbook(output_file)
    except FileNotFoundError:
        logger.error(f'Error: {output_folder} does not exist or is not a directory')
    worksheets = keys
    list_of_lists = generate_export_list(nested_dict, keys, headings)
    count = 0
    for list_stats in list_of_lists:
        worksheet = workbook.add_worksheet(worksheets[count])
        count += 1
        row = 0
        col = 0
        worksheet.set_column(0, 5, 25)
        for list_item in list_stats:
            local_col = col
            for item in list_item:
                worksheet.write(row, local_col, item)
                local_col += 1
            row += 1
    workbook.close()


def z_score_analysis(sample_dict,
                     run_list,
                     output_folder,
                     logger):
    """

    Performs a z-score analysis on a set of samples. The z-score is a measure that describes
    a value's relationship to the mean of a group of values. This function computes the z-score
    for both genus and species level data of each sample, and writes the results to an Excel
    spreadsheet.

    :param sample_dict: dict, A dictionary containing information about each sample. The keys are sample names,
        and the values are another dictionary that contains 'genus_dict' and 'species_dict'
        for each sample. Each of these inner dictionaries has genus/species names as keys which are generated by the
        process_read_data function
    :param run_list: list, A list containing the names of the runs to be included in the analysis.
    :param output_folder: str, The directory where the resulting Excel files will be saved. Must be a valid system path.
    :param logger: The logging object to generate the logfile
    :return: None, The function does not return any value. It writes the z-score analysis results to Excel files.
    """
    # Data containers for the z-score analysis
    all_dict_species = defaultdict(list)
    all_dict_genus = defaultdict(list)
    run_dict_genus = {}
    run_dict_species = {}
    z_score_dict_genus = {}
    z_score_dict_species = {}

    for val in run_list:
        vals_add_genus = defaultdict(list)
        vals_add_species = defaultdict(list)
        run_dict_species[val] = vals_add_species
        run_dict_genus[val] = vals_add_genus

    for key in sample_dict.keys():
        run = sample_dict[key]['metadata']['run']

        for genus_item_key in sample_dict[key]['genus_dict'].keys():
            read_corrected = sample_dict[key]['genus_dict'][genus_item_key][1]
            all_dict_genus[genus_item_key].append(read_corrected)
            run_dict_genus[run][genus_item_key].append(read_corrected)

        for species_item_key in sample_dict[key]['species_dict'].keys():
            read_corrected = sample_dict[key]['species_dict'][species_item_key][1]
            all_dict_species[species_item_key].append(read_corrected)
            run_dict_species[run][species_item_key].append(read_corrected)

    z_score_dict_genus['all'] = {}
    z_score_dict_species['all'] = {}

    for key, vals in all_dict_genus.items():
        all_dict = {'mean': np.mean(vals), 'std': np.std(vals)}
        z_score_dict_genus['all'][key] = all_dict

    for key, vals in all_dict_species.items():
        all_dict = {'mean': np.mean(vals), 'std': np.std(vals)}
        z_score_dict_species['all'][key] = all_dict

    for key in run_dict_genus.keys():
        z_score_dict_genus[key] = {}
        for genus, items in run_dict_genus[key].items():
            gen_dict = {'mean': np.mean(items), 'std': np.std(items)}
            z_score_dict_genus[key][genus] = gen_dict

    for key in run_dict_species.keys():
        z_score_dict_species[key] = {}
        for species, items in run_dict_species[key].items():
            gen_dict = {'mean': np.mean(items), 'std': np.std(items)}
            z_score_dict_species[key][species] = gen_dict

    for key in sample_dict.keys():
        run = sample_dict[key]['metadata']['run']
        for genus, vals in sample_dict[key]['genus_dict'].items():
            neg_val = sample_dict['NCrun' + run]['genus_dict'][genus][1]
            mean_all = z_score_dict_genus['all'][genus]['mean']
            std_all = z_score_dict_genus['all'][genus]['std']
            mean_run = z_score_dict_genus[run][genus]['mean']
            std_run = z_score_dict_genus[run][genus]['std']
            if std_all != 0:
                z_score_all = (vals[1] - mean_all) / std_all
            else:
                z_score_all = None
            if std_run != 0:
                z_score_run = (vals[1] - mean_run) / std_run
                z_score_neg = (neg_val - mean_run) / std_run
                z_score_higher_than_nc_for_run = z_score_run > z_score_neg
            else:
                z_score_run = None
                z_score_higher_than_nc_for_run = None
            sample_dict[key]['genus_dict'][genus].extend([z_score_all, z_score_run, z_score_higher_than_nc_for_run])
        for species, vals in sample_dict[key]['species_dict'].items():
            neg_val = sample_dict['NCrun' + run]['species_dict'][species][1]
            mean_all = z_score_dict_species['all'][species]['mean']
            std_all = z_score_dict_species['all'][species]['std']
            mean_run = z_score_dict_species[run][species]['mean']
            std_run = z_score_dict_species[run][species]['std']
            if std_all != 0:
                z_score_all = (vals[1] - mean_all) / std_all
            else:
                z_score_all = None
            if std_run != 0:
                z_score_run = (vals[1] - mean_run) / std_run
                z_score_neg = (neg_val - mean_run) / std_run
                z_score_higher_than_nc_for_run = z_score_run > z_score_neg
            else:
                z_score_run = None
                z_score_higher_than_nc_for_run = None
            sample_dict[key]['species_dict'][species].extend([z_score_all, z_score_run, z_score_higher_than_nc_for_run])

    for key in sample_dict.keys():
        likely_pathogen = {}
        possible_pathogen = {}
        genus_of_note = {}
        for species_key in sample_dict[key]['species_dict'].keys():
            species_results = sample_dict[key]['species_dict'][species_key]
            if species_key not in exclude_species:
                if species_results[3] is not None:
                    if species_key in flag_species:
                        if species_results[4] is not None:
                            if species_results[4]:
                                likely_pathogen[species_key] = species_results
                    elif species_results[2] > z_score_threshold and species_results[3] > z_score_threshold and \
                            species_results[0] >= minimum_reads:
                        if species_results[4] is not None:
                            if species_results[4]:
                                likely_pathogen[species_key] = species_results
                    elif species_results[2] > z_score_threshold or species_results[3] > z_score_threshold:
                        if species_results[4] is not None:
                            if species_results[4]:
                                possible_pathogen[species_key] = species_results

        for genus_key in sample_dict[key]['genus_dict'].keys():
            genus_results = sample_dict[key]['genus_dict'][genus_key]
            if genus_key not in exclude_genus:
                if genus_results[3] is not None:
                    if genus_key in flag_genus:
                        if genus_results[0] > 0:
                            genus_of_note[genus_key] = genus_results
                    elif genus_results[2] > z_score_threshold or genus_results[3] > z_score_threshold:
                        if genus_results[0] >= minimum_reads_genus:
                            if genus_results[4] is not None:
                                if genus_results[4]:
                                    genus_of_note[genus_key] = genus_results

        sample_dict[key]['results']['likely_pathogen'] = likely_pathogen
        sample_dict[key]['results']['possible_pathogen'] = possible_pathogen
        sample_dict[key]['results']['genus_of_note'] = genus_of_note

    keys_print = ['likely_pathogen', 'possible_pathogen', 'genus_of_note']
    column_headings = ['Species or Genus', 'Reads', 'Reads per nonhuman read', 'z score all runs', 'z score within run',
                       'z score greater than negative control z score']

    # Generate excel files for each sample with the pathogen z-score data
    for participant in sample_dict.keys():
        write_to_excel(sample_dict[participant]['results'],
                       f'{participant}_results.xlsx',
                       keys_print,
                       column_headings,
                       output_folder,
                       logger)

    # Generate a general excel file with the read counts for each sample
    gen_data_output = []
    general_results_output = os.path.join(output_folder, 'general_results.xlsx')
    column_headers_general = ['Sample', 'Total reads', 'Non-human reads']
    gen_data_output.append(column_headers_general)
    workbook = xlsxwriter.Workbook(general_results_output)
    worksheet = workbook.add_worksheet()

    for sample, data in sample_dict.items():
        total_reads = data['metadata']['total_reads']
        nonhuman_reads = data['metadata']['nonhuman_reads']
        output_item = [sample, total_reads, nonhuman_reads]
        gen_data_output.append(output_item)

    for i, row_data in enumerate(gen_data_output):
        for j, cell_data in enumerate(row_data):
            worksheet.write(i, j, cell_data)

    workbook.close()
