import unicodecsv
import xlsxwriter
from collections import defaultdict
import numpy as np

# The list of runs to be included in the analysis. The run name should match the data in the metadata csv

# run_list = ['1', '2', '4', '3', '5', '6', '7a', '7b']
run_list = ['4', '5', '6', '7a', '7b']
# run_list = ['3']

# Settings to be used for the analysis

z_score_threshold = 2
minimum_reads = 2
minimum_reads_genus = 2
flag_species = ['Human immunodeficiency virus 1']
flag_genus = ['Mycobacterium']
exclude_species = ['Homo sapiens']
exclude_genus = ['Homo']

# The input files are flagged here as lists. Each item is a csv file with species/genus in the
# first column and each sample is in a new column. Read counts are given. If no reads are present
# the cell contains 'NA'
# The positive control should be excluded from this analysis

genus_reads = ['Run4_genus.csv', 'Run5_genus.csv', 'Run6_genus.csv', 'Run7a_genus.csv',
               'Run7b_genus_minus12.csv']
species_reads = ['Run4_species.csv', 'Run5_species.csv', 'Run6_species.csv', 'Run7a_species.csv',
                 'Run7b_species_minus12.csv']

# genus_reads = ['Run1_genus.csv', 'Run2_genus.csv', 'Run3_genus.csv', 'Run4_genus.csv',
#                'Run5_genus.csv', 'Run6_genus.csv', 'Run7a_genus.csv', 'Run7b_genus.csv']
# species_reads = ['Run1_species.csv', 'Run2_species.csv', 'Run3_species.csv', 'Run4_species.csv',
#                  'Run5_species.csv', 'Run6_species.csv', 'Run7a_species.csv', 'Run7b_species.csv']
# genus_reads = ['Run3_genus.csv']
# species_reads = ['Run3_species.csv']

# The metadata needed for the run
metadata_raw = ['metadata_new.csv']

# This function takes csv files and outputs a list object. The first row of the csv is used as the headings.
# The files ('.csv') are added to the file_list variable above.


def importer(file_list_csvs):
    output = []
    for file in file_list_csvs:
        with open(file, 'rb') as f:
            reader = unicodecsv.DictReader(f)
            output.extend(list(reader))
    return output


# The variable that will hold the episodes is defined here
genus_data = importer(genus_reads)
species_data = importer(species_reads)
metadata = importer(metadata_raw)

sample_dict = {}

for val in metadata:
    sample = val['Sample']
    sample_dict[sample] = {}
    genus_dict = {}
    species_dict = {}
    result_dict = {}
    sample_dict[sample]['metadata'] = val
    sample_dict[sample]['genus_dict'] = genus_dict
    sample_dict[sample]['species_dict'] = species_dict
    sample_dict[sample]['results'] = result_dict

for genuses in genus_data:
    genus = genuses['name']
    for key in sample_dict.keys():
        if key in genuses:
            value = genuses[key]
            if value == 'NA':
                float_val = 0
            else:
                try:
                    float_val = float(value)
                except:
                    print(key)
                    print('val:' + value)

            nonhuman_reads = float(sample_dict[key]['metadata']['nonhuman reads'])
            ratio = float_val / nonhuman_reads
            sample_dict[key]['genus_dict'][genus] = [float_val, ratio]

for species_vals in species_data:

    species = species_vals['name']
    for key in sample_dict.keys():

        if key in species_vals:
            value = species_vals[key]
            if value == 'NA':
                float_val = 0
            else:
                try:
                    float_val = float(value)
                except:
                    print(key)
                    print('val:' + value)

            nonhuman_reads = float(sample_dict[key]['metadata']['nonhuman reads'])
            ratio = float_val / nonhuman_reads
            sample_dict[key]['species_dict'][species] = [float_val, ratio]

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
    run = sample_dict[key]['metadata']['Run']

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
    run = sample_dict[key]['metadata']['Run']
    for genus, vals in sample_dict[key]['genus_dict'].items():
        neg_val = sample_dict['NC_r' + run]['genus_dict'][genus][1]
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
            z_score_neg = None
            z_score_higher_than_nc_for_run = None
        sample_dict[key]['genus_dict'][genus].extend([z_score_all, z_score_run, z_score_higher_than_nc_for_run])
    for species, vals in sample_dict[key]['species_dict'].items():
        neg_val = sample_dict['NC_r' + run]['species_dict'][species][1]
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
            z_score_neg = None
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
                if species_results[2] > z_score_threshold and species_results[3] > z_score_threshold and species_results[0] >= minimum_reads:
                    if species_results[4] is not None:
                        if species_results[4]:
                            likely_pathogen[species_key] = species_results
                elif species_results[2] > z_score_threshold or species_results[3] > z_score_threshold:
                    if species_results[4] is not None:
                        if species_results[4]:
                            possible_pathogen[species_key] = species_results
                elif species_key in flag_species:
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


def generate_export_list(nested_dict, keys, headings):
    output = []
    for key in keys:
        output_add = [headings]
        for org, val in nested_dict[key].items():
            add_val = [org]
            add_val.extend(val)
            output_add.append(add_val)
        output.append(output_add)
    return output


def write_to_excel(nested_dict, file_name, keys, headings):
    workbook = xlsxwriter.Workbook(file_name)
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


for participant in sample_dict.keys():
    write_to_excel(sample_dict[participant]['results'], participant + '.xlsx', keys_print, column_headings)