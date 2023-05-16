"""
******************************************************************************************************************
This script extracts reads from a krakenfile matching a specific taxid
******************************************************************************************************************
"""
import os


def extract_reads_command(krakenfile_path, org, taxid, input_read_1_path, input_read_2_path, output_dir):
    org_output = os.path.join(output_dir, org)
    if not os.path.exists(org_output):
        os.makedirs(org_output)
    output_paths = (os.path.join(org_output, org + '_1'), os.path.join(org_output, org + '_2'))
    return ['extract_kraken_reads.py', '-k', krakenfile_path, '-t', taxid, '--fastq-output', '-s', input_read_1_path,
            '-s2', input_read_2_path, '-o', os.path.join(org_output, org + '_1.fastq'), '-o2',
            os.path.join(org_output, org + '_2.fastq')], output_paths, org_output


def extract_targeted_reads(krakenfile_dict_1, krakenfile_dict_2, target_dict, raw_reads_dict, output_dir):
    output_data = []
    for id_item, reference_data in target_dict.items():
        id_output = os.path.join(output_dir, id_item)
        if not os.path.exists(id_output):
            os.makedirs(id_output)
        taxid = reference_data[1]
        kraken_dict = reference_data[4]
        organism = reference_data[0]
        raw_read_1 = raw_reads_dict[id_item][0]
        raw_read_2 = raw_reads_dict[id_item][1]
        reference = reference_data[2]

        krakenfile = None
        if kraken_dict == 1:
            krakenfile = krakenfile_dict_1[id_item]
        elif kraken_dict == 2:
            krakenfile = krakenfile_dict_2[id_item]

        if krakenfile is not None:
            extract_data = extract_reads_command(krakenfile, organism, taxid, raw_read_1, raw_read_2, output_dir)[0]
            command = extract_data[0]
            output_paths = extract_data[1]
            out_dir = extract_data[2]
            output = [id_item, command, output_paths, organism, reference, out_dir]
            output_data.append(output)
    return output_data

