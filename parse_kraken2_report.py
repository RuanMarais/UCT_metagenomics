"""
******************************************************************************************************************
parse_kraken2_report.py
Kraken2 report interpretation
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""

import csv
import os


def generate_kraken2_report_csv(report, output_directory, logger):
    """
    This function parses Kraken2 reports and generates a CSV file as well as a dictionary data object for further
    analysis

    :param report: The kraken2 report filepath
    :param output_directory: The output directory filepath for the csv file
    :return: The dictionary with keys: Percentage of reads covered, Number of reads covered,
    Number of reads assigned directly, Rank code, NCBI Taxonomy ID, Scientific name
    """
    # Open the Kraken 2 report file
    try:
        with open(report, 'r') as report_file:
            # Open the output CSV file
            try:
                report_csv = os.path.join(output_directory, 'kraken2_report_csv.csv')
                with open(report_csv, 'w', newline='') as csv_file:
                    # Create a CSV writer
                    writer = csv.writer(csv_file)
                    # Write the header
                    headers = ['Percentage of reads covered', 'Number of reads covered', 'Number of reads assigned directly',
                               'Rank code', 'NCBI Taxonomy ID', 'Scientific name']
                    writer.writerow(headers)
                    # Read the Kraken 2 report line by line
                    data = []
                    for line in report_file:
                        # Split the line into fields
                        fields = line.strip().split('\t')
                        # Remove leading spaces from the scientific name
                        fields[-1] = fields[-1].lstrip()
                        # Write the fields to the CSV file
                        writer.writerow(fields)
                        # Append data to list as dictionary
                        row_dict = dict(zip(headers, fields))
                        data.append(row_dict)
                    return data
            except FileNotFoundError:
                logger.error(f'Could not find kraken2 report output directory at {output_directory}')
                return None
    except FileNotFoundError:
        logger.error(f'Kraken2 report not found at {report}')
        return None


