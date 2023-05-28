"""
******************************************************************************************************************
diamond_assign.py
Functions to generate commands for the DIAMOND classification tool
UCT_metagenomics pipeline
Author: Gert Marais, University of Cape Town, 2023
******************************************************************************************************************
"""
import os


def diamond_build_db(db_fasta,
                     output_dir,
                     output_prefix):
    """
    This function uses the DIAMOND 'makedb' command to create a DIAMOND database from a fasta file

    :param db_fasta: The fasta file that will be converted into the DIAMOND database
    :param output_dir: The directory filepath to save the created database
    :param output_prefix: The created database name prefix
    :return: The DIAMOND 'makedb' command to run using subprocess
    """
    db_path = os.path.join(output_dir, output_prefix)
    return ['diamond', 'makedb', '--in', db_fasta, '-d', db_path], db_path


def diamond_classify(database,
                     query,
                     output_dir,
                     output_prefix):
    """
    This function uses the DIAMOND 'blastx' command to classify a query file using a DIAMOND database

    :param database: The DIAMOND database filepath
    :param query: The query filepath to be classified
    :param output_dir: The output directory filepath for the classification output
    :param output_prefix: The naming prefix for the classification output
    :return: The DIAMOND 'blastx' command to run using subprocess
    """
    return ['diamond', 'blastx', '-d', database, '-q', query, '-o',
            os.path.join(output_dir, output_prefix + '_diamond.tsv')]
