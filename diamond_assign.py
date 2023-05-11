"""
******************************************************************************************************************
This this script uses DIAMOND to classify reads
******************************************************************************************************************
"""

import os


def diamond_build_db(db_fasta, output_dir, output_prefix):
    db_path = os.path.join(output_dir, output_prefix)
    return ['diamond', 'makedb', '--in', db_fasta, '-d', db_path], db_path


def diamond_classify(database, query, output_dir, output_prefix):
    return ['diamond', 'blastx', '-d', database, '-q', query, '-o',
            os.path.join(output_dir, output_prefix + '_diamond.tsv')]
