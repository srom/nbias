"""
Split single protein fasta file into one fasta file per assembly.
"""
import argparse
import os
import logging
import gzip
import re

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


logger = logging.getLogger(__name__)


ASSEMBLY_RE = r'^.*@([^\$]+)\$(GCA_[0-9\.]+).*$'


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s (%(levelname)s) %(message)s")

    parser = argparse.ArgumentParser()
    parser.add_argument('input_fasta', type=str)
    parser.add_argument('--output_folder', type=str, default=None)
    args = parser.parse_args()

    input_fasta = args.input_fasta
    output_folder = args.output_folder

    if output_folder is None:
        output_folder = os.path.join(os.getcwd(), 'data/sequences')

    logger.info(f'Splitting fasta file {input_fasta}')

    n_files = 0
    current_assembly = None
    current_records = []
    for i, record in enumerate(SeqIO.parse(input_fasta, 'fasta')):
        if i == 0 or (i+1) % int(10e5) == 0:
            logger.info(f'Processing record {i+1:,} | {n_files:,} files saved')

        try:
            assembly, record_id = parse_assembly(record)
        except ValueError as e:
            logger.warning(e)
            continue

        if assembly != current_assembly:
            n_files += save_fasta(output_folder, current_assembly, current_records)
            current_assembly = assembly
            current_records = []

        current_records.append(
            make_new_record(current_assembly, record_id, record)
        )

    n_files += save_fasta(output_folder, current_assembly, current_records)

    logger.info(f'Final number of files: {n_files:,}')
    logger.info('DONE')


def make_new_record(assembly, record_id, record):
    return SeqRecord(
        record.seq,
        id=record_id,
        description=assembly
    )


def save_fasta(output_folder, assembly, records):
    if len(records) == 0:
        return 0

    output_path = os.path.join(
        output_folder, 
        f'{assembly}/{assembly}_protein.faa.gz',
    )

    with gzip.open(output_path, 'wt') as f_gz:
        SeqIO.write(records, f_gz, 'fasta')

    return 1


def parse_assembly(record):
    m = re.match(ASSEMBLY_RE, record.id)
    if m is None:
        raise ValueError(f'No assembly or protein id found for record {record.id}')

    record_id, assembly = m[1], m[2]
    return (
        assembly, 
        record_id.replace('-', '_'),
    )


if __name__ == '__main__':
    main()
