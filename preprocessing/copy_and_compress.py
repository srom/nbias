"""
Copy and compress assembly files matching a pathlib glob pattern
into their named folder on data/sequences.
"""
import argparse
import os
import logging
import gzip
import shutil
import re
from pathlib import Path

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


ADMISSIBLE_EXTENSIONS_RE = r'^.+\.(fna|faa|gff)$'


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s (%(levelname)s) %(message)s")

    parser = argparse.ArgumentParser()
    parser.add_argument('file_type', type=str, choices=['genome', 'gff', 'CDS', 'proteome'])
    parser.add_argument('pattern', type=str)
    parser.add_argument('--base_folder', type=str, default=None)
    args = parser.parse_args()

    file_type = args.file_type
    pattern = args.pattern
    base_folder = args.base_folder

    if base_folder is None:
        base_folder = os.path.join(os.getcwd(), 'data')

    logger.info('Copy & Compress')
    logger.info(f'File type: {file_type}')
    logger.info(f'Base folder: {base_folder}')
    logger.info(f'Pattern: {pattern}')

    paths = list(Path(base_folder).glob(pattern))

    logger.info(f'Copying and compressing {len(paths):,} files')

    for i, path in enumerate(paths):
        if i == 0 or (i+1) % 100 == 0:
            logger.info(f'Processing file {i+1:,} / {len(paths):,}')

        if re.match(ADMISSIBLE_EXTENSIONS_RE, path.name) is None:
            logger.warning(f'Invalid extension for file {path.name}')
            continue

        try:
            assembly = parse_assembly(path)
        except ValueError as e:
            logger.warning(e)
            continue

        output_folder = os.path.join(os.getcwd(), f'data/sequences/{assembly}')
        Path(output_folder).mkdir(parents=False, exist_ok=True)

        if file_type == 'genome':
            output_path = os.path.join(output_folder, f'{assembly}_genomic.fna.gz')
        elif file_type == 'gff':
            output_path = os.path.join(output_folder, f'{assembly}_genomic.gff.gz')
        elif file_type == 'CDS':
            output_path = os.path.join(output_folder, f'{assembly}_cds_from_genomic.fna.gz')
        elif file_type == 'proteome':
            output_path = os.path.join(output_folder, f'{assembly}_protein.faa.gz')
        else:
            raise ValueError(f'Unknown file type {file_type}')

        with path.open('rb') as f_in:
            with gzip.open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    logger.info('DONE')


def parse_assembly(path):
    m = re.match(r'^(GCA_[^_]+)_.+$', path.name)
    if m is not None:
        return m[1]
    else:
        raise ValueError(f'No assembly found in filename: {path.name}')


if __name__ == '__main__':
    main()
