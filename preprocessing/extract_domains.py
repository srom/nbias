"""
Extract protein domains from Pfam or TIGR from raw file into
the assembly folder @ data/sequences/<assembly>/<assembly>_<pfam|tigr>.csv.gz
"""
import argparse
import os
import logging
import re
import gzip
import pathlib
import subprocess

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s (%(levelname)s) %(message)s")

    parser = argparse.ArgumentParser()
    parser.add_argument('domain_type', type=str, choices=['pfam', 'tigr'])
    parser.add_argument('--batch_size', type=int, default=int(1e4))
    args = parser.parse_args()

    domain_type = args.domain_type
    batch_size = args.batch_size

    logger.info(f'Extracting {domain_type} protein domains (batch size: {batch_size:,})')

    assemblies_metadata_path = os.path.join(os.getcwd(), 'data/assemblies.csv')
    domain_folder = os.path.join(os.getcwd(), f'data/{domain_type}')
    output_base_path = os.path.join(os.getcwd(), 'data/sequences')

    assemblies_set = set(pd.read_csv(assemblies_metadata_path)['assembly_accession'].unique())

    process_batch = process_batch_fn(domain_type, output_base_path)

    for batch_df in process_files(domain_folder, assemblies_set, batch_size):
        process_batch(batch_df)

    logger.info('Post processing: compress domain files')
    for i, assembly in enumerate(sorted(assemblies_set)):
        if i == 0 or (i+1) % 1000 == 0:
            logger.info(f'Compressing file {i+1:,} / {len(assemblies_set):,}')

        path = os.path.join(output_base_path, f'{assembly}/{assembly}_{domain_type}.csv')
        if os.path.isfile(path):
            subprocess.run(['gzip', '-f', path], check=True)

    logger.info('DONE')


def process_files(folder, assemblies_set, batch_size, skiplines=4, n_cols=19):
    """
    Read alignment output for pfam or tigr.
    The first 4 lines need to be skipped.
    """
    files = []
    for p in pathlib.Path(folder).glob('*.txt.gz'):
        if p.is_file():
            files.append(p)

    logger.info(f'Processing {len(files)} files')

    p = r'\s+'.join([r'([^\s]+)' for _ in range(n_cols)])
    pattern = f'^{p}$'
    
    batch = []
    n_records = 0
    for i, path in enumerate(sorted(files, key=lambda p: p.name)):
        logger.info(f'Processing file {i+1} / {len(files)}: {path.name}')

        line_nb = 0
        with gzip.open(str(path), 'rt') as f:
            for line in f:
                line_nb += 1
                if line_nb < skiplines:
                    continue

                m = re.match(pattern, line)
                if m is None:
                    continue

                n_records += 1

                if n_records % int(1e5) == 0:
                    logger.info(f'{n_records:,} records processed')

                row = [m[i+1] for i in range(n_cols)]
                
                first_el = row[0]
                
                a, accession = tuple(first_el.split('$'))

                if accession not in assemblies_set:
                    continue

                _, protein_id_raw = tuple(a.split('@'))
                label = row[-1] if row[-1] != '-' else None
                
                protein_id = protein_id_raw.replace('-', '_')
                query = row[2]
                query_id = row[3]
                
                data_row = [
                    accession,
                    protein_id,
                    query,
                    query_id,
                    label,
                ]
                batch.append(data_row)
                
                if len(batch) >= batch_size:
                    yield prepare_batch(batch)
                    batch = []

        if len(batch) > 0:
            yield prepare_batch(batch)
            batch = []

    if len(batch) > 0:
        yield prepare_batch(batch)

    logger.info(f'Total number of records: {n_records:,}')

    return


def prepare_batch(batch):
    return pd.DataFrame(batch, columns=[
        'assembly_accession',
        'protein_id',
        'query',
        'query_id',
        'label',
    ]).set_index('assembly_accession')


def process_batch_fn(domain_type, output_base_path):
    seen_assemblies = set()

    def fn(batch_df):
        assemblies = sorted(set(batch_df.index.tolist()))
        for assembly in assemblies:
            assembly_df = batch_df.loc[[assembly]]

            if len(assembly_df) == 0:
                continue

            path = os.path.join(output_base_path, f'{assembly}/{assembly}_{domain_type}.csv')

            mode = 'a'
            header = False
            if assembly not in seen_assemblies:
                mode = 'w'
                header = True
                seen_assemblies.add(assembly)

            assembly_df.to_csv(path, index=True, mode=mode, header=header)

    return fn


if __name__ == '__main__':
    main()
