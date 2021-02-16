"""
Compute the average Jensen-Shannon distance from genome-wide to tRNA synthetase genes
across all complete genomes in our dataset.
"""
import os
from os.path import join
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s (%(levelname)s) %(message)s")

    logger.info('Fetching assemblies')

    cwd = os.getcwd()
    assemblies_path = join(cwd, 'data/assemblies.csv')
    assembly_df = pd.read_csv(assemblies_path, index_col='assembly_accession')

    logger.info(f'Number of assemblies: {len(assembly_df):,}')

    complete_genome_df = assembly_df[
        (assembly_df['phylum'].notnull()) & 
        (assembly_df['assembly_level'] == 'Complete Genome')
    ]

    perc_total = 100 * len(complete_genome_df) / len(assembly_df) 
    logger.info(f'Number of complete genomes: {len(complete_genome_df):,} ({perc_total:.0f}%)')

    logger.info('Identify tRNA synthetases based on TIGR labels')

    tigr_master = pd.read_csv(join(cwd, 'data/tigr_master.csv'), index_col='id')
    tigr_trna_synt = tigr_master[
        tigr_master['description'].str.contains('--tRNA ligase')
    ].index.tolist()

    distances = []
    stds = []
    data = []
    for i, assembly in enumerate(complete_genome_df.index):
        if i == 0 or (i+1) % 100 == 0:
            logger.info(f'{i+1} / {len(complete_genome_df)}')

        row = complete_genome_df.loc[assembly]
        species = row['organism_name']
        phylum = row['phylum']

        tigr_map = pd.read_csv(
            join(cwd, f'data/sequences/{assembly}/{assembly}_tigr.csv.gz'),
            index_col='query'
        )
        protein_ids = tigr_map.loc[
            sorted(set(tigr_map.index.tolist()) & set(tigr_trna_synt))
        ]['protein_id'].unique()

        protein_dist = pd.read_csv(
            join(cwd, f'data/sequences/{assembly}/{assembly}_tri_nucleotide_distance_to_mean.csv'),
            index_col='protein_id',
        )

        protein_ids_present = set(protein_ids) & set(protein_dist.index.tolist())

        assembly_distances = [
            v for v in protein_dist.loc[protein_ids_present]['distance'].values
            if not pd.isnull(v)
        ]
        if len(assembly_distances) < 10:
            continue

        std = np.std(assembly_distances)

        distances += assembly_distances

        if pd.isnull(std):
            continue

        stds.append(std)

        n_prot = len(assembly_distances)
        m_dist = np.mean(assembly_distances)

        data.append([assembly, species, phylum, n_prot, m_dist, std])

    res_df = pd.DataFrame(data, columns=[
        'assembly_accession',
        'organism_name',
        'phylum',
        'n_trna_synt_genes',
        'distance_mean',
        'distance_std',
    ])

    logger.info('Results:')
    logger.info(f'Average number of proteins: {len(distances) / len(complete_genome_df):.2f}')
    logger.info(f'Average distance: {np.mean(distances):.3f}')
    logger.info(f'Standard deviation: {np.std(distances):.3f}')
    logger.info(f'Average assembly standard deviation: {np.mean(stds):.3f}')

    out_path = join(cwd, 'data/trna_synt_stats.csv')
    logger.info(f'Exporting to {out_path}')
    res_df.to_csv(out_path, index=False)

    logger.info('DONE')


if __name__ == '__main__':
    main()
