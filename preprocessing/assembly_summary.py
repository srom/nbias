"""
Gather metadata for all assemblies in the dataset 
and output to assemblies.csv
"""
import argparse
import os
import logging
import re

import numpy as np
import pandas as pd


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s (%(levelname)s) %(message)s")

    data_path = os.path.join(os.getcwd(), 'data')
    sequences_path = os.path.join(data_path, 'sequences')
    ncbi_assembly_summary_path = os.path.join(data_path, 'NCBI_assembly_summary.txt')
    arc_metadata_path = os.path.join(data_path, 'archaea/DB_Archaea95_update012020.info.txt')
    bac_metadata_path = os.path.join(data_path, 'bacteria/DB_BACT95_HCLUST0.5/DB_BACT95_HCLUST0.5.info.txt')
    output_file = os.path.join(data_path, 'assemblies.csv')

    logger.info('Loading NCBI assembly summary')
    ncbi_summary_df = pd.read_csv(
        ncbi_assembly_summary_path, 
        sep='\t', 
        skiprows=1,
    ).set_index('assembly_accession')

    logger.info('Loading Archea metadata file')
    arc_metadata_df = read_metadata_file(arc_metadata_path)

    logger.info('Loading Bacteria metadata file')
    bac_metadata_df = read_metadata_file(bac_metadata_path)

    logger.info('Concatenating metadata files')
    metadata_df = pd.concat(
        [arc_metadata_df, bac_metadata_df],
        ignore_index=True,
    ).set_index('assembly_accession')

    logger.info('Merging metadata files')
    output_df = pd.merge(
        ncbi_summary_df[[
            'taxid',
            'species_taxid',
            'organism_name',
            'assembly_level',
        ]],
        metadata_df[[
            'domain',
            'phylum',
            'class',
            'order',
            'family',
            'genus',
            'species',
            'strain',
        ]],
        how='inner',
        on='assembly_accession',
    )

    output_df['taxid'] = pd.to_numeric(output_df['taxid'])
    output_df['species_taxid'] = pd.to_numeric(output_df['species_taxid'])

    output_df = output_df.reset_index(drop=False).sort_values(
        ['species_taxid', 'taxid', 'assembly_accession', 'organism_name']
    ).set_index('assembly_accession')

    logger.info('Check that no key attributes are missing')
    assert len(output_df[output_df['organism_name'].isnull()]) == 0
    assert len(output_df[output_df['taxid'].isnull()]) == 0
    assert len(output_df[output_df['species_taxid'].isnull()]) == 0
    assert len(output_df[output_df['domain'].isnull()]) == 0

    logger.info('Check that there are no duplicates')
    assert len(output_df.reset_index()['assembly_accession'].unique()) == len(output_df)

    logger.info('Filter out assemblies with missing files')
    assemblies_to_keep = []
    for i, assembly in enumerate(output_df.index):
        if i == 0 or (i+1) % 1000 == 0:
            logger.info(f'Checking assembly {i+1} / {len(output_df)}')

        genome_path = os.path.join(sequences_path, f'{assembly}/{assembly}_genomic.fna.gz')
        if not os.path.isfile(genome_path):
            continue

        cds_path = os.path.join(sequences_path, f'{assembly}/{assembly}_cds_from_genomic.fna.gz')
        if not os.path.isfile(cds_path):
            continue

        gff_path = os.path.join(sequences_path, f'{assembly}/{assembly}_genomic.gff.gz')
        if not os.path.isfile(gff_path):
            continue

        proteome_path = os.path.join(sequences_path, f'{assembly}/{assembly}_protein.faa.gz')
        if not os.path.isfile(proteome_path):
            continue

        assemblies_to_keep.append(assembly)

    logger.info(f'Writing output file with {len(assemblies_to_keep):,} rows')
    output_columns = [
        'taxid',
        'species_taxid',
        'organism_name',
        'domain',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'strain',
        'assembly_level',
    ]
    output_df.loc[assemblies_to_keep][output_columns].to_csv(output_file)

    logger.info('DONE')


def read_metadata_file(path):
    column_names = [
        'taxid',
        'organism_name',
        'assembly_accession',
        'refseq_category',
        'assembly_level',
        'unused_column_1',
        'unused_column_2',
        'taxonomy',
    ]
    df = pd.read_csv(path, sep='\t', header=None, names=column_names)
    classification_levels = [
        'domain',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'strain',
    ]

    def set_classification_level(ix):

        def fn(classification):
            levels = [
                l.strip() if l.strip() != '' and l.strip().lower() != 'na'
                else None
                for l in classification.split(';')
            ]
            try:
                return levels[i]
            except IndexError:
                return None

        return fn

    for i, c in enumerate(classification_levels):
        df[c] = df['taxonomy'].apply(set_classification_level(i))

    return df


if __name__ == '__main__':
    main()
