"""
Extract CDS from GFF annotations
"""
import argparse
import os
import logging
import re
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils


logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s (%(levelname)s) %(message)s")

    parser = argparse.ArgumentParser()
    parser.add_argument('--pattern', type=str, default='archaea/genomes_assemblies_gff/*.gff')
    args = parser.parse_args()

    pattern = args.pattern

    output_base_path = os.path.join(os.getcwd(), 'data/sequences')

    logger.info(f'Extracting CDS for pattern: {pattern}')

    base_folder = os.path.join(os.getcwd(), 'data')

    paths = list(Path(base_folder).glob(pattern))

    logger.info(f'Number of gff files found: {len(paths):,}')

    for i, path in enumerate(paths):
        if i == 0 or (i+1) % 100 == 0:
            logger.info(f'Processing annotations {i+1:,} / {len(paths):,}')

        try:
            assembly, version = parse_assembly(path)
        except ValueError as e:
            logger.warning(e)
            continue

        gff_path = str(path)

        output_folder = os.path.join(output_base_path, f'{assembly}')
        Path(output_folder).mkdir(parents=False, exist_ok=True)

        output_path = os.path.join(output_folder, f'{assembly}_cds_from_genomic.fna.gz')
        genome_path = os.path.join(
            base_folder,
            f'archaea/genomes_assemblies_DNA_fasta/ncbi-genomes-2020-05-15/{assembly}_{version}_genomic.fna',
        )

        cds_records = extract_cds_records(assembly, version, gff_path, genome_path)

        if len(cds_records) == 0:
            logger.warning(f'No CDS for assembly {assembly}_{version}')
            continue

        with gzip.open(output_path, 'wt') as f_gz:
            SeqIO.write(cds_records, f_gz, 'fasta')

    logger.info('DONE')


def extract_cds_records(assembly, version, gff_path, genome_path):
    with open(genome_path) as f:
        chromosome_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

    source = None
    annotations = gffutils.create_db(gff_path, ':memory:', merge_strategy='replace')

    records = []
    for f in annotations.features_of_type('CDS'):
        if source is None:
            source = 'prodigal' if 'prodigal' in f.source.lower().strip() else 'ncbi'

        start_ix = f.start - 1
        end_ix = f.end

        seq = chromosome_dict[f.seqid].seq

        if f.strand == '+':
            cds_seq = seq[start_ix:end_ix]
        else:
            cds_seq = seq[start_ix:end_ix].reverse_complement()
            
        cds_id = make_id(
            f.id, 
            f.seqid,
            f.attributes, 
            source,
        )
        description = make_description(
            assembly, 
            version, 
            source, 
            f.seqid, 
            f.start, 
            f.end, 
            f.strand,
        )
        seq_record = SeqRecord(
            cds_seq,
            id=cds_id,
            description=description,
        )
        records.append(seq_record)

    return records


def make_id(feature_id, chromosome_id, attributes, source):
    if source == 'prodigal':
        suffix = feature_id.replace('gene-', '').split('_')[-1]
        return f'{chromosome_id}_{suffix}'
    else:
        if 'protein_id' in attributes and len(attributes['protein_id']) > 0:
            return attributes['protein_id'][0]
        elif feature_id.startswith('cds-'):
            return feature_id.replace('cds-', '')
        elif feature_id.startswith('cds_'):
            return feature_id.replace('cds-', '')
        else:
            return feature_id


def make_description(
    assembly, 
    version, 
    source, 
    chromosome_id, 
    start, 
    end, 
    strand,
):
    return ';'.join([
        assembly, 
        version, 
        source, 
        chromosome_id, 
        str(start), 
        str(end), 
        strand,
    ])


def get_annotation_paths(gff_folder):
    metadata_list = []
    for path in Path(gff_folder).glob('*.gff'):
        try:
            annotation_metadata = AnnotationMetadata(path)
        except ValueError as e:
            logger.warning(e)
            continue

        metadata_list.append(annotation_metadata)

    return metadata_list


def parse_assembly(path):
    m = re.match(r'^(GCA_[^_]+)_(.+)_genomic.*\.gff$', path.name)
    if m is not None:
        assembly = m[1]
        version = m[2]
        return assembly, version
    else:
        raise ValueError(f'No assembly found in gff filename: {path.name}')


if __name__ == '__main__':
    main()
