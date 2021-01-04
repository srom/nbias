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
    parser.add_argument('source', type=str, choices=['prodigal', 'ncbi'])
    args = parser.parse_args()

    source = args.source

    output_base_path = os.path.join(os.getcwd(), 'data/sequences')
    metadata_output_path = os.path.join(os.getcwd(), 'data/assemblies.csv')

    output_metadata = {
        'assembly_accession': [],
        'asm_name': [],
        'source': [],
        'n_cds': [],
    }

    logger.info(f'Extracting CDS from {source} annotations')

    base_folder = os.path.join(os.getcwd(), 'data/bacteria')

    genomes_folder = os.path.join(base_folder, 'GENOMES_FASTA')
    if source == 'prodigal':
        gff_folder = genomes_folder
    else:
        gff_folder = os.path.join(base_folder, 'DB_BACT95_HCLUST0.5/GENOMES_GFF')

    metadata_list = get_annotation_paths(gff_folder)

    logger.info(f'Number of gff files found for source {source}: {len(metadata_list):,}')

    for i, annotation_metadata in enumerate(metadata_list):
        if i == 0 or (i+1) % 100 == 0:
            logger.info(f'Processing {source} annotations {i+1:,} / {len(metadata_list):,}')

        assembly = annotation_metadata.assembly
        version = annotation_metadata.version
        gff_path = str(annotation_metadata.path)

        output_folder = os.path.join(output_base_path, f'{assembly}')
        Path(output_folder).mkdir(parents=False, exist_ok=True)

        output_path = os.path.join(output_folder, f'{assembly}_cds_from_genomic.fna.gz')
        genome_path = os.path.join(genomes_folder, f'{assembly}_{version}_genomic.fna')

        cds_records = extract_cds_records(assembly, version, source, gff_path, genome_path)

        if len(cds_records) == 0:
            logger.warning(f'No CDS for assembly {assembly}_{version}')
            continue

        output_metadata['assembly_accession'].append(assembly)
        output_metadata['asm_name'].append(version.replace('_', ' '))
        output_metadata['source'].append(source)
        output_metadata['n_cds'].append(len(cds_records))

        with gzip.open(output_path, 'wt') as f_gz:
            SeqIO.write(cds_records, f_gz, 'fasta')

    logger.info(f'Exporting metadata to {metadata_output_path}')
    pd.DataFrame(output_metadata).to_csv(metadata_output_path, index=False)

    logger.info('DONE')


def extract_cds_records(assembly, version, source, gff_path, genome_path):
    with open(genome_path) as f:
        chromosome_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

    annotations = gffutils.create_db(gff_path, ':memory:', merge_strategy='replace')

    records = []
    for f in annotations.features_of_type('CDS'): 
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
        suffix = feature_id.split('_')[-1]
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


class AnnotationMetadata(object):

    def __init__(self, path):
        self.path = path
        self.assembly, self.version = parse_assembly(self.path)


if __name__ == '__main__':
    main()
