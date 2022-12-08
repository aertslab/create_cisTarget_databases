#!/usr/bin/env python3

"""
Purpose :      Create cisTarget cross-species motif rankings databases.

Copyright (C): 2020-2021 - Gert Hulselmans
"""


import argparse
import glob
import os
import re
import sys
import time
from typing import List

from cistarget_db import CisTargetDatabase, DatabaseTypes

ct_rankings_db_motifs_vs_regions_or_genes_species_pattern = re.compile(
    r'''^(?P<db_prefix_minimal>.*?)(\.(?P<species>(?!(part[^.]+[0-9]+)|cross_?species)[^.]+))\.(?P<db_type>(?P<motifs_or_tracks_type>motifs)_vs_(?P<regions_or_genes_type>regions|genes)\.rankings\.feather)$'''
)


def create_cross_species_motifs_rankings_db(cross_species_db_prefix: str,
                                            ct_rankings_db_motifs_vs_regions_or_genes_species_filenames: List[str]):

    cross_species_db_prefix = (
        cross_species_db_prefix
        if cross_species_db_prefix.endswith('.cross_species')
        else cross_species_db_prefix + '.cross_species'
    )

    start_creating_ct_cross_species_rankings_db_motifs_vs_regions_or_genes_time = time.monotonic()
    print(
        'Create cisTarget cross-species motifs rankings db from:\n  - "'
        + '"\n  - "'.join(ct_rankings_db_motifs_vs_regions_or_genes_species_filenames) + '"\n'
    )

    ct_cross_species_rankings_db_motifs_vs_regions_or_genes = CisTargetDatabase.create_cross_species_rankings_db(
        species_rankings_db_filenames=ct_rankings_db_motifs_vs_regions_or_genes_species_filenames
    )

    elapsed_creating_ct_cross_species_rankings_db_motifs_vs_regions_or_genes_time = (
            time.monotonic() - start_creating_ct_cross_species_rankings_db_motifs_vs_regions_or_genes_time
    )
    print(
        f'Creating cisTarget cross-species '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.column_kind} vs '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.row_kind} rankings db took: '
        f'{elapsed_creating_ct_cross_species_rankings_db_motifs_vs_regions_or_genes_time:.06f} seconds'
    )

    start_writing_ct_cross_species_rankings_db_motifs_vs_regions_or_genes_time = time.monotonic()
    ct_cross_species_rankings_db_motifs_vs_regions_or_genes_filename = (
        ct_cross_species_rankings_db_motifs_vs_regions_or_genes.create_db_filename_from_db_prefix(
            db_prefix=cross_species_db_prefix,
            extension='feather'
        )
    )
    print(
        f'Writing cisTarget cross-species '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.column_kind} vs '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.row_kind} rankings db: '
        f'"{ct_cross_species_rankings_db_motifs_vs_regions_or_genes_filename}"'
    )

    ct_cross_species_rankings_db_motifs_vs_regions_or_genes.write_db(
        db_prefix=cross_species_db_prefix,
        version=2,
    )

    elapsed_writing_ct_cross_species_rankings_db_motifs_vs_regions_or_genes_time = (
            time.monotonic() - start_writing_ct_cross_species_rankings_db_motifs_vs_regions_or_genes_time
    )
    print(
        f'Writing cisTarget cross-species '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.column_kind} '
        f'vs {ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.row_kind} rankings db '
        f'took: {elapsed_writing_ct_cross_species_rankings_db_motifs_vs_regions_or_genes_time:.06f} seconds\n'
    )

    start_writing_ct_cross_species_rankings_db_regions_or_genes_vs_motifs_time = time.monotonic()
    ct_cross_species_rankings_db_regions_or_genes_vs_motifs_filename = \
        DatabaseTypes.from_strings(
            scores_or_rankings='rankings',
            column_kind=ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.row_kind,
            row_kind=ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.column_kind
        ).create_db_filename(
            db_prefix=cross_species_db_prefix,
            extension='feather'
        )
    print(
        f'Writing cisTarget cross-species '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.row_kind} vs '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.column_kind} rankings db: '
        f'"{ct_cross_species_rankings_db_regions_or_genes_vs_motifs_filename}"'
    )

    ct_cross_species_rankings_db_motifs_vs_regions_or_genes.transpose().write_db(
        db_prefix=cross_species_db_prefix,
        version=2,
    )

    elapsed_writing_ct_cross_species_rankings_db_regions_or_genes_vs_motifs_time = (
            time.monotonic() - start_writing_ct_cross_species_rankings_db_regions_or_genes_vs_motifs_time
    )
    print(
        f'Writing cisTarget cross-species '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.row_kind} vs '
        f'{ct_cross_species_rankings_db_motifs_vs_regions_or_genes.db_type.column_kind} rankings db '
        f'took: {elapsed_writing_ct_cross_species_rankings_db_regions_or_genes_vs_motifs_time:.06f} seconds\n'
    )


def main():
    parser = argparse.ArgumentParser(
        description='Create cisTarget cross-species motifs rankings databases.'
    )

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        action='store',
        type=str,
        required=True,
        help='Input directory or database prefix with cisTarget motifs vs regions or genes rankings databases per '
             'species.'
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='output_dir',
        action='store',
        type=str,
        required=True,
        help='Output directory to which the cisTarget cross-species motifs rankings database files will be written.'
    )

    args = parser.parse_args()

    # Construct glob string to find "*.*.motifs_vs_*.rankings.feather" cisTarget rankings database per species.
    ct_rankings_db_motifs_vs_regions_or_genes_species_glob_str = (
        os.path.join(args.input, '*.*.motifs_vs_*.rankings.feather')
        if os.path.isdir(args.input)
        else args.input + '*.*.motifs_vs_*.rankings.feather'
    )

    # Get all "*.*.motifs_vs_*.rankings.feather" cisTarget databases in the input directory or the ones that start with
    # the requested prefix.
    ct_rankings_db_motifs_vs_regions_or_genes_species_filenames = sorted(
        glob.glob(ct_rankings_db_motifs_vs_regions_or_genes_species_glob_str)
    )

    if len(ct_rankings_db_motifs_vs_regions_or_genes_species_filenames) == 0:
        print(
            f'Error: No cisTarget motifs vs regions or genes rankings databases per species found matching '
            f'glob: "{ct_rankings_db_motifs_vs_regions_or_genes_species_glob_str}".',
            file=sys.stderr
        )
        sys.exit(1)

    ct_dbs_hierarchical_dict = dict()

    for ct_rankings_db_motifs_vs_regions_or_genes_species_filename in ct_rankings_db_motifs_vs_regions_or_genes_species_filenames:
        # Parse cisTarget database filename, so databases made from the same regions/genes, but lifted over
        # to different species can be grouped together.
        match = ct_rankings_db_motifs_vs_regions_or_genes_species_pattern.match(
            ct_rankings_db_motifs_vs_regions_or_genes_species_filename
        )

        if match:
            db_prefix_minimal = match.group('db_prefix_minimal')
            db_type = match.group('db_type')
            species = match.group('species')

            ct_dbs_hierarchical_dict.setdefault(db_prefix_minimal, dict())
            ct_dbs_hierarchical_dict[db_prefix_minimal].setdefault(db_type, dict())
            ct_dbs_hierarchical_dict[db_prefix_minimal][db_type].setdefault(species, dict())
            ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species] = ct_rankings_db_motifs_vs_regions_or_genes_species_filename

    for db_prefix_minimal in ct_dbs_hierarchical_dict:
        for db_type in ct_dbs_hierarchical_dict[db_prefix_minimal]:
            # Get all cisTarget motifs vs regions or genes species rankings databases for the same regions/genes.
            ct_rankings_db_motifs_vs_regions_or_genes_species_filenames = [
                ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species]
                for species in sorted(ct_dbs_hierarchical_dict[db_prefix_minimal][db_type])
            ]

            create_cross_species_motifs_rankings_db(
                cross_species_db_prefix=db_prefix_minimal,
                ct_rankings_db_motifs_vs_regions_or_genes_species_filenames=ct_rankings_db_motifs_vs_regions_or_genes_species_filenames
            )


if __name__ == '__main__':
    main()
