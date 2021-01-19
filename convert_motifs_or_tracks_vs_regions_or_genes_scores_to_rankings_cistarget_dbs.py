#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Convert motifs or track vs regions or genes (features) scores cisTarget database to
               cisTarget rankings database.

Copyright (C): 2020-2021 - Gert Hulselmans
"""


import argparse
import random
import sys
import time

from cistarget_db import DatabaseTypes, CisTargetDatabase


def main():
    parser = argparse.ArgumentParser(
        description='Convert motifs or tracks vs regions or genes (features) scores cisTarget database to cisTarget '
                    'rankings database.'
    )

    parser.add_argument(
        '-i',
        '--db',
        dest='ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename',
        action='store',
        type=str,
        required=True,
        help='Motifs or tracks vs regions or genes (features) scores cisTarget database filename. The rankings '
             'cisTarget database feather file will be written to the same directory.'
    )

    parser.add_argument(
        '-s',
        '--seed',
        dest='seed',
        action='store',
        type=int,
        required=False,
        help='Random seed used for breaking ties when creating rankings for a range of tied scores. '
             'When setting this seed to a specific value and running this script with the same input, will result in '
             'the same cisTarget rankings databases as output.'
    )

    args = parser.parse_args()

    # Get cisTarget db type, db prefix and extension.
    db_type, db_prefix, extension = DatabaseTypes.create_database_type_and_db_prefix_and_extension_from_db_filename(
        db_filename=args.ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename
    )

    if db_type not in {DatabaseTypes.SCORES_DB_MOTIFS_VS_REGIONS,
                       DatabaseTypes.SCORES_DB_MOTIFS_VS_GENES,
                       DatabaseTypes.SCORES_DB_TRACKS_VS_REGIONS,
                       DatabaseTypes.SCORES_DB_TRACKS_VS_GENES}:
        print(
            f'Error: cisTarget database "{args.ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename}" is not a '
            'motifs or tracks vs regions or genes (features) scores cisTarget database.',
            file=sys.stderr
        )
        sys.exit(1)

    def write_db(ct_db: CisTargetDatabase, db_prefix: str):
        """
        Write cisTarget database to a Feather file and print database location and elapsed time.

        :param ct_db: cisTarget database object.
        :param db_prefix: Feather database file prefix.
        :return:
        """
        db_filename = ct_db.create_db_filename_from_db_prefix(db_prefix=db_prefix, extension='feather')

        print(
            f'Write {ct_db.db_type.row_kind} vs {ct_db.db_type.column_kind} {ct_db.db_type.scores_or_rankings} db: '
            f'{db_filename}.'
        )

        start_time = time.monotonic()
        ct_db.write_db(
            db_prefix=db_prefix,
            version=1
        )
        elapsed_time = time.monotonic() - start_time

        print(
            f'Writing {ct_db.db_type.row_kind} vs {ct_db.db_type.column_kind} {ct_db.db_type.scores_or_rankings} db '
            f'took: {elapsed_time} seconds.\n'
        )

    # Create cisTarget rankings database (motifs or tracks vs regions or genes) from cisTarget scores database filename
    # (motifs or tracks vs regions or genes).
    start_reading_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_db_time = time.monotonic()
    print(

        f'\nReading {db_type.column_kind} vs {db_type.row_kind} cisTarget '
        f'scores database "{args.ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename}".'
    )

    ct_scores_db_motifs_or_tracks_vs_regions_or_genes = CisTargetDatabase.read_db(
        db_filename_or_dbs_filenames=args.ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename
    )

    elapsed_reading_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_db_time = \
        time.monotonic() - start_reading_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_db_time
    print(
        f'Reading {db_type.column_kind} vs {db_type.row_kind} cisTarget scores databases took: '
        f'{elapsed_reading_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_db_time} seconds.\n'
    )

    # Set random seed to provided input value or a random integer.
    seed = args.seed if args.seed else random.randint(0, 2 ** 64)

    # Create cisTarget rankings database (motifs or tracks vs regions or genes) from cisTarget scores database filename
    # (motifs or tracks vs regions or genes).
    start_creating_rankings_ct_rankings_db_motifs_or_tracks_vs_regions_or_genes_db_time = time.monotonic()
    print(
        f'Create rankings from "{args.ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename}" with random seed '
        f'set to {seed}.',
        file=sys.stderr
    )

    ct_rankings_db_motifs_or_tracks_vs_regions_or_genes = \
        ct_scores_db_motifs_or_tracks_vs_regions_or_genes.convert_scores_db_to_rankings_db(seed=seed)

    elapsed_creating_rankings_ct_rankings_db_motifs_or_tracks_vs_regions_or_genes_db_time = \
        time.monotonic() - start_creating_rankings_ct_rankings_db_motifs_or_tracks_vs_regions_or_genes_db_time
    print(
        f'Creating rankings from scores database took '
        f'{elapsed_creating_rankings_ct_rankings_db_motifs_or_tracks_vs_regions_or_genes_db_time:.06f} seconds.\n',
        file=sys.stderr
    )

    # Reclaim memory occupied by motifs or tracks vs regions or genes cisTarget scores database.
    del ct_scores_db_motifs_or_tracks_vs_regions_or_genes

    # Write cisTarget rankings database (motifs or tracks vs regions or genes) to Feather file.
    write_db(ct_db=ct_rankings_db_motifs_or_tracks_vs_regions_or_genes, db_prefix=db_prefix)

    # Create cisTarget rankings database (regions or genes vs motifs or tracks) from (motifs or tracks vs regions or
    # genes) version.
    print(
        f'Convert {db_type.column_kind} vs {db_type.row_kind} cisTarget rankings database to {db_type.row_kind} vs '
        f'{db_type.column_kind} cisTarget rankings databases.',
        file=sys.stderr
    )
    ct_rankings_db_regions_or_genes_vs_motifs_or_tracks = ct_rankings_db_motifs_or_tracks_vs_regions_or_genes.transpose()

    # Write cisTarget rankings database (regions or genes vs motifs or tracks) to Feather file.
    write_db(ct_db=ct_rankings_db_regions_or_genes_vs_motifs_or_tracks, db_prefix=db_prefix)


if __name__ == '__main__':
    main()
