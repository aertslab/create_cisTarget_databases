#!/usr/bin/env python3

"""
Purpose :      Combine partial motif or track vs regions or genes cisTarget scores databases to a complete
               cisTarget scores database:
                 1) Combine partial cisTarget motif or track vs regions or genes scores databases
                    to a complete cisTarget motif or tracks vs regions or genes database.
                 2) Transpose a complete cisTarget motif or tracks vs regions or genes scores database
                    to a complete cisTarget region or genes vs motif or tracks scores database.

Copyright (C): 2022 - Gert Hulselmans
"""


import argparse
import glob
import os
import re
import sys
import time
from typing import List

from cistarget_db import CisTargetDatabase, DatabaseTypes

partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_pattern = re.compile(
    r'''^(?P<db_prefix_minimal>.*?)(\.(?P<species>(?!part[^.]+[0-9]+)[^.]+))?\.part_(?P<current_part>[0-9]+)_of_(?P<nbr_total_parts>[0-9]+)(\.(?P<min_max_motifs>min_[0-9]+_to_max_([0-9]+_)?motifs))?\.(?P<db_type>(?P<motifs_or_tracks_type>motifs|tracks)_vs_(?P<regions_or_genes_type>regions|genes)\.scores\.feather)$'''
)


def combine_partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_cistarget_dbs(
        partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames: List[str],
        output_db_prefix: str):
    """
    Combine partial cisTarget motif or track vs regions or genes score databases to a full cisTarget database.

    :param partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames:
        partial cisTarget motifs or tracks vs regions or genes scores databases filenames (for one species)
    :param output_db_prefix:
        Output database prefixes used for writing complete cisTarget motifs or tracks vs regions or genes scores
        database and complete cisTarget regions or genes vs motifs or tracks scores database.
    :return:
    """

    start_reading_partial_dbs_to_complete_db_time = time.monotonic()
    nbr_partial_dbs = len(partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames)
    partial_db_type = DatabaseTypes.create_database_type_and_db_prefix_and_extension_from_db_filename(
        db_filename=partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames[0]
    )[0]
    print(
        f'\nReading cisTarget {nbr_partial_dbs} partial {partial_db_type.column_kind} vs {partial_db_type.row_kind} '
        f'scores databases to one complete database:\n  - "'
        + '"\n  - "'.join(partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames) + '"\n'
    )

    ct_scores_db_motifs_or_tracks_vs_regions_or_genes = CisTargetDatabase.read_db(
        db_filename_or_dbs_filenames=partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames
    )

    elapsed_reading_partial_dbs_to_complete_db_time = time.monotonic() - start_reading_partial_dbs_to_complete_db_time
    print(
        f'Reading cisTarget {nbr_partial_dbs} partial {partial_db_type.column_kind} vs {partial_db_type.row_kind} '
        f'scores databases to one complete database took: '
        f'{elapsed_reading_partial_dbs_to_complete_db_time:.06f} seconds\n'
    )

    start_writing_motifs_or_tracks_vs_regions_or_genes__scores_db_time = time.monotonic()
    ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename = \
        ct_scores_db_motifs_or_tracks_vs_regions_or_genes.create_db_filename_from_db_prefix(
            db_prefix=output_db_prefix,
            extension='feather'
        )
    print(
        f'Writing full cisTarget {ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.column_kind} vs '
        f'{ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.row_kind} scores db: '
        f'"{ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename}"'
    )
    
    ct_scores_db_motifs_or_tracks_vs_regions_or_genes.write_db(
        db_prefix=output_db_prefix,
        version=2,
    )

    elapsed_writing_motifs_or_tracks_vs_regions_or_genes_scores_db_time = (
            time.monotonic() - start_writing_motifs_or_tracks_vs_regions_or_genes__scores_db_time
    )
    print(
        f'Writing full cisTarget {ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.column_kind} vs '
        f'{ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.row_kind} scores db took: '
        f'{elapsed_writing_motifs_or_tracks_vs_regions_or_genes_scores_db_time:.06f} seconds\n'
    )

    start_writing_motifs_or_tracks_vs_regions_or_genes_scores_db_time = time.monotonic()
    ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename = \
        DatabaseTypes.from_strings(
            scores_or_rankings='scores',
            column_kind=ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.row_kind,
            row_kind=ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.column_kind
        ).create_db_filename(
            db_prefix=output_db_prefix,
            extension='feather'
        )
    print(
        f'Writing full cisTarget {ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.row_kind} vs '
        f'{ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.column_kind} scores db: '
        f'"{ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename}"'
    )

    ct_scores_db_motifs_or_tracks_vs_regions_or_genes.transpose(order='F').write_db(
        db_prefix=output_db_prefix,
        version=2,
    )

    elapsed_writing_motifs_or_tracks_vs_regions_or_genes_scores_db_time = (
            time.monotonic() - start_writing_motifs_or_tracks_vs_regions_or_genes_scores_db_time
    )
    print(
        f'Writing full cisTarget {ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.row_kind} vs '
        f'{ct_scores_db_motifs_or_tracks_vs_regions_or_genes.db_type.column_kind} scores db took: '
        f'{elapsed_writing_motifs_or_tracks_vs_regions_or_genes_scores_db_time:.06f} seconds\n'
    )


def main():
    parser = argparse.ArgumentParser(
        description='Combine partial cisTarget motifs or tracks vs regions or genes scores databases to: '
                    '1) a complete cisTarget motifs or tracks vs regions or genes scores database and'
                    '2) a complete cisTarget regions or genes vs motifs or tracks scores database.'
    )

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        action='store',
        type=str,
        required=True,
        help='Input directory or database prefix with partial cisTarget motif or track vs regions or genes '
             'scores database Feather files.'
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='output_dir',
        action='store',
        type=str,
        required=True,
        help='Output directory to which the 1) complete cisTarget motifs or tracks vs regions or genes scores '
             'database Feather files and 2) complete cisTarget regions or genes vs motif or track scores '
             'database Feather files will be written.'
    )

    args = parser.parse_args()

    # Construct glob string to find "*.part*_*_of_*.*_vs_*.scores.feather" cisTarget databases.
    partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_glob_str = (
        os.path.join(args.input, '*.part*_*_of_*.*_vs_*.scores.feather')
        if os.path.isdir(args.input)
        else args.input + '*.part*_*_of_*.*_vs_*.scores.feather'
    )

    # Get all "*.*_vs_*.scores.feather" cisTarget databases in the input directory or the ones that start with the
    # requested prefix and sort them.
    partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_filenames = sorted(
        glob.glob(partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_glob_str)
    )

    # Construct glob string to find "*.part*_*_of_*.min_*_to_max_*motifs.*_vs_*.scores.feather" cisTarget databases.
    partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_with_min_max_motifs_glob_str = (
        os.path.join(args.input, '*.part*_*_of_*.min_*_to_max_*motifs.*_vs_*.scores.feather')
        if os.path.isdir(args.input)
        else args.input + '*.part*_*_of_*.min_*_to_max_*motifs.*_vs_*.scores.feather'
    )

    # Get all "*.part*_*_of_*.min_*_to_max_*motifs.*_vs_*.scores.feather" cisTarget databases in the input directory or
    # the ones that start with the requested prefix and sort them.
    partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_filenames.extend(
        sorted(
            glob.glob(partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_with_min_max_motifs_glob_str)
        )
    )

    if len(partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_filenames) == 0:
        print(
            f'Error: No partial cisTarget motifs or tracks or regions or genes scores databases found '
            f'matching glob: "{partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_glob_str}" '
            f'or "{partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_with_min_max_motifs_glob_str}".',
            file=sys.stderr
        )
        sys.exit(1)

    ct_dbs_hierarchical_dict = dict()

    for partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename in partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_filenames:
        # Parse partial cisTarget motifs or tracks vs regions or genes scores database filename, so partial databases
        # made from the same regions/genes and/or species but different motifs/tracks that belong together can be
        # grouped and combined to one database in a later step.
        match = partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_pattern.match(
            partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename
        )

        if match:
            db_prefix_minimal = match.group('db_prefix_minimal')
            db_type = match.group('db_type')
            # Specify "no_species" if species was not specified in cisTarget database filename.
            species = match.group('species') if match.group('species') else 'no_species'
            current_part = int(match.group('current_part'))
            nbr_total_parts = int(match.group('nbr_total_parts'))
            min_max_motifs = match.group('min_max_motifs') if match.group('min_max_motifs') else 'no_min_max_motifs'

            # Group databases in a hierarchical structure:
            #   - db_prefix_minimal
            #   - db_type
            #   - species
            #   - min and max number of motifs per Cluster-Buster file, current part number, total number of parts
            ct_dbs_hierarchical_dict.setdefault(db_prefix_minimal, dict())
            ct_dbs_hierarchical_dict[db_prefix_minimal].setdefault(db_type, dict())
            ct_dbs_hierarchical_dict[db_prefix_minimal][db_type].setdefault(species, dict())

            ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species][
                (min_max_motifs, current_part, nbr_total_parts)
            ] = partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_filename

    for db_prefix_minimal in ct_dbs_hierarchical_dict:
        for db_type in ct_dbs_hierarchical_dict[db_prefix_minimal]:
            for species in ct_dbs_hierarchical_dict[db_prefix_minimal][db_type]:
                # Get all partial cisTarget motifs or tracks vs regions or genes scores database filenames for the
                # current species and sort them by min and max number of motifs per Cluster-Buster file, current part
                # and nbr of total parts number.
                partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames = [
                    ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species][
                        (min_max_motifs, current_part, nbr_total_parts)
                    ]
                    for min_max_motifs, current_part, nbr_total_parts in sorted(
                        ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species]
                    )
                ]

                # Construct database prefix name for the complete cisTarget scores database.
                output_db_prefix = f'{os.path.join(args.output_dir, os.path.basename(db_prefix_minimal))}' \
                                   f'{"." + str(species) if species != "no_species" else ""}'

                combine_partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_cistarget_dbs(
                    partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames=
                    partial_ct_scores_db_motifs_or_tracks_vs_regions_or_genes_vs_one_species_filenames,
                    output_db_prefix=output_db_prefix
                )


if __name__ == '__main__':
    main()
