#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Combine partial motif or track vs regions or genes (features) cisTarget scores databases to a complete
               cisTarget scores database:
                 1) Transpose partial cisTarget motif or track vs regions or genes (features) scores databases
                    to partial cisTarget regions or genes (features) vs motif or tracks databases.
                 2) Combine those partial cisTarget regions or genes (features) vs motif or tracks scores databases
                    to a complete cisTarget regions or genes (features) vs motif or tracks database.
                 3) Transpose a complete cisTarget regions or genes (features) vs motif or tracks scores database
                    to a complete cisTarget motif or tracks vs region or genes scores database.

Copyright (C): 2020-2021 - Gert Hulselmans
"""


import argparse
import glob
import os
import re
import sys
import time

from typing import List

from cistarget_db import DatabaseTypes, CisTargetDatabase


partial_ct_scores_db_features_vs_motifs_or_tracks_pattern = re.compile(
    r'''^(?P<db_prefix_minimal>.*?)(\.(?P<species>(?!part[^.]+[0-9]+)[^.]+))?\.part_(?P<current_part>[0-9]+)_of_(?P<nbr_total_parts>[0-9]+)\.(?P<db_type>(?P<feature_type>regions|genes)_vs_(?P<motifs_or_tracks_type>motifs|tracks)\.scores\.feather)$'''
)


def combine_partial_cistarget_dbs(
        partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames: List[str],
        output_db_prefix: str):
    """
    Combine partial cisTarget databases to a full cisTarget database.

    :param partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames:
        partial cisTarget features vs motifs or tracks scores databases filenames (for one species)
    :param output_db_prefix:
        Output database prefixes used for writing complete cisTarget features vs motifs or tracks scores database and
        complete cisTarget motifs or tracks vs features scores database.
    :return:
    """

    start_reading_partial_dbs_to_complete_db_time = time.monotonic()
    nbr_partial_dbs = len(partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames)
    partial_db_type = DatabaseTypes.create_database_type_and_db_prefix_and_extension_from_db_filename(
        db_filename=partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames[0]
    )[0]
    print(
        f'\nReading cisTarget {nbr_partial_dbs} partial {partial_db_type.column_kind} vs {partial_db_type.row_kind} '
        f'scores databases to one complete database:\n  - "'
        + '"\n  - "'.join(partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames) + '"\n'
    )

    ct_scores_db_features_vs_motifs_or_tracks = CisTargetDatabase.read_db(
        db_filename_or_dbs_filenames=partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames
    )

    elapsed_reading_partial_dbs_to_complete_db_time = time.monotonic() - start_reading_partial_dbs_to_complete_db_time
    print(
        f'Reading cisTarget {nbr_partial_dbs} partial {partial_db_type.column_kind} vs {partial_db_type.row_kind} '
        f'scores databases to one complete database took: '
        f'{elapsed_reading_partial_dbs_to_complete_db_time:.06f} seconds\n'
    )

    start_writing_features_vs_motifs_or_tracks_scores_db_time = time.monotonic()
    ct_scores_db_features_vs_motifs_or_tracks_filename = \
        ct_scores_db_features_vs_motifs_or_tracks.create_db_filename_from_db_prefix(
            db_prefix=output_db_prefix,
            extension='feather'
        )
    print(
        f'Writing full cisTarget {ct_scores_db_features_vs_motifs_or_tracks.db_type.column_kind} vs '
        f'{ct_scores_db_features_vs_motifs_or_tracks.db_type.row_kind} scores db: '
        f'"{ct_scores_db_features_vs_motifs_or_tracks_filename}"'
    )
    
    ct_scores_db_features_vs_motifs_or_tracks.write_db(
        db_prefix=output_db_prefix,
        version=1
    )

    elapsed_writing_features_vs_motifs_or_tracks_scores_db_time = (
            time.monotonic() - start_writing_features_vs_motifs_or_tracks_scores_db_time
    )
    print(
        f'Writing full cisTarget {ct_scores_db_features_vs_motifs_or_tracks.db_type.column_kind} vs '
        f'{ct_scores_db_features_vs_motifs_or_tracks.db_type.row_kind} scores db took: '
        f'{elapsed_writing_features_vs_motifs_or_tracks_scores_db_time:.06f} seconds\n'
    )

    start_writing_motifs_or_tracks_vs_features_scores_db_time = time.monotonic()
    ct_scores_db_motifs_or_tracks_vs_features_filename = \
        DatabaseTypes.from_strings(
            scores_or_rankings='scores',
            column_kind=ct_scores_db_features_vs_motifs_or_tracks.db_type.row_kind,
            row_kind=ct_scores_db_features_vs_motifs_or_tracks.db_type.column_kind
        ).create_db_filename(
            db_prefix=output_db_prefix,
            extension='feather'
        )
    print(
        f'Writing full cisTarget {ct_scores_db_features_vs_motifs_or_tracks.db_type.row_kind} vs '
        f'{ct_scores_db_features_vs_motifs_or_tracks.db_type.column_kind} scores db: '
        f'"{ct_scores_db_motifs_or_tracks_vs_features_filename}"'
    )

    ct_scores_db_features_vs_motifs_or_tracks.transpose(order='F').write_db(
        db_prefix=output_db_prefix,
        version=1
    )

    elapsed_writing_motifs_or_tracks_vs_features_scores_db_time = (
            time.monotonic() - start_writing_motifs_or_tracks_vs_features_scores_db_time
    )
    print(
        f'Writing full cisTarget {ct_scores_db_features_vs_motifs_or_tracks.db_type.row_kind} vs '
        f'{ct_scores_db_features_vs_motifs_or_tracks.db_type.column_kind} scores db took: '
        f'{elapsed_writing_motifs_or_tracks_vs_features_scores_db_time:.06f} seconds\n'
    )


def main():
    parser = argparse.ArgumentParser(
        description='Combine partial cisTarget regions or genes (features) vs motifs or tracks scores databases to: '
                    '1) a complete cisTarget regions or genes (features) vs motifs or tracks scores database and '
                    '2) a complete cisTarget motifs or tracks vs regions or genes (features) scores database.'
    )

    parser.add_argument(
        '-i',
        '--input',
        dest='input',
        action='store',
        type=str,
        required=True,
        help='Input directory or database prefix with partial cisTarget regions or genes (features) vs motif or track '
             'scores database Feather files.'
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='output_dir',
        action='store',
        type=str,
        required=True,
        help='Output directory to which the 1) complete cisTarget regions or genes (features) vs motif or track scores '
             'database Feather files and 2) complete cisTarget motifs or tracks vs regions or genes (features) scores '
             'database Feather files will be written.'
    )

    args = parser.parse_args()

    # Construct glob string to find '*.part*_*_of_*.*_vs_*.scores.feather' cisTarget databases.
    partial_ct_scores_db_features_vs_motifs_or_tracks_glob_str = (
        os.path.join(args.input, '*.part*_*_of_*.*_vs_*.scores.feather')
        if os.path.isdir(args.input)
        else args.input + '*.part*_*_of_*.*_vs_*.scores.feather'
    )

    # Get all "*.*_vs_*.scores.feather" cisTarget databases in the input directory or the ones that start with the
    # requested prefix and sort them.
    partial_ct_scores_db_features_vs_motifs_or_tracks_vs_filenames = sorted(
        glob.glob(partial_ct_scores_db_features_vs_motifs_or_tracks_glob_str)
    )

    if len(partial_ct_scores_db_features_vs_motifs_or_tracks_vs_filenames) == 0:
        print(
            f'Error: No partial cisTarget regions or genes (features) vs motifs or tracks scores databases found '
            f'matching glob: "{partial_ct_scores_db_features_vs_motifs_or_tracks_glob_str}".',
            file=sys.stderr
        )
        sys.exit(1)

    ct_dbs_hierarchical_dict = dict()

    for partial_ct_scores_db_features_vs_motifs_or_tracks_filename in partial_ct_scores_db_features_vs_motifs_or_tracks_vs_filenames:
        # Parse partial cisTarget features vs motifs or tracks scores database filename, so partial databases made from
        # the same regions/genes and/or species but different motifs/tracks that belong together can be grouped and
        # combined to one database in a later step.
        match = partial_ct_scores_db_features_vs_motifs_or_tracks_pattern.match(
            partial_ct_scores_db_features_vs_motifs_or_tracks_filename
        )

        if match:
            db_prefix_minimal = match.group('db_prefix_minimal')
            db_type = match.group('db_type')
            # Specify "no_species" if species was not specified in cisTarget database filename.
            species = match.group('species') if match.group('species') else 'no_species'
            current_part = int(match.group('current_part'))
            nbr_total_parts = int(match.group('nbr_total_parts'))

            # Group databases in a hierarchical structure:
            #   - db_prefix_minimal
            #   - db_type
            #   - species
            #   - current part number, total number of parts
            ct_dbs_hierarchical_dict.setdefault(db_prefix_minimal, dict())
            ct_dbs_hierarchical_dict[db_prefix_minimal].setdefault(db_type, dict())
            ct_dbs_hierarchical_dict[db_prefix_minimal][db_type].setdefault(species, dict())
            ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species].setdefault(
                (current_part, nbr_total_parts), dict()
            )
            ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species][(current_part, nbr_total_parts)] = \
                partial_ct_scores_db_features_vs_motifs_or_tracks_filename

    for db_prefix_minimal in ct_dbs_hierarchical_dict:
        for db_type in ct_dbs_hierarchical_dict[db_prefix_minimal]:
            for species in ct_dbs_hierarchical_dict[db_prefix_minimal][db_type]:
                # Get all partial cisTarget features vs motifs or tracks scores database filenames for the current
                # species and sort them by current part and nbr of total parts number, so when reading them to
                # a complete cisTarget scores database, the motifs/tracks are already in the correct order.
                partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames = [
                    ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species][(current_part, nbr_total_parts)]
                    for current_part, nbr_total_parts in sorted(
                        ct_dbs_hierarchical_dict[db_prefix_minimal][db_type][species]
                    )
                ]

                # Construct database prefix name for the complete cisTarget scores database.
                output_db_prefix = f'{os.path.join(args.output_dir, os.path.basename(db_prefix_minimal))}' \
                                   f'{"." + str(species) if species != "no_species" else ""}'

                combine_partial_cistarget_dbs(
                    partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames=partial_ct_scores_db_features_vs_motifs_or_tracks_vs_one_species_filenames,
                    output_db_prefix=output_db_prefix
                )


if __name__ == '__main__':
    main()
