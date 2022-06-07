#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Create cisTarget motif databases.

Copyright (C): 2019-2022 - Gert Hulselmans
"""


import argparse
import multiprocessing as mp
import os
import random
import shutil
import sys
import time

from typing import Tuple

import pandas as pd

from cistarget_db import MotifsOrTracksType, RegionOrGeneIDs, MotifOrTrackIDs, DatabaseTypes, CisTargetDatabase
from clusterbuster import get_motif_id_to_filename_and_nbr_motifs_dict, run_cluster_buster_for_motif


def main():
    parser = argparse.ArgumentParser(
        description='Create cisTarget motif databases.'
    )

    parser.add_argument(
        '-f',
        '--fasta',
        dest='fasta_filename',
        action='store',
        type=str,
        required=True,
        help='FASTA filename which contains the regions/genes to score with Cluster-Buster for each motif. When '
             'creating a cisTarget species database from regions/genes lifted over from a different species, provide '
             'the original FASTA file for that species to -F.'
    )

    parser.add_argument(
        '-F',
        '--fasta-original-species',
        dest='original_species_fasta_filename',
        action='store',
        type=str,
        required=False,
        help='FASTA filename which contains all the regions/genes of the original species. The fasta file provided to '
             '-f can contain less regions (not all regions could be lifted over) than the one provided to -F, but to '
             'create a cisTarget cross-species database later, all individual cisTarget species databases need to '
             'contain the same amount of regions/genes.'
    )

    parser.add_argument(
        '-M',
        '--motifs_dir',
        dest='motifs_dir',
        action='store',
        type=str,
        required=True,
        help='Path to directory with Cluster-Buster motifs.'
    )

    parser.add_argument(
        '-m',
        '--motifs',
        dest='motifs_list_filename',
        action='store',
        type=str,
        required=True,
        help='Filename with list of motif IDs or motif MD5 names to be scored from directory specified by '
             '"--motifs_dir".'
    )

    parser.add_argument(
        '-5',
        '--md5',
        dest='motif_md5_to_motif_id_filename',
        action='store',
        type=str,
        required=False,
        help='Filename with motif MD5 to motif ID mappings to map Cluster-Buster motif MD5 filenames to motif IDs.'
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='db_prefix',
        action='store',
        type=str,
        required=True,
        help='Feather database prefix output filename.'
    )

    parser.add_argument(
        '-c',
        '--cbust',
        dest='cluster_buster_path',
        action='store',
        type=str,
        required=False,
        default='cbust',
        help='Path to Cluster-Buster (https://github.com/weng-lab/cluster-buster/). Default: "cbust".'
    )

    parser.add_argument(
        '-t',
        '--threads',
        dest='nbr_threads',
        action='store',
        type=int,
        required=False,
        default=1,
        help='Number of threads to use when scoring motifs. Default: 1.'
    )

    parser.add_argument(
        '-p',
        '--partial',
        dest='partial',
        nargs=2,
        metavar=('CURRENT_PART', 'NBR_TOTAL_PARTS'),
        action='store',
        type=int,
        required=False,
        help='Divide the motif list in a number of total parts (of similar size) and score only the part defined by '
             'current_part. This allows creating partial databases on machines which do not have enough RAM to score '
             'all motifs in one iteration. This will only create a partial regions/genes vs motifs scoring database '
             '({db_prefix}.part_000{current_part}_of_000{nbr_total_parts}.regions_vs_motifs.scores.feather or '
             '{db_prefix}.part_000{current_part}_of_000{nbr_total_parts}.genes_vs_motifs.scores.feather).'
    )

    parser.add_argument(
        '-g',
        '--genes',
        dest='extract_gene_id_from_region_id_regex_replace',
        action='store',
        type=str,
        required=False,
        default=None,
        help='Take top CRM score for a gene by taking the maximum CRM score of multiple regions for that gene. '
        'Define a regex which will remove the non-gene part of the region ID, so only the gene ID remains. '
        'Examples: "gene_id#some_number": "#[0-9]+$" or "region_id@@gene_id": "^.+@@".'
    )

    parser.add_argument(
        '-b',
        '--bgpadding',
        dest='bg_padding',
        action='store',
        type=int,
        required=False,
        default=0,
        help='Background padding in bp that was added for each sequence in FASTA file. Default: 0.'
    )

    parser.add_argument(
        '--min',
        dest='min_nbr_motifs',
        action='store',
        type=int,
        required=False,
        default=1,
        help='Minimum number of motifs needed per Cluster-Buster motif file to be considered for scoring '
             '(filters motifs list). Default: 1.'
    )

    parser.add_argument(
        '--max',
        dest='max_nbr_motifs',
        action='store',
        type=int,
        required=False,
        default=None,
        help='Maximum number of motifs needed per Cluster-Buster motif file to be considered for scoring '
             '(filters motifs list). Default: None.'
    )

    parser.add_argument(
        '-l',
        '--mask',
        dest='mask',
        action='store_true',
        help='Consider masked (lowercase) nucleotides as Ns.'
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
             'the same rankings databases as output.'
    )

    parser.add_argument(
        '-r',
        '--ssh',
        dest='ssh_command',
        action='store',
        type=str,
        required=False,
        help='If defined, run Cluster-Buster over ssh by running the provided command to make the connection before '
             'running Cluster-Buster itself. '
             "Example: 'ssh -o ControlMaster=auto -o ControlPath=/tmp/ssh-control-path-%%l-%%h-%%p-%%r -o ControlPersist=600 <hostname>'"
    )

    args = parser.parse_args()

    if not os.path.exists(args.fasta_filename):
        print(
            f'Error: FASTA filename "{args.fasta_filename}" does not exist.',
            file=sys.stderr
        )
        sys.exit(1)

    if args.original_species_fasta_filename and not os.path.exists(args.original_species_fasta_filename):
        print(
            f'Error: Original species FASTA filename "{args.original_species_fasta_filename}" does not exist.',
            file=sys.stderr
        )
        sys.exit(1)

    if not os.path.exists(args.motifs_dir):
        print(
            f'Error: Motif directory "{args.motifs_dir}" does not exist.',
            file=sys.stderr
        )
        sys.exit(1)

    if args.motif_md5_to_motif_id_filename and not os.path.exists(args.motif_md5_to_motif_id_filename):
        print(
            f'Error: Motif MD5 to motif ID mappings filename "{args.motif_md5_to_motif_id_filename}" does not exist.',
            file=sys.stderr
        )
        sys.exit(1)

    if not os.path.exists(args.motifs_list_filename):
        print(
            f'Error: Motifs list filename "{args.motifs_list_filename}" does not exist.',
            file=sys.stderr
        )
        sys.exit(1)

    if os.path.dirname(args.db_prefix) and not os.path.exists(os.path.dirname(args.db_prefix)):
        print(
            f'Error: Parent directory "{os.path.dirname(args.db_prefix)}" for Feather database prefix output filename '
            'does not exist.',
            file=sys.stderr
        )
        sys.exit(1)

    if args.partial:
        current_part, nbr_total_parts = args.partial

        if current_part < 1 or current_part > nbr_total_parts:
            print(
                f'Error: Current part ({current_part}) should be between 1 and the number of total parts '
                f'({nbr_total_parts}).',
                file=sys.stderr
            )
            sys.exit(1)

        # Add info about which part of the database this wil be.
        db_prefix = f'{args.db_prefix}.part_{current_part:04d}_of_{nbr_total_parts:04d}'
    else:
        db_prefix = args.db_prefix

    if not (args.min_nbr_motifs == 1 and not args.max_nbr_motifs):
        if not args.partial:
            # Add ".part_0001_of_0001" if partial was not specified but if min and max number of motifs per
            # Cluster-Buster motif file is specified, so it is easier to combine those databases with different subsets
            # of motifs in a later step.
            db_prefix = f'{args.db_prefix}.part_0001_of_0001'

        # Add info about number of motifs in Cluster-Buster motif file which will be scored (if min or max is set to a
        # non-default value).
        if args.max_nbr_motifs:
            db_prefix = f'{db_prefix}.min_{args.min_nbr_motifs:d}_to_max_{args.max_nbr_motifs:d}_motifs'
        else:
            db_prefix = f'{db_prefix}.min_{args.min_nbr_motifs:d}_to_max_motifs'

    # Get absolute path to Cluster-Buster binary and see if it can be executed.
    cluster_buster_path = shutil.which(args.cluster_buster_path)

    if not cluster_buster_path:
        print(f'Error: Cluster-Buster binary ("{args.cluster_buster_path}") could not be found or is not executable.')
        sys.exit(1)

    # Set random seed to provided input value or a random integer.
    seed = args.seed if args.seed else random.randint(0, 2**64)

    if args.original_species_fasta_filename:
        # When creating cisTarget databases for a species with lifted over regions, check if the regions/genes in the
        # current species FASTA file are available in the original species FASTA file. Due to liftover, regions might
        # be lost in the current species, but the cisTarget database needs to contain all regions/genes from the
        # original species to create the cisTarget cross-species database later.

        # Get all region or gene IDs from current species FASTA sequence names as a RegionOrGeneIDs object.
        region_or_gene_ids_current_species = RegionOrGeneIDs.get_region_or_gene_ids_from_fasta(
            fasta_filename=args.fasta_filename,
            extract_gene_id_from_region_id_regex_replace=args.extract_gene_id_from_region_id_regex_replace
        )

        # Get all region or gene IDs from the original FASTA sequence names as a RegionOrGeneIDs object.
        region_or_gene_ids = RegionOrGeneIDs.get_region_or_gene_ids_from_fasta(
            fasta_filename=args.original_species_fasta_filename,
            extract_gene_id_from_region_id_regex_replace=args.extract_gene_id_from_region_id_regex_replace
        )

        if not region_or_gene_ids_current_species.issubset(region_or_gene_ids):
            print(
                f'Error: Region IDs/gene IDs in "{args.fasta_filename}" are not all present in '
                f'"{args.original_species_fasta_filename}".',
                file=sys.stderr
            )
            sys.exit(1)
    else:
        # Get all region or gene IDs from the FASTA sequence names as a RegionOrGeneIDs object.
        region_or_gene_ids = RegionOrGeneIDs.get_region_or_gene_ids_from_fasta(
            fasta_filename=args.fasta_filename,
            extract_gene_id_from_region_id_regex_replace=args.extract_gene_id_from_region_id_regex_replace
        )

    # Get absolute path name for FASTA filename so in case Cluster-Buster is ran over ssh, the FASTA file can be found.
    fasta_filename = os.path.abspath(args.fasta_filename)

    # Get motif ID to motif file name mapping and motif ID to number of motifs per motif file mapping for
    # and (optionally) filtered list of motif IDs:
    #   - if partial is set
    #   - if min_nbr_motifs is set
    #   - if max_nbr_motifs is set
    motif_id_to_filename_dict, motif_id_to_nbr_motifs_dict = get_motif_id_to_filename_and_nbr_motifs_dict(
        motifs_dir=os.path.abspath(args.motifs_dir),
        motifs_list_filename=args.motifs_list_filename,
        motif_md5_to_motif_id_filename=args.motif_md5_to_motif_id_filename,
        partial=(current_part, nbr_total_parts) if args.partial else None,
        min_nbr_motifs=args.min_nbr_motifs,
        max_nbr_motifs=args.max_nbr_motifs
    )

    # Create MotifOrTracksIDs object from plain motif IDs.
    motif_ids = MotifOrTrackIDs(
        motif_or_track_ids=set(motif_id_to_filename_dict),
        motifs_or_tracks_type=MotifsOrTracksType.MOTIFS
    )

    nbr_region_or_gene_ids = len(region_or_gene_ids)
    nbr_motifs = len(motif_id_to_filename_dict)

    if nbr_region_or_gene_ids == 0:
        print(f'Error: No {region_or_gene_ids.type.value} provided.', file=sys.stderr)
        sys.exit(1)

    if nbr_motifs == 0:
        print('Error: No motifs provided.', file=sys.stderr)
        sys.exit(1)

    print(
        f'Initialize dataframe ({nbr_region_or_gene_ids} {region_or_gene_ids.type.value} '
        f'x {nbr_motifs} motifs) for storing CRM scores for each {region_or_gene_ids.type.value} per motif.',
        file=sys.stderr
    )

    ct_scores_db_motifs_vs_regions_or_genes = CisTargetDatabase.create_db(
        db_type=DatabaseTypes.from_strings(
            scores_or_rankings='scores',
            column_kind='motifs',
            row_kind=region_or_gene_ids.type.value
        ),
        region_or_gene_ids=region_or_gene_ids,
        motif_or_track_ids=motif_ids,
        order='F'
    )

    def write_crm_scores_for_motif_to_ct_scores_db(df_motif_id_and_crm_scores: Tuple[str, pd.DataFrame]) -> None:
        if 'nbr_of_scored_motifs' not in write_crm_scores_for_motif_to_ct_scores_db.__dict__:
            write_crm_scores_for_motif_to_ct_scores_db.nbr_of_scored_motifs = 0

        motif_id, df_crm_scores = df_motif_id_and_crm_scores

        start_time = time.monotonic()
        ct_scores_db_motifs_vs_regions_or_genes.update_scores_for_motif_or_track(
            motif_or_track_id=motif_id,
            df_scores_for_motif_or_track=df_crm_scores['crm_score']
        )
        elapsed_time = time.monotonic() - start_time

        write_crm_scores_for_motif_to_ct_scores_db.nbr_of_scored_motifs += 1

        print(
            f'Adding Cluster-Buster CRM scores ({write_crm_scores_for_motif_to_ct_scores_db.nbr_of_scored_motifs:d} of '
            f'{nbr_motifs:d}) for motif "{motif_id:s}" took {elapsed_time:0.6f} seconds.',
            file=sys.stderr
        )

    def report_error(exception: BaseException) -> None:
        print(exception, file=sys.stderr)

    start_time = time.monotonic()

    with mp.Pool(processes=args.nbr_threads) as pool:
        # Motif IDs are sorted by number of motifs in motif ID Cluster-Buster file (high to low) and then by motif ID
        # name (in get_motif_id_to_filename_and_nbr_motifs_dict()), so motif IDs which have a lot of motifs in their
        # Cluster-Buster motif file are scored first, before singletons are scored to prevent that a few Cluster-Buster
        # motifs files with a huge number of motifs got sheduled near the end of the motif scoring and underutilizing
        # the compute node.
        for motif_id, motif_filename in motif_id_to_filename_dict.items():

            # Score all regions/genes in the FASTA file for the current motif and write the result in the
            # ct_scores_db_motifs_vs_regions_or_genes CisTargetDatabase object.
            pool.apply_async(
                func=run_cluster_buster_for_motif,
                args=[
                    cluster_buster_path,
                    fasta_filename,
                    motif_filename,
                    motif_id,
                    args.extract_gene_id_from_region_id_regex_replace,
                    args.bg_padding,
                    args.mask,
                    args.ssh_command
                ],
                callback=write_crm_scores_for_motif_to_ct_scores_db,
                error_callback=report_error
            )

        # Prevents any more task from being submitted to the pool.
        pool.close()

        # Wait for worker processes to exit.
        pool.join()

    elapsed_time = time.monotonic() - start_time

    print(
        f'\nScoring {nbr_motifs} motifs with Cluster-Buster took: {elapsed_time:.06f} seconds\n'
    )

    if 'nbr_of_scored_motifs' not in write_crm_scores_for_motif_to_ct_scores_db.__dict__:
        print(
            f'Error: None of {nbr_motifs:d} motifs were scored successfully.',
            file=sys.stderr
        )
        sys.exit(1)
    elif write_crm_scores_for_motif_to_ct_scores_db.nbr_of_scored_motifs != nbr_motifs:
        print(
            f'Error: Only {write_crm_scores_for_motif_to_ct_scores_db.nbr_of_scored_motifs:d} out of {nbr_motifs:d} '
            f'motifs were scored successfully.',
            file=sys.stderr
        )
        sys.exit(1)

    print('', file=sys.stderr)

    def write_db(ct_db: CisTargetDatabase, db_prefix: str):
        """
        Write cisTarget database to a Feather file and print database location and elapsed time.

        :param ct_db: cisTarget database object.
        :param db_prefix: Feather database file prefix.
        :return:
        """
        db_filename = ct_db.create_db_filename_from_db_prefix(db_prefix=db_prefix, extension='feather')

        print(
            f'Writing cisTarget {ct_db.db_type.row_kind} vs {ct_db.db_type.column_kind} '
            f'{ct_db.db_type.scores_or_rankings} db: "{db_filename}"'
        )

        start_time = time.monotonic()
        ct_db.write_db(
            db_prefix=db_prefix,
            version=1
        )
        elapsed_time = time.monotonic() - start_time

        print(
            f'Writing cisTarget {ct_db.db_type.row_kind} vs {ct_db.db_type.column_kind} '
            f'{ct_db.db_type.scores_or_rankings} db took: {elapsed_time:.06f} seconds\n'
        )

    # Write cisTarget scores database (motifs vs regions or genes) to Feather file.
    write_db(ct_db=ct_scores_db_motifs_vs_regions_or_genes, db_prefix=db_prefix)

    if not args.partial:
        # Create cisTarget scores database (regions or genes vs motifs) from (motifs vs regions or genes) version.
        ct_scores_db_regions_or_genes_vs_motifs = ct_scores_db_motifs_vs_regions_or_genes.transpose()

        # Write cisTarget scores database (regions or genes vs motifs) to Feather file.
        write_db(ct_db=ct_scores_db_regions_or_genes_vs_motifs, db_prefix=db_prefix)

        # Create cisTarget rankings database (motifs vs regions or genes) from cisTarget scores database filename
        # (motifs vs regions or genes).
        print(
            f'''Create rankings from "{
                ct_scores_db_motifs_vs_regions_or_genes.create_db_filename_from_db_prefix(
                    db_prefix=db_prefix,
                    extension='feather'
                )
            }" with random seed set to {seed}.''',
            file=sys.stderr
        )

        start_time = time.monotonic()
        ct_rankings_db_motifs_vs_regions_or_genes = \
            ct_scores_db_motifs_vs_regions_or_genes.convert_scores_db_to_rankings_db(seed=seed)
        elapsed_time = time.monotonic() - start_time

        print(
            f'Creating cisTarget rankings db from cisTarget scores db took: '
            f'{elapsed_time:.06f} seconds\n'
        )

        # Reclaim memory occupied by cisTarget scores databases.
        del ct_scores_db_motifs_vs_regions_or_genes
        del ct_scores_db_regions_or_genes_vs_motifs

        # Do not write cisTarget rankings database (motifs vs regions or genes) to Feather file
        # as it can take a very long time to write it (1.5 hours for 1 milion regions) as the
        # rankings database numpy array is in "C" order and writing a Feather database requires
        # travering the numpy array in column order.
        #write_db(ct_db=ct_rankings_db_motifs_vs_regions_or_genes, db_prefix=db_prefix)

        # Create cisTarget rankings database (regions or genes vs motifs) from (motifs vs regions or genes) version.
        ct_rankings_db_regions_or_genes_vs_motifs = ct_rankings_db_motifs_vs_regions_or_genes.transpose()

        # Write cisTarget rankings database (regions or genes vs motifs) to Feather file.
        write_db(ct_db=ct_rankings_db_regions_or_genes_vs_motifs, db_prefix=db_prefix)


if __name__ == '__main__':
    main()
