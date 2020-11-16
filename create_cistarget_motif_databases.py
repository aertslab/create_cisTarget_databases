#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Create cisTarget motif databases.

Copyright (C): 2019-2020 - Gert Hulselmans
"""


import argparse
import io
import os
import random
import re
import shutil
import subprocess
import sys
import time
import multiprocessing as mp

from typing import Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd

from cistarget_db import FeaturesType, MotifsOrTracksType, FeatureIDs, MotifOrTrackIDs, DatabaseTypes, CisTargetDatabase


def get_motif_id_to_filename_dict(motifs_dir: str,
                                  motifs_list_filename: str,
                                  motif_md5_to_motif_id_filename: Optional[str] = None
                                  ) -> Dict[str, str]:
    """
    Create motif ID to Cluster-Buster motif file mapping.

    :param motifs_dir: Directory with Cluster-Buster motif files (with motif MD5 name or motif ID motif files).
    :param motifs_list_filename: File with Cluster-Buster motif MD5 names or motif IDs.
    :param motif_md5_to_motif_id_filename: TSV file with motif MD5 names to motif IDs mapping (optional).
    :return: motif_id_to_filename_dict: motif ID to CLuster-Buster motif filename mapping.
    """

    motif_id_to_filename_dict = dict()
    motif_md5_to_motif_id_dict = dict()
    motif_id_to_motif_md5_dict = dict()

    if motif_md5_to_motif_id_filename:
        # Get motif MD5 name to motif ID mapping if motif_md5_to_motif_id_filename was provided.
        with open(motif_md5_to_motif_id_filename, 'r') as fh:
            for line in fh:
                line = line.rstrip()

                if line and not line.startswith('#'):
                    motif_md5, motif_id = line.rstrip().split('\t')[0:2]

                    # Store motif MD5 name to motif ID mapping and vice versa.
                    motif_md5_to_motif_id_dict[motif_md5] = motif_id
                    motif_id_to_motif_md5_dict[motif_id] = motif_md5

    # Create motif ID to Cluster-Buster motif filename mapping.
    with open(motifs_list_filename, 'r') as fh:
        for line in fh:
            motif_md5_or_id = line.rstrip()

            if motif_md5_or_id and not motif_md5_or_id.startswith('#'):
                if motif_md5_or_id.endswith('.cb'):
                    # Remove ".cb" extension from motif MD5 name or motif ID.
                    motif_md5_or_id = motif_md5_or_id[:-3]

                if motif_md5_to_motif_id_dict:
                    # A motif_md5_to_motif_id_filename was provided, so assume Cluster-Buster motif filenames in
                    # motifs_dir have motif MD5 names.
                    if motif_md5_or_id in motif_md5_to_motif_id_dict:
                        # Get associated motif ID for motif MD5 name if a motif MD5 name was provided.
                        motif_id = motif_md5_to_motif_id_dict[motif_md5_or_id]
                        motif_md5 = motif_md5_or_id
                    elif motif_md5_or_id in motif_id_to_motif_md5_dict:
                        # Get associated motif MD5 name for motif ID if a motif ID was provided.
                        motif_id = motif_md5_or_id
                        motif_md5 = motif_id_to_motif_md5_dict[motif_md5_or_id]
                    else:
                        print(
                            f'Error: Could not find motif MD5 name <=> motif ID association for "{motif_md5_or_id}".',
                            file=sys.stderr
                        )
                        sys.exit(1)

                    # Cluster-Buster motif MD5 name filename.
                    motif_filename = os.path.join(motifs_dir, motif_md5 + '.cb')
                else:
                    # No motif_md5_to_motif_id_filename was provided, so assume Cluster-Buster motif filenames in
                    # motifs_dir have motif IDs.

                    motif_id = motif_md5_or_id
                    # Cluster-Buster motif ID filename.
                    motif_filename = os.path.join(motifs_dir, motif_id + '.cb')

                if not os.path.exists(motif_filename):
                    print(
                        f'Error: Cluster-Buster motif filename "{motif_filename}" does not exist for motif {motif_id}.',
                        file=sys.stderr
                    )
                    sys.exit(1)

                motif_id_to_filename_dict[motif_id] = motif_filename

    return motif_id_to_filename_dict


def get_region_ids_or_gene_ids_from_fasta(fasta_filename: str,
                                          extract_gene_id_from_region_id_regex_replace: Optional[str] = None
                                          ) -> FeatureIDs:
    """
    Get all region IDs or gene IDs from FASTA filename:
      - When extract_gene_id_from_region_id_regex_replace=None, region IDs are returned and each region ID is only
        allowed once in the FASTA file.
      - When extract_gene_id_from_region_id_regex_replace is set to a regex to remove the non gene ID part from the
        region IDs, gene IDs are returned and each gene is allowed to appear more than once in the FASTA file.

    :param fasta_filename:
         FASTA filename with sequences for region IDs or gene IDs.
    :param extract_gene_id_from_region_id_regex_replace:
         regex for removing unwanted parts from the region ID to extract the gene ID.
    :return: FeatureIDs object for regions or genes.
    """

    gene_ids = set()
    region_ids = set()
    duplicated_region_ids = False

    with open(fasta_filename, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                # Get region ID by getting everything after '>' up till the first whitespace.
                region_id = line[1:].split(maxsplit=1)[0]

                if extract_gene_id_from_region_id_regex_replace:
                    # Extract gene ID from region ID.
                    gene_id = re.sub(
                        extract_gene_id_from_region_id_regex_replace,
                        '',
                        region_id
                    )

                    gene_ids.add(gene_id)
                else:
                    # Check if all region IDs only appear once.
                    if region_id in region_ids:
                        print(
                            f'Error: region ID "{region_id:s}" is not unique in FASTA file "{fasta_filename:s}".',
                            file=sys.stderr
                        )

                        duplicated_region_ids = True

                    region_ids.add(region_id)

    if extract_gene_id_from_region_id_regex_replace:
        return FeatureIDs(feature_ids=gene_ids, features_type=FeaturesType.GENES)
    else:
        if duplicated_region_ids:
            sys.exit(1)

        return FeatureIDs(feature_ids=region_ids, features_type=FeaturesType.REGIONS)


def run_cluster_buster_for_motif(cluster_buster_path: str, fasta_filename: str, motif_filename: str, motif_id: str,
                                 extract_gene_id_from_region_id_regex_replace: Optional[str] = None,
                                 bg_padding: int = 0, mask: bool = False,
                                 ssh_command: Optional[Union[str, list]] = None
                                 ) -> Tuple[str, pd.DataFrame]:
    """
    Score each sequence in the FASTA file with Cluster-Buster and only keep the top CRM score per region ID/gene ID.

    :param cluster_buster_path: Path to Cluster-Buster binary.
    :param fasta_filename:      FASTA filename with regions to score.
    :param motif_filename:      Cluster-Buster motif filename which contains the motif to score all regions with.
    :param motif_id:            Motif ID.
    :param extract_gene_id_from_region_id_regex_replace:
                                Define a regex which will remove the non-gene part of the region ID of each sequence
                                name in the FASTA file, so only the gene ID remains. If set to None the whole region ID
                                will be kept instead. In case of region IDs, the best CRM score per region is kept.
                                In case of gene IDs, the best CRM score from multiple regions is kept.
    :param bg_padding:          Use X bp at start and end of each sequence only for calculating the background
                                nucleotide frequency, but not for scoring the motif itself.
    :param mask:                Consider masked (lowercase) nucleotides as Ns.
    :param ssh_command:         If defined, run Cluster-Buster over ssh by running the provided command to make the
                                connection before running Cluster-Buster.
                                Example : 'ssh -o ControlMaster=auto -o ControlPath=/tmp/ssh-control-path-%l-%h-%p-%r -o ControlPersist=600 <hostname>'
    :return:                    (motif_id, df_crm_scores): motif ID and dataframe with top CRM score per region/gene ID.
    """

    clusterbuster_command = []

    if ssh_command:
        # Add SSH command to the start of the Cluster-Buster command.
        if isinstance(ssh_command, str):
            clusterbuster_command.extend(ssh_command.split())
        elif isinstance(ssh_command, list):
            clusterbuster_command.extend(ssh_command)

    # Construct Cluster-Buster command line.
    clusterbuster_command.extend([
        cluster_buster_path,
        '-f', '4',
        '-c', '0.0',
        '-r', '10000',
        '-b', str(bg_padding),
        '-t', '1'
    ])

    if mask:
        clusterbuster_command.append('-l')

    clusterbuster_command.extend([
        motif_filename,
        fasta_filename
    ])

    # Score each region in FASTA file with Cluster-Buster for the provided motif and get top CRM score for each region.
    try:
        pid = subprocess.Popen(args=clusterbuster_command,
                               bufsize=-1,
                               executable=None,
                               stdin=None,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               preexec_fn=None,
                               close_fds=False,
                               shell=False,
                               cwd=None,
                               env=None,
                               universal_newlines=False,
                               startupinfo=None,
                               creationflags=0)
        stdout_data, stderr_data = pid.communicate()
    except OSError as msg:
        raise RuntimeError("Execution error for: '" + ' '.join(clusterbuster_command) + "': " + str(msg))

    if pid.returncode != 0:
        raise RuntimeError("Error: Non-zero exit status for: '" + ' '.join(clusterbuster_command) + "'")

    # Read Cluster-Buster standard out as a pandas dataframe.
    df_crm_scores = pd.read_csv(
        filepath_or_buffer=io.BytesIO(stdout_data),
        sep='\t',
        header=0,
        names=['seq_name', 'crm_score', 'seq_number', 'rank'],
        index_col='seq_name',
        usecols=['seq_name', 'crm_score'],
        dtype={'seq_name': str, 'crm_score': np.float32},
        engine='c'
    )

    if extract_gene_id_from_region_id_regex_replace:
        # Extract gene ID from the region ID by removing the non-gene part.
        #
        # Take the top CRM score for each gene ID by taking the maximum CRM
        # score of all region IDs that belong to the same gene ID.
        #
        # Examples:
        #   - extract_gene_id_from_region_id_regex_replace = '#[0-9]+$'
        #       - region ID (input):
        #           "geneA#1", "geneB#1", "geneC#1", "geneA#2", "geneC#2"
        #       - gene IDs (output):
        #           "geneA", "geneB", "geneC"
        #   - extract_gene_id_from_region_id_regex_replace = '^.+@@'
        #       - region IDs (input):
        #           "region1@@geneA", "region2@@geneB", "region3@@geneA"
        #       - gene IDs (output):
        #           "geneA", "geneB"
        df_crm_scores = df_crm_scores.assign(
            gene_ids=df_crm_scores.index.str.replace(
                extract_gene_id_from_region_id_regex_replace,
                '',
                regex=True
            )
        ).groupby('gene_ids').max()

    return motif_id, df_crm_scores


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
             'creating a species CisTarget database from regions/genes lifted over from a different species, provide '
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
             'create a cross species CisTarget database later, all individual species CisTarget databases need to '
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
        help='Filename with list of motif IDs or motif MD5 names to be scored from directory specified by "--motifs_dir".'
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

    # Get absolute path to Cluster-Buster binary and see if it can be executed.
    cluster_buster_path = shutil.which(args.cluster_buster_path)

    if not cluster_buster_path:
        print(f'Error: Cluster-Buster binary ("{args.cluster_buster_path}") could not be found or is not executable.')
        sys.exit(1)

    # Set random seed to provided input value or a random integer.
    seed = args.seed if args.seed else random.randint(0, 2**64)

    if args.original_species_fasta_filename:
        # When creating CisTarget databases for a species with lifted over regions, check if the regions/genes in the
        # current species FASTA file are available in the original species FASTA file. Due to liftover, regions might
        # be lost in the current species, but the CisTarget database needs to contain all regions/genes from the
        # original species to create the cross species CisTarget database later.

        # Get all region IDs or gene IDs from current species FASTA sequence names as a FeaturesIDs object.
        region_ids_or_gene_ids_current_species = get_region_ids_or_gene_ids_from_fasta(
            fasta_filename=args.fasta_filename,
            extract_gene_id_from_region_id_regex_replace=args.extract_gene_id_from_region_id_regex_replace
        )

        # Get all region IDs or gene IDs from the original FASTA sequence names as a FeaturesIDs object.
        region_ids_or_gene_ids = get_region_ids_or_gene_ids_from_fasta(
            fasta_filename=args.original_species_fasta_filename,
            extract_gene_id_from_region_id_regex_replace=args.extract_gene_id_from_region_id_regex_replace
        )

        if not region_ids_or_gene_ids_current_species.issubset(region_ids_or_gene_ids):
            print(
                f'Error: Region IDs/gene IDs in "{args.fasta_filename}" are not all present in '
                f'"{args.original_species_fasta_filename}".',
                file=sys.stderr
            )
            sys.exit(1)
    else:
        # Get all region IDs or gene IDs from the FASTA sequence names as a FeaturesIDs object.
        region_ids_or_gene_ids = get_region_ids_or_gene_ids_from_fasta(
            fasta_filename=args.fasta_filename,
            extract_gene_id_from_region_id_regex_replace=args.extract_gene_id_from_region_id_regex_replace
        )

    # Get absolute path name for FASTA filename so in case Cluster-Buster is ran over ssh, the FASTA file can be found.
    fasta_filename = os.path.abspath(args.fasta_filename)

    motif_id_to_filename_dict = get_motif_id_to_filename_dict(
        motifs_dir=os.path.abspath(args.motifs_dir),
        motifs_list_filename=args.motifs_list_filename,
        motif_md5_to_motif_id_filename=args.motif_md5_to_motif_id_filename
    )

    # Create MotifOrTracksIDs object from plain motif IDs.
    motif_ids = MotifOrTrackIDs(
        motif_or_track_ids=set(motif_id_to_filename_dict),
        motifs_or_tracks_type=MotifsOrTracksType.MOTIFS
    )

    nbr_region_ids_or_gene_ids = len(region_ids_or_gene_ids)
    nbr_motifs = len(motif_id_to_filename_dict)

    if nbr_region_ids_or_gene_ids == 0:
        print(f'Error: No {region_ids_or_gene_ids.type.value} provided.', file=sys.stderr)
        sys.exit(1)

    if nbr_motifs == 0:
        print('Error: No motifs provided.', file=sys.stderr)
        sys.exit(1)

    if args.partial:
        def divide_in_chunks(list_to_chunk, n):
            """
            Yield n number of sequential chunks from list_to_chunk, where each chunk is similar in size (max diff = 1).
            """

            d, r = divmod(len(list_to_chunk), n)

            for i in range(n):
                si = (d + 1) * (i if i < r else r) + d * (0 if i < r else i - r)
                yield list_to_chunk[si:si + (d + 1 if i < r else d)]

        # Create a new MotifOrTracksIDs object with only the motif IDs for the requested part (defined by current_part).
        motif_ids = MotifOrTrackIDs(
            motif_or_track_ids=list(divide_in_chunks(list_to_chunk=motif_ids.ids, n=nbr_total_parts))[current_part - 1],
            motifs_or_tracks_type=MotifsOrTracksType.MOTIFS
        )

        motif_id_to_filename_dict = {
            motif_id: filename
            for motif_id, filename in motif_id_to_filename_dict.items()
            if motif_id in motif_ids.ids
        }

        nbr_motifs = len(motif_ids)

    print(
        f'Initialize dataframe ({nbr_region_ids_or_gene_ids} {region_ids_or_gene_ids.type.value} '
        f'x {nbr_motifs} motifs) for storing CRM scores for each {region_ids_or_gene_ids.type.value} per motif.',
        file=sys.stderr
    )

    ct_scores_db_motifs_vs_regions_or_genes = CisTargetDatabase.create_db(
        db_type=DatabaseTypes.from_strings(
            scores_or_rankings='scores',
            column_kind='motifs',
            row_kind=region_ids_or_gene_ids.type.value
        ),
        feature_ids=region_ids_or_gene_ids,
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

    with mp.Pool(processes=args.nbr_threads) as pool:
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
            f'Write {ct_db.db_type.scores_or_rankings} of {ct_db.feature_ids.type.value} for each motif in '
            f'{ct_db.db_type.column_kind} vs {ct_db.db_type.row_kind} style to '
            f'"{db_filename}".',
            file=sys.stderr
        )

        start_time = time.monotonic()
        ct_db.write_db(
            db_prefix=db_prefix,
            version=1
        )
        elapsed_time = time.monotonic() - start_time

        print(f'Database written in {elapsed_time:0.6f} seconds.\n', file=sys.stderr)

    if not args.partial:
        # Write cisTarget scores database (motifs vs regions or genes) to Feather file.
        write_db(ct_db=ct_scores_db_motifs_vs_regions_or_genes, db_prefix=db_prefix)

    # Create cisTarget scores database (regions or genes vs motifs) from (motifs vs regions or genes) version.
    ct_scores_db_regions_or_genes_vs_motifs = ct_scores_db_motifs_vs_regions_or_genes.transpose()

    # Write cisTarget scores database (regions or genes vs motifs) to Feather file.
    write_db(ct_db=ct_scores_db_regions_or_genes_vs_motifs, db_prefix=db_prefix)

    if not args.partial:
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

        print(f'Creating rankings from scores database took {elapsed_time:.06f} seconds.\n', file=sys.stderr)

        # Reclaim memory occupied by cisTarget scores databases.
        del ct_scores_db_motifs_vs_regions_or_genes
        del ct_scores_db_regions_or_genes_vs_motifs

        # Write cisTarget rankings database (motifs vs regions or genes) to Feather file.
        write_db(ct_db=ct_rankings_db_motifs_vs_regions_or_genes, db_prefix=db_prefix)

        # Create cisTarget rankings database (regions or genes vs motifs) from (motifs vs regions or genes) version.
        ct_rankings_db_regions_or_genes_vs_motifs = ct_rankings_db_motifs_vs_regions_or_genes.transpose()

        # Write cisTarget rankings database (regions or genes vs motifs) to Feather file.
        write_db(ct_db=ct_rankings_db_regions_or_genes_vs_motifs, db_prefix=db_prefix)


if __name__ == '__main__':
    main()
