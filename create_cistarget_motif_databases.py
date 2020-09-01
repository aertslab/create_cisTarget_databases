#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Create cisTarget motif databases.

Copyright (C): 2019-2020 - Gert Hulselmans
"""


import argparse
import io
import os
import re
import subprocess
import sys
import multiprocessing as mp

import numpy as np
import pandas as pd


def get_motif_id_to_filename_dict(motifs_dir, motifs_list_filename, motif_md5_to_motif_id_filename=None):
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


def get_region_ids_or_gene_ids_from_fasta(fasta_filename, extract_gene_id_from_region_id_regex_replace=None):
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
    :return: ('regions', region_ids) or ('genes', gene_ids).
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
                            'Error: region ID "{0:s}" is not unique in FASTA file "{1:s}".'.format(
                                region_id,
                                fasta_filename
                            ),
                            file=sys.stderr
                        )

                        duplicated_region_ids = True

                    region_ids.add(region_id)

    if extract_gene_id_from_region_id_regex_replace:
        # Return ('genes', gene IDs).
        return 'genes', sorted(gene_ids)
    else:
        if duplicated_region_ids:
            sys.exit(1)

        # Return ('regions', region IDs).
        return 'regions', sorted(region_ids)


def run_cluster_buster_for_motif(cluster_buster_path, fasta_filename, motif_filename, motif_id,
                                 extract_gene_id_from_region_id_regex_replace=None, bg_padding=0):
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
    :return:                    (motif_id, df_crm_scores): motif ID and dataframe with top CRM score per region/gene ID.
    """

    # Score each region in FASTA file with Cluster-Buster for the provided motif and get top CRM score for each region.
    clusterbuster_command = [cluster_buster_path,
                             '-f', '4',
                             '-c', '0.0',
                             '-r', '10000',
                             '-b', str(bg_padding),
                             '-t', '1',
                             motif_filename,
                             fasta_filename]

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
        help='FASTA filename which contains the regions to score with Cluster-Buster for each motif.'
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
        dest='feature_table_output_filename',
        action='store',
        type=str,
        required=True,
        help='Feature table output filename.'
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

    args = parser.parse_args()

    if not os.path.exists(args.fasta_filename):
        print(
            'Error: FASTA filename "{0:s}" does not exist.'.format(args.fasta_filename),
            file=sys.stderr
        )
        sys.exit(1)

    if not os.path.exists(args.motifs_dir):
        print(
            'Error: Motif directory "{0:s}" does not exist.'.format(args.motifs_dir),
            file=sys.stderr
        )
        sys.exit(1)

    if args.motif_md5_to_motif_id_filename and not os.path.exists(args.motif_md5_to_motif_id_filename):
        print(
            'Error: Motif MD5 to motif ID mappings filename "{0:s}" does not exist.'.format(args.motif_md5_to_motif_id_filename),
            file=sys.stderr
        )
        sys.exit(1)

    if not os.path.exists(args.fasta_filename):
        print(
            'Error: Motifs list filename "{0:s}" does not exist.'.format(args.motifs_list_filename),
            file=sys.stderr
        )
        sys.exit(1)

    motif_id_to_filename_dict = get_motif_id_to_filename_dict(
        motifs_dir=args.motifs_dir,
        motifs_list_filename=args.motifs_list_filename,
        motif_md5_to_motif_id_filename=args.motif_md5_to_motif_id_filename
    )

    # Get type of sequences ("regions" or "genes")
    # and sorted list of region IDs or gene IDs.
    (regions_or_genes_type,
     region_ids_or_gene_ids) = get_region_ids_or_gene_ids_from_fasta(
        args.fasta_filename,
        args.extract_gene_id_from_region_id_regex_replace
    )

    nbr_region_ids_or_gene_ids = len(region_ids_or_gene_ids)
    nbr_motifs = len(motif_id_to_filename_dict)

    if nbr_region_ids_or_gene_ids == 0:
        print(f'Error: No {regions_or_genes_type} provided.', file=sys.stderr)
        sys.exit(1)

    if nbr_motifs == 0:
        print('Error: No motifs provided.', file=sys.stderr)
        sys.exit(1)

    print(
        f'Initialize dataframe ({nbr_region_ids_or_gene_ids} {regions_or_genes_type} ' \
        f'x {nbr_motifs} motifs) for storing CRM scores for each {regions_or_genes_type} ' \
        'per motif.',
        file=sys.stderr
    )

    # Create zeroed dataframe for all region IDs or gene IDs vs all motif IDs.
    df_scores__motifs_vs_regions_or_genes = pd.DataFrame(
        data=np.zeros((nbr_region_ids_or_gene_ids, nbr_motifs),
                      dtype=np.float32),
        index=region_ids_or_gene_ids,
        columns=sorted(motif_id_to_filename_dict.keys())
    )

    # Add index name: 'regions' or 'genes'.
    df_scores__motifs_vs_regions_or_genes.rename_axis(
        index=regions_or_genes_type,
        axis='index',
        copy=False,
        inplace=True
    )

    # Add column axis name: 'motifs'.
    df_scores__motifs_vs_regions_or_genes.rename_axis(
        columns='motifs',
        axis='columns',
        copy=False,
        inplace=True
    )

    def write_crm_scores_for_motif_to_df_scores(df_motif_id_and_crm_scores):
        if 'nbr_of_scored_motifs' not in write_crm_scores_for_motif_to_df_scores.__dict__:
            write_crm_scores_for_motif_to_df_scores.nbr_of_scored_motifs = 0

        motif_id, df_crm_scores = df_motif_id_and_crm_scores

        df_scores__motifs_vs_regions_or_genes.loc[df_crm_scores.index, motif_id] = df_crm_scores['crm_score']

        write_crm_scores_for_motif_to_df_scores.nbr_of_scored_motifs += 1

        print(
            'Adding Cluster-Buster CRM scores ({0:d} of {1:d}) for motif "{2:s}".'.format(
                write_crm_scores_for_motif_to_df_scores.nbr_of_scored_motifs,
                nbr_motifs,
                motif_id
            ),
            file=sys.stderr
        )

    def report_error(exception):
        print(exception, file=sys.stderr)

    with mp.Pool(processes=args.nbr_threads) as pool:
        for motif_id, motif_filename in motif_id_to_filename_dict.items():
            pool.apply_async(
                func=run_cluster_buster_for_motif,
                args=[
                    args.cluster_buster_path,
                    args.fasta_filename,
                    motif_filename,
                    motif_id,
                    args.extract_gene_id_from_region_id_regex_replace,
                    args.bg_padding
                ],
                callback=write_crm_scores_for_motif_to_df_scores,
                error_callback=report_error
            )

        # Prevents any more task from being submitted to the pool.
        pool.close()

        # Wait for worker processes to exit.
        pool.join()

    if 'nbr_of_scored_motifs' not in write_crm_scores_for_motif_to_df_scores.__dict__:
        print(
            'Error: None of {0:d} motifs was scored successfully'.format(nbr_motifs),
            file=sys.stderr
        )
        sys.exit(1)
    elif write_crm_scores_for_motif_to_df_scores.nbr_of_scored_motifs != nbr_motifs:
        print(
            'Error: Only {0:d} of {1:d} motifs was scored successfully'.format(
                write_crm_scores_for_motif_to_df_scores.nbr_of_scored_motifs,
                nbr_motifs
            ),
            file=sys.stderr
        )
        sys.exit(1)

    print(
        'Write cisTarget motif databases in feather format: ' \
        f'"{args.feature_table_output_filename}.*.*.feather".',
        file=sys.stderr
    )

    print(
        f'Write CRM scores of {regions_or_genes_type} for each motif in ' \
        f'motifs vs {regions_or_genes_type} style to ' \
        f'{args.feature_table_output_filename}.motifs_vs_{regions_or_genes_type}.scores.feather',
        file=sys.stderr
    )

    df_scores__motifs_vs_regions_or_genes.reset_index(inplace=True)
    df_scores__motifs_vs_regions_or_genes.to_feather(
        path=f'{args.feature_table_output_filename}.motifs_vs_{regions_or_genes_type}.scores.feather'
    )
    df_scores__motifs_vs_regions_or_genes.set_index(regions_or_genes_type, inplace=True)
    # Add column axis name: 'motifs'.
    df_scores__motifs_vs_regions_or_genes.rename_axis(
        columns='motifs',
        axis='columns',
        copy=False,
        inplace=True
    )

    print(
        f'Write CRM scores of {regions_or_genes_type} for each motif in ' \
        f'{regions_or_genes_type} vs motifs style to ' \
        f'{args.feature_table_output_filename}.{regions_or_genes_type}_vs_motifs.scores.feather',
        file=sys.stderr
    )

    df_scores__regions_or_genes_vs_genes = df_scores__motifs_vs_regions_or_genes.transpose()
    df_scores__regions_or_genes_vs_genes.reset_index(inplace=True)
    df_scores__regions_or_genes_vs_genes.to_feather(
        path=f'{args.feature_table_output_filename}.{regions_or_genes_type}_vs_motifs.scores.feather'
    )
    del df_scores__regions_or_genes_vs_genes

    def rank_CRM_scores_and_assign_random_ranking_in_range_for_ties_func(crm_scores_with_ties_for_motif_numpy):
        rng = np.random.default_rng()
        # Create random permutation so tied scores will have a different ranking each time.
        random_permutations_to_break_ties_numpy = rng.permutation(crm_scores_with_ties_for_motif_numpy.shape[0])

        # Rank CRM scores for each region/gene for the current motif and break ties:
        #   - Get current column with CRM scores for a certain motif and
        #     multiply by -1 (CRM scores >= 0) so sorting later will result
        #     in ranking the highest CRM score first (descending):
        #
        #       (-crm_scores_with_ties_for_motif_numpy)
        #
        #   - Access the negated CRM scores in a random order:
        #
        #       (-crm_scores_with_ties_for_motif_numpy)[random_permutations_to_break_ties_numpy]
        #
        #     so when sorting it, CRM scores for regions/genes at the start
        #     of the array do not get artificially better rankings than
        #     regions/genes more at the bottom of the array (as argsort
        #     works on a first come, first served basis).
        #
        #   - Sort the negated CRM scores (accessed in a random order) and
        #     return an array with indices that would sort those CRM scores
        #     from high to low (first position in the returned array
        #     contains the index to the value in
        #     crm_scores_with_ties_for_motif_numpy with the highest CRM
        #     score):
        #
        #       (-crm_scores_with_ties_for_motif_numpy)[random_permutations_to_break_ties_numpy].argsort()
        #
        #   - Undo the random order access of the array created in the
        #     previous step, so the indices that would sort
        #     crm_scores_with_ties_for_motif_numpy from high CRM scores to
        #     low CRM scores correspond again with the input array:
        #
        #       random_permutations_to_break_ties_numpy[
        #           (-crm_scores_with_ties_for_motif_numpy)[random_permutations_to_break_ties_numpy].argsort()
        #       ]
        #
        #   - Finally convert the array (previous step) which contains
        #     indices which would sort crm_scores_with_ties_for_motif_numpy
        #     from high CRM scores to low CRM scores and which would break
        #     tied scores in a fair (random) way to a ranking (int32):
        #
        #       ... .argsort().astype(np.int32)
        #
        rank_column_with_broken_ties_numpy = random_permutations_to_break_ties_numpy[
            (-crm_scores_with_ties_for_motif_numpy)[random_permutations_to_break_ties_numpy].argsort()
        ].argsort().astype(np.int32)

        return rank_column_with_broken_ties_numpy

    print(
        f'Create rankings from {args.feature_table_output_filename}.motifs_vs_{regions_or_genes_type}.scores.feather',
        file=sys.stderr
    )

    # Create feature table ranking.
    df_ranking__motifs_vs_regions_or_genes = df_scores__motifs_vs_regions_or_genes.apply(
        rank_CRM_scores_and_assign_random_ranking_in_range_for_ties_func,
        axis='index',
        raw=True
    )

    print(
        f'Write rankings of {regions_or_genes_type} for each motif in ' \
        f'motifs vs {regions_or_genes_type} style to ' \
        f'{args.feature_table_output_filename}.motifs_vs_{regions_or_genes_type}.rankings.feather',
        file=sys.stderr
    )

    df_ranking__motifs_vs_regions_or_genes.reset_index(inplace=True)
    df_ranking__motifs_vs_regions_or_genes.to_feather(
        path=f'{args.feature_table_output_filename}.motifs_vs_{regions_or_genes_type}.rankings.feather'
    )
    df_ranking__motifs_vs_regions_or_genes.set_index(regions_or_genes_type, inplace=True)

    del df_scores__motifs_vs_regions_or_genes

    print(
        f'Write rankings of {regions_or_genes_type} for each motif in ' \
        f'{regions_or_genes_type} vs motifs style to ' \
        f'{args.feature_table_output_filename}.{regions_or_genes_type}_vs_motifs.rankings.feather',
        file=sys.stderr
    )

    df_ranking__regions_or_genes_vs_motifs = df_ranking__motifs_vs_regions_or_genes.transpose()
    df_ranking__regions_or_genes_vs_motifs.reset_index(inplace=True)
    df_ranking__regions_or_genes_vs_motifs.to_feather(
        path=f'{args.feature_table_output_filename}.{regions_or_genes_type}_vs_motifs.rankings.feather'
    )


if __name__ == '__main__':
    main()
