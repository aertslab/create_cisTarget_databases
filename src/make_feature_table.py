#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Make feature table.

Copyright (C): 2019-2020 - Gert Hulselmans
"""


import argparse

import io
import os
import re
import sys
import subprocess
import numpy as np
import pandas as pd
import multiprocessing as mp


def get_motif_id_to_filename_dict(motifs_dir, motifs_list_filename, motif_md5_to_motif_id_filename=None):
    motif_id_to_filename_dict = dict()
    motif_md5_to_motif_id_dict = dict()
    motif_id_to_motif_md5_dict = dict()

    # Create motif MD5 names to motif ID names mapping and vice versa
    # if a motif MD5 to motif ID file was provided. This will be used
    # later to map Cluster-Buster motif files in motifs_dir with motif
    # MD5 names to motif IDs.
    if motif_md5_to_motif_id_filename:
        with open(motif_md5_to_motif_id_filename, 'r') as fh:
            for line in fh:
                line = line.rstrip()

                if line and not line.startswith('#'):
                    motif_md5, motif_id = line.rstrip().split('\t')[0:2]

                    # Store motif MD5 to motif ID mapping as vice versa.
                    motif_md5_to_motif_id_dict[motif_md5] = motif_id
                    motif_id_to_motif_md5_dict[motif_id] = motif_md5

    # Create motif ID to Cluster-Buster motif filename mapping.
    with open(motifs_list_filename, 'r') as fh:
        for line in fh:
            motif_md5_or_id = line.rstrip()

            if motif_md5_or_id and not motif_md5_or_id.startswith('#'):
                if motif_md5_or_id.endswith('.cb'):
                    # Remove ".cb" extension from motif ID.
                    motif_md5_or_id = motif_md5_or_id[:-3]

                if motif_md5_to_motif_id_dict:
                    if motif_md5_or_id in motif_md5_to_motif_id_dict:
                        # Get associated motif ID for motif MD5.
                        motif_id = motif_md5_to_motif_id_dict[motif_md5_or_id]
                        motif_md5 = motif_md5_or_id
                    elif motif_md5_or_id in motif_id_to_motif_md5_dict:
                        # Get associated motif MD5 for motif ID.
                        motif_id = motif_md5_or_id
                        motif_md5 = motif_id_to_motif_md5_dict[motif_md5_or_id]
                    else:
                        print(
                            f'Error: Could not find motif MD5 <=> motif ID association for "{motif_md5_or_id}".',
                            file=sys.stderr
                        )
                        sys.exit(1)

                    # Cluster-Buster motif MD5 filename.
                    motif_filename = os.path.join(motifs_dir, motif_md5 + '.cb')
                else:
                    motif_id = motif_md5_or_id
                    motif_filename = os.path.join(motifs_dir, motif_id + '.cb')

                if not os.path.exists(motif_filename):
                    print(
                        f'Error: Cluster-Buster motif filename "{motif_filename}" does not exist for motif {motif_id}.',
                        file=sys.stderr
                    )
                    sys.exit(1)

                motif_id_to_filename_dict[motif_id] = motif_filename

    return motif_id_to_filename_dict


def get_region_ids_or_gene_ids_from_fasta(fasta_filename, extract_gene_id_from_region_id_regex_replace):
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


def run_cluster_buster_for_motif(cluster_buster_path, fasta_filename, motif_filename, motif_id, extract_gene_id_from_region_id_regex_replace=None):
    # Score each region in FASTA file with Cluster-Buster
    # for motif and get top CRM score for each region.
    clusterbuster_command = [cluster_buster_path,
                             '-f', '4',
                             '-c', '0.0',
                             '-r', '10000',
                             '-t', '1',
                             motif_filename,
                             fasta_filename]

    try:
        pid = subprocess.Popen(args=clusterbuster_command,
                               bufsize=1,
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
        print("\nExecution error for: '" + ' '.join(clusterbuster_command) + "': " + str(msg),
              file=sys.stderr)
        sys.exit(1)

    if pid.returncode != 0:
        print("\nError: Non-zero exit status for: " + ' '.join(clusterbuster_command) + "'",
              file=sys.stderr)
        sys.exit(1)

    crm_scores_df = pd.read_csv(
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
        crm_scores_df = crm_scores_df.assign(
            gene_ids=crm_scores_df.index.str.replace(
                extract_gene_id_from_region_id_regex_replace,
                '',
                regex=True
            )
        ).groupby('gene_ids').max()

    return motif_id, crm_scores_df


def main():
    parser = argparse.ArgumentParser(
        description='Make feature table.'
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
        '-O',
        '--output-format',
        dest='feature_table_output_format',
        action='store',
        type=str,
        required=False,
        choices=['feather', 'tsv', 'tsv_for_R'],
        default='feather',
        help='Feature table output format. Default: "feather" (fast).'
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

    # Create zeroed dataframe for all region IDs or gene IDs vs all motif IDs.
    df_feature_table = pd.DataFrame(
        data=np.zeros((len(region_ids_or_gene_ids),
                       len(motif_id_to_filename_dict)),
                      dtype=np.float32),
        index=region_ids_or_gene_ids,
        columns=sorted(motif_id_to_filename_dict.keys())
    )

    # Add index name: 'regions' or 'genes'.
    df_feature_table.rename_axis(
        index=regions_or_genes_type,
        axis='index',
        copy=False,
        inplace=True
    )

    # Add column axis name: 'motifs'.
    df_feature_table.rename_axis(
        columns='motifs',
        axis='columns',
        copy=False,
        inplace=True
    )

    nbr_motifs = len(motif_id_to_filename_dict)

    def add_crm_scores_for_motif_to_df_feature_table(motif_id_and_crm_scores_df):
        if 'nbr_of_scored_motifs' not in add_crm_scores_for_motif_to_df_feature_table.__dict__:
            add_crm_scores_for_motif_to_df_feature_table.nbr_of_scored_motifs = 0

        motif_id, crm_scores_df = motif_id_and_crm_scores_df

        df_feature_table.loc[crm_scores_df.index.tolist(), motif_id] = crm_scores_df['crm_score']

        add_crm_scores_for_motif_to_df_feature_table.nbr_of_scored_motifs += 1

        print(
            'Adding Cluster-Buster CRM scores ({0:d} of {1:d}) for motif "{2:s}".'.format(
                add_crm_scores_for_motif_to_df_feature_table.nbr_of_scored_motifs,
                nbr_motifs,
                motif_id
            ),
            file=sys.stderr
        )

    with mp.Pool(processes=args.nbr_threads) as pool:
        for motif_id, motif_filename in motif_id_to_filename_dict.items():
            pool.apply_async(
                func=run_cluster_buster_for_motif,
                args=[
                    args.cluster_buster_path,
                    args.fasta_filename,
                    motif_filename,
                    motif_id,
                    args.extract_gene_id_from_region_id_regex_replace
                ],
                callback=add_crm_scores_for_motif_to_df_feature_table
            )

        # Prevents any more task from being submitted to the pool.
        pool.close()

        # Wait for worker processes to exit.
        pool.join()

    def feature_table_write_tsv(df_feature_table, feature_table_output_filename, feature_table_output_format, regions_or_genes_type):
        """ Write feature table TSV file manually instead of with pandas. """

        # Get column names with all features.
        column_names = df_feature_table.columns.tolist()

        # Get row names with all region IDs or gene IDs.
        row_names = df_feature_table.index.tolist()

        # Get all CRM scores (numpy array).
        crm_scores = df_feature_table.get_values()

        # Write feature table TSV file.
        with open(feature_table_output_filename, 'w') as feature_table_output_fh:
            # Write header line.
            if feature_table_output_format == 'tsv':
                print(regions_or_genes_type, end='\t', file=feature_table_output_fh)

            print('\t'.join(column_names), end='\n', file=feature_table_output_fh)

            # Write each row.
            for row_idx in range(len(row_names)):
                # Write region ID or gene ID.
                print(row_names[row_idx], end='\t', file=feature_table_output_fh)

                # Write CRM scores for all features for the current region or gene.
                crm_scores[row_idx].tofile(file=feature_table_output_fh, sep='\t')

                # Write newline.
                print(end='\n', file=feature_table_output_fh)

    if args.feature_table_output_format == 'tsv':
        print(
            'Write feature table result TSV table: "{0:s}".'.format(
                args.feature_table_output_filename
            ),
            file=sys.stderr
        )

        feature_table_write_tsv(
            df_feature_table=df_feature_table,
            feature_table_output_filename=args.feature_table_output_filename,
            feature_table_output_format=args.feature_table_output_format,
            regions_or_genes_type=regions_or_genes_type
        )
        # Faster than the following code:
        # df_feature_table.to_csv(
        #     path_or_buf=args.feature_table_output_filename,
        #     sep='\t',
        #     header=True,
        #     index=True,
        #     index_label=regions_or_genes_type,
        #     quoting=None,
        #     chunksize=1000,
        #)
    elif args.feature_table_output_format == 'tsv_for_R':
        print(
            'Write feature table result TSV table for R: "{0:s}".'.format(
                args.feature_table_output_filename
            ),
            file=sys.stderr
        )

        feature_table_write_tsv(
            df_feature_table=df_feature_table,
            feature_table_output_filename=args.feature_table_output_filename,
            feature_table_output_format=args.feature_table_output_format,
            regions_or_genes_type=regions_or_genes_type
        )
        # Faster than the following code:
        # df_feature_table.to_csv(
        #     path_or_buf=args.feature_table_output_filename,
        #     sep='\t',
        #     header=True,
        #     index=True,
        #     index_label=False,
        #     quoting=None,
        #     chunksize=1000,
        # )
    elif args.feature_table_output_format == 'feather':
        print(
            'Write feature table result table in feather format: "{0:s}".'.format(
                args.feature_table_output_filename
            ),
            file=sys.stderr
        )
        df_feature_table.reset_index(inplace=True)
        df_feature_table.to_feather(
            fname=args.feature_table_output_filename
        )


if __name__ == '__main__':
    main()
