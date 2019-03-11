#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Make feature table.

Copyright (C): 2019 - Gert Hulselmans
"""


import argparse

import io
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import multiprocessing as mp


def get_motif_name_to_filename_dict(motifs_dir, motifs_list_filename):
    motif_name_to_filename_dict = dict()

    with open(motifs_list_filename, 'r') as fh:
        for line in fh:
            motif_name=line.rstrip()

            if motif_name and not motif_name.startswith('#'):
                if motif_name.endswith('.cb'):
                    # Remove ".cb" extension from motif name.
                    motif_name = motif_name[:-3]

                motif_filename = os.path.join(motifs_dir, motif_name + '.cb')

                if not os.path.exists(motif_filename):
                    print(
                        'Error: Motif filename "{0:s}" does not exist.'.format(motif_filename),
                        file=sys.stderr
                    )
                    sys.exit(1)

                motif_name_to_filename_dict[motif_name] = motif_filename

    return motif_name_to_filename_dict


def get_sequence_names_from_fasta(fasta_filename):
    sequence_names_list = list()
    sequence_names_set = set()
    duplicated_sequences = False

    with open(fasta_filename, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                # Get sequence name by getting everything after '>' up till the first whitespace.
                sequence_name = line[1:].split(maxsplit=1)[0]

                # Check if all sequence names only appear once.
                if sequence_name in sequence_names_set:
                    print(
                        'Error: Sequence name "{0:s}" is not unique in FASTA file "{1:s}".'.format(
                            sequence_name,
                            fasta_filename
                        ),
                        file=sys.stderr
                    )
                    duplicated_sequences = True

                sequence_names_list.append(sequence_name)
                sequence_names_set.add(sequence_name)

    if duplicated_sequences:
        sys.exit(1)

    return sequence_names_list


def run_cluster_buster_for_motif(cluster_buster_path, fasta_filename, motif_filename, motif_name):
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

    return motif_name, crm_scores_df


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
        help='Path to motif directory.'
    )

    parser.add_argument(
        '-m',
        '--motifs',
        dest='motifs_list_filename',
        action='store',
        type=str,
        required=True,
        help='Filename with list of motif names relative to the directory specified by "--motifs".'
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

    if not os.path.exists(args.fasta_filename):
        print(
            'Error: Motifs list filename "{0:s}" does not exist.'.format(args.motifs_list_filename),
            file=sys.stderr
        )
        sys.exit(1)

    motif_name_to_filename_dict = get_motif_name_to_filename_dict(
        motifs_dir=args.motifs_dir,
        motifs_list_filename=args.motifs_list_filename
    )

    sequence_names = get_sequence_names_from_fasta(args.fasta_filename)

    # Create zeroed dataframe for all sequences vs all motif names.
    df_feature_table = pd.DataFrame(
        data=np.zeros((len(sequence_names),
                       len(motif_name_to_filename_dict)),
                      dtype=np.float32),
        index=sequence_names,
        columns=sorted(motif_name_to_filename_dict.keys())
    )

    nbr_motifs = len(motif_name_to_filename_dict)

    def add_crm_scores_for_motif_to_df_feature_table(motif_name_and_crm_scores_df):
        if 'nbr_of_scored_motifs' not in add_crm_scores_for_motif_to_df_feature_table.__dict__:
            add_crm_scores_for_motif_to_df_feature_table.nbr_of_scored_motifs = 0

        motif_name, crm_scores_df = motif_name_and_crm_scores_df
        df_feature_table.loc[crm_scores_df.index.tolist(), motif_name] = crm_scores_df['crm_score']

        add_crm_scores_for_motif_to_df_feature_table.nbr_of_scored_motifs += 1

        print(
            'Adding Cluster-Buster CRM scores ({0:d} of {1:d}) for motif "{2:s}".'.format(
                add_crm_scores_for_motif_to_df_feature_table.nbr_of_scored_motifs,
                nbr_motifs,
                motif_name
            ),
            file=sys.stderr
        )

    with mp.Pool(processes=args.nbr_threads) as pool:
        for motif_name, motif_filename in motif_name_to_filename_dict.items():
            pool.apply_async(func=run_cluster_buster_for_motif,
                             args=[args.cluster_buster_path,
                                   args.fasta_filename,
                                   motif_filename,
                                   motif_name],
                             callback=add_crm_scores_for_motif_to_df_feature_table)

        # Prevents any more task from being submitted to the pool.
        pool.close()

        # Wait for worker processes to exit.
        pool.join()

    def feature_table_write_tsv(df_feature_table, feature_table_output_filename, feature_table_output_format):
        """ Write feature table TSV file manually instead of with pandas. """

        # Get column names with all features.
        column_names = df_feature_table.columns.tolist()

        # Get row names with all regions.
        row_names = df_feature_table.index.tolist()

        # Get all CRM scores (numpy array).
        crm_scores = df_feature_table.get_values()

        # Write feature table TSV file.
        with open(feature_table_output_filename, 'w') as feature_table_output_fh:
            # Write header line.
            if feature_table_output_format == 'tsv':
                print('regions', end='\t', file=feature_table_output_fh)

            print('\t'.join(column_names), end='\n', file=feature_table_output_fh)

            # Write each row.
            for row_idx in range(len(row_names)):
                # Write region name.
                print(row_names[row_idx], end='\t', file=feature_table_output_fh)

                # Write CRM scores for all features for the current region.
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

        feature_table_write_tsv(df_feature_table=df_feature_table,
                                feature_table_output_filename=args.feature_table_output_filename,
                                feature_table_output_format=args.feature_table_output_format)
        # Faster than the following code:
        # df_feature_table.to_csv(
        #     path_or_buf=args.feature_table_output_filename,
        #     sep='\t',
        #     header=True,
        #     index=True,
        #     index_label="regions",
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

        feature_table_write_tsv(df_feature_table=df_feature_table,
                                feature_table_output_filename=args.feature_table_output_filename,
                                feature_table_output_format=args.feature_table_output_format)
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
