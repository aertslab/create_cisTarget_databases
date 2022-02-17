#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Convert Feather cisTarget regions or genes vs motifs or tracks rankings database to sqlite3 (legacy)
               database.

Copyright (C): 2021 - Gert Hulselmans
"""

import argparse
import sqlite3

import polars as pl

CREATE_TABLE_STATEMENTS = r"""
DROP TABLE IF EXISTS rankings;
DROP TABLE IF EXISTS motifs;
CREATE TABLE rankings (geneID VARCHAR(255), ranking BLOB);
CREATE TABLE motifs (motifName VARCHAR(255), idx INTEGER);
"""

CREATE_INDEX_STATEMENT = r"CREATE UNIQUE INDEX id ON rankings (geneID)"
INSERT_FEATURE_STATEMENT = r"INSERT INTO motifs (idx, motifName) VALUES (?, ?);"
INSERT_RANKING_STATEMENT = r"INSERT INTO rankings (geneID, ranking) VALUES (?, ?);"


def convert_feather_db_to_sqlite3_db(feather_db_filename: str, sqlite3_db_filename: str):
    """
    Convert Feather cisTarget regions or genes vs motifs or tracks rankings database to sqlite3 (legacy) database.

    :param feather_db_filename: Feather cisTarget regions or genes vs motifs or tracks rankings database filename.
    :param sqlite3_db_filename:
    :return:
    """

    # Read Feather file in Polars dataframe.
    ct_db_df = pl.read_ipc(feather_db_filename)

    with sqlite3.connect(sqlite3_db_filename) as db:
        db.text_factory = str
        cursor = db.cursor()
        cursor.executescript(CREATE_TABLE_STATEMENTS)

        # Get feature ID column name.
        feature_ids_column_name = [
            column_name
            for column_name in ct_db_df.columns
            if column_name in ('motifs', 'tracks', 'features')
        ]

        # Get region or gene IDs.
        feature_ids = ct_db_df[feature_ids_column_name].to_numpy().flatten().tolist()

        # Get gene IDs or region IDs.
        gene_or_region_ids = [
            gene_or_region_id
            for gene_or_region_id in ct_db_df.columns
            if gene_or_region_id not in ('motifs', 'tracks', 'features')
        ]

        for feature_idx, feature_id in enumerate(feature_ids):
            cursor.execute(INSERT_FEATURE_STATEMENT, (feature_idx, feature_id))

        for gene_or_region_id in gene_or_region_ids:
            cursor.execute(INSERT_RANKING_STATEMENT,
                           (gene_or_region_id, ct_db_df[gene_or_region_id].to_numpy().tobytes()))

        cursor.execute(CREATE_INDEX_STATEMENT)
        cursor.close()


def main():
    parser = argparse.ArgumentParser(
        description='Convert cisTarget regions or genes vs motifs or tracks rankings Feather database to '
                    'legacy sqlite3 format.'
    )

    parser.add_argument(
        '-i',
        '--input',
        dest='ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_feather_filename',
        action='store',
        type=str,
        required=True,
        help='Input cisTarget regions or genes vs motifs or tracks rankings Feather database filename.'
    )

    parser.add_argument(
        '-o',
        '--output',
        dest='ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_sqlite3_filename',
        action='store',
        type=str,
        required=True,
        help='Output cisTarget regions or genes vs motifs or tracks rankings sqlite3 database filename.'
    )

    args = parser.parse_args()

    convert_feather_db_to_sqlite3_db(
        feather_db_filename=args.ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_feather_filename,
        sqlite3_db_filename=args.ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_sqlite3_filename
    )


if __name__ == '__main__':
    main()
