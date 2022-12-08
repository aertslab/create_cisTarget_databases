#!/usr/bin/env python3

"""
Purpose :      Convert Feather cisTarget regions or genes vs motifs or tracks rankings database to SQLite3 (legacy)
               database and vice versa.

Copyright (C): 2021-2022 - Gert Hulselmans
"""

import argparse
import sqlite3
from operator import itemgetter

import numpy as np
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

# SQL query that retrieves the ordered list of features in the database.
FEATURE_IDS_QUERY = r"SELECT motifName FROM motifs ORDER BY idx;"
# SQL query for retrieving the full list of genes scored in this database.
ALL_GENE_IDS_QUERY = r"SELECT geneID FROM rankings ORDER BY geneID;"
# SQL query for retrieving the whole database.
ALL_RANKINGS_QUERY = r"SELECT geneID, ranking FROM rankings ORDER BY geneID;"


def convert_feather_db_to_sqlite3_db(
    feather_db_filename: str, sqlite3_db_filename: str
) -> None:
    """
    Convert Feather cisTarget regions or genes vs motifs or tracks rankings database to sqlite3 (legacy) database.

    :param feather_db_filename: Feather cisTarget regions or genes vs motifs or tracks rankings database filename.
    :param sqlite3_db_filename: SQLite3 cisTarget regions or genes vs motifs or tracks rankings database filename.
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
            if column_name in ("motifs", "tracks", "features")
        ]

        # Get region or gene IDs.
        feature_ids = ct_db_df[feature_ids_column_name].to_numpy().flatten().tolist()

        # Get gene IDs or region IDs.
        gene_or_region_ids = [
            gene_or_region_id
            for gene_or_region_id in ct_db_df.columns
            if gene_or_region_id not in ("motifs", "tracks", "features")
        ]

        for feature_idx, feature_id in enumerate(feature_ids):
            cursor.execute(INSERT_FEATURE_STATEMENT, (feature_idx, feature_id))

        for gene_or_region_id in gene_or_region_ids:
            cursor.execute(
                INSERT_RANKING_STATEMENT,
                (gene_or_region_id, ct_db_df[gene_or_region_id].to_numpy().tobytes()),
            )

        cursor.execute(CREATE_INDEX_STATEMENT)
        cursor.close()


def convert_sqlite3_db_to_feather_db(
    sqlite3_db_filename: str, feather_db_filename: str
) -> None:
    """
    Convert sqlite3 (legacy) cisTarget regions or genes vs motifs or tracks rankings database to Feather database.

    :param sqlite3_db_filename: SQLite3 cisTarget regions or genes vs motifs or tracks rankings database filename.
    :param feather_db_filename: Feather cisTarget regions or genes vs motifs or tracks rankings database filename.
    :return:
    """

    with sqlite3.connect(sqlite3_db_filename) as db:
        cursor = db.cursor()

        # Get all gene or region IDS from the SQLite3 cisTarget database.
        region_or_gene_ids = list(
            map(itemgetter(0), cursor.execute(ALL_GENE_IDS_QUERY).fetchall())
        )
        # Get all motifs or tracks IDS from the SQLite3 cisTarget database.
        motif_or_track_ids = list(
            map(itemgetter(0), cursor.execute(FEATURE_IDS_QUERY).fetchall())
        )

        # Get dtype of numpy array buffers depending on the number of gene or region IDs.
        db_dtype = np.int32 if len(region_or_gene_ids) > 2**15 else np.int16

        # Read each gene or region rankings numpy buffer and convert
        # to a Polars Series and build a Polars DataFrame from those Series.
        ct_db_df = pl.DataFrame(
            [
                pl.Series(gene_or_region_id, np.frombuffer(ranking, dtype=db_dtype))
                for gene_or_region_id, ranking in cursor.execute(ALL_RANKINGS_QUERY)
            ]
        )

        cursor.close()

    # Assume the SQLite3 database contains genes unless all region_or_gene_ids contain
    # both ":" and "-" (e.g.: "chr1:1234-5678)".
    regions_or_genes = "genes"

    for region_or_gene_id in region_or_gene_ids:
        region_char1_pos = region_or_gene_id.find(":")

        if region_char1_pos >= 0:
            region_char2_pos = region_or_gene_id.find("-")

            if region_char2_pos >= 0:
                regions_or_genes = "regions"
            else:
                regions_or_genes = "genes"
                break
        else:
            regions_or_genes = "genes"
            break

    # Assume the SQLite3 database contains motifs if motif_or_track_ids contains a
    # feature that starts with "jaspar".
    motifs_or_tracks = "tracks"

    for motif_or_track_id in motif_or_track_ids:
        if motif_or_track_id.startswith("jaspar"):
            motifs_or_tracks = "motifs"

    # Add "motifs" or "tracks" column with motif and track IDs.
    ct_db_df.hstack(
        [pl.Series(motifs_or_tracks, motif_or_track_ids)],
        in_place=True,
    )

    # Write cisTarget Feather file.
    ct_db_df.write_ipc(
        f"{feather_db_filename}.{regions_or_genes}_vs_{motifs_or_tracks}.rankings.feather",
        compression="zstd",
    )


def main():
    parser = argparse.ArgumentParser(
        description="Convert cisTarget regions or genes vs motifs or tracks rankings Feather database to "
        "legacy SQLite3 format and vice versa."
    )

    parser.add_argument(
        "-f",
        "--feather",
        dest="ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_feather_filename",
        action="store",
        type=str,
        required=True,
        help="cisTarget regions or genes vs motifs or tracks rankings Feather database filename.",
    )

    parser.add_argument(
        "-s",
        "--sqlite3",
        dest="ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_sqlite3_filename",
        action="store",
        type=str,
        required=True,
        help="cisTarget regions or genes vs motifs or tracks rankings SQLite3 database filename.",
    )

    parser.add_argument(
        "-t",
        "--to",
        dest="convert_to",
        action="store",
        type=str,
        choices=["feather", "sqlite3"],
        required=True,
        default="sqlite3",
        help='Database type ("feather" or "sqlite3") to convert to. Default: "sqlite3".',
    )

    args = parser.parse_args()

    if args.convert_to == "sqlite3":
        convert_feather_db_to_sqlite3_db(
            feather_db_filename=args.ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_feather_filename,
            sqlite3_db_filename=args.ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_sqlite3_filename,
        )
    else:
        convert_sqlite3_db_to_feather_db(
            sqlite3_db_filename=args.ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_sqlite3_filename,
            feather_db_filename=args.ct_rankings_db_regions_or_genes_vs_motifs_or_tracks_feather_filename,
        )


if __name__ == "__main__":
    main()
