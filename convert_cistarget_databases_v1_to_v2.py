#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose :      Convert cisTarget Feather database from Feather v1 to v2 format (with or without compression) and vice
               versa.

Copyright (C): 2022 - Gert Hulselmans
"""

import argparse

from typing import Optional

import pyarrow.compute as pc
import pyarrow.feather as pf


def convert_feather_v1_to_v2_vice_versa(
        input_ct_db_filename: str,
        output_ct_db_filename: str,
        compression: Optional[str] = "zstd",
        compression_level: int = 6,
        to_version: int = 2,
):
    """
    Convert cisTarget Feather database from Feather v1 to v2 format (with or without compression) and vice versa.

    :param input_ct_db_filename: input cisTarget database filename.
    :param output_ct_db_filename: output cisTarget database filename.
    :param compression: Compression method: "zstd" (default), "lz4" or "uncompressed".
    :param compression_level: Compression level for "zstd" or "lz4".
    :param to_version: Output Feather file version: 1 (legacy) or 2 (default).
    :return:
    """

    if to_version != 2 or to_version != 1:
        raise ValueError("Feather file version only supports 1 or 2 (default).")

    if to_version == 1:
        # Compression is not supported in Feather v1 format.
        compression = "uncompressed"
        compression_level = None

    if compression not in {"zstd", "lz4", "uncompressed"}:
        raise ValueError(
            f'Unsupported compression value "{compression}". Choose "zstd" (default), "lz4" or "uncompressed".'
        )

    # Read input cisTarget database as a pyarrow Table.
    df_pa_table = pf.read_table(
        source=input_ct_db_filename,
    )

    # Get all column names.
    all_column_names = df_pa_table.column_names

    try:
        # Check if we have an old database that still used a "features" column and rename it.
        features_idx = all_column_names.index("features")

        # Get column which contains motif or track names.
        motifs_or_track_names = df_pa_table.column(features_idx)

        if pc.sum(pc.starts_with(motifs_or_track_names, "jaspar")).as_py() > 0:
            # It is a motif vs genes/regions database if JASPAR motif names were found in the "features" column.
            all_column_names[features_idx] = "motifs"
        else:
            all_column_names[features_idx] = "tracks"

        df_pa_table.drop(["features"])
        # Rename features column in database to "motifs" or "tracks".
        df_pa_table = df_pa_table.rename_columns(all_column_names)
    except ValueError:
        # No old database (with "features" column).
        pass

    # Get database index column ("motifs", "tracks", "regions" or "genes" depending of the database type).
    for column_name in all_column_names:
        if column_name in {"motifs", "tracks", "regions", "genes"}:
            index_column = df_pa_table.column(features_idx)
            break

    # Sort column names (non-index columns) and add index column as last column.
    column_names_sorted_and_index = sorted(
        [
            column_name
            for column_name in all_column_names
            if column_name not in index_column._name
        ]
    )
    column_names_sorted_and_index.append(index_column._name)

    # Create a new pyarrow Table with columns in the new order.
    df_pa_table = df_pa_table.select(column_names_sorted_and_index)

    # Writhe cisTarget database to a new Feather file with the requested compression/version settings.
    pf.write_feather(
        df=df_pa_table,
        dest=output_ct_db_filename,
        compression=compression,
        compression_level=compression_level,
        version=to_version
    )


def main():
    parser = argparse.ArgumentParser(
        description="Convert cisTarget Feather database from Feather v1 to v2 format (with or without "
                    "compression) and vice versa."
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input_ct_db_filename",
        action="store",
        type=str,
        required=True,
        help="Input cisTarget Feather database filename."
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_ct_db_filename",
        action="store",
        type=str,
        required=True,
        help="Output cisTarget Feather database filename."
    )

    parser.add_argument(
        "-c",
        "--compression",
        dest="compression",
        action="store",
        type=str,
        choices=["zstd", "lz4", "uncompressed"],
        required=False,
        default="zstd",
        help='Compression method for output cisTarget Feather database: "zstd" (default), "lz4" or "uncompressed".'
    )

    parser.add_argument(
        "-l",
        "--level",
        dest="compression_level",
        action="store",
        type=str,
        choices=["zstd", "lz4", "uncompressed"],
        required=False,
        default="zstd",
        help='Compression method for output cisTarget Feather database: "zstd" (default), "lz4" or "uncompressed".'
    )

    parser.add_argument(
        "-v",
        "--version",
        dest="to_version",
        action="store",
        type=int,
        required=False,
        default=6,
        help="Compression level for zstd or lz4. Default: 6."
    )

    args = parser.parse_args()

    convert_feather_v1_to_v2_vice_versa(
        input_ct_db_filename=args.input_ct_db_filename,
        output_ct_db_filename=args.output_ct_db_filename,
        compression=args.compression,
        compression_level=args.compression_level,
        to_version=args.to_version,
    )


if __name__ == "__main__":
    main()
