#!/usr/bin/env python3

"""
Purpose :      Create cisTarget track databases.

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

from bigwigaverageoverbed import (
    get_track_id_to_filename_dict,
    run_bigwig_average_over_bed_for_track,
)
from cistarget_db import (
    CisTargetDatabase,
    DatabaseTypes,
    MotifOrTrackIDs,
    MotifsOrTracksType,
    RegionOrGeneIDs,
)


def main():
    parser = argparse.ArgumentParser(description="Create cisTarget track databases.")

    parser.add_argument(
        "-b",
        "--bed",
        dest="bed_filename",
        action="store",
        type=str,
        required=True,
        help="BED filename which contains the regions/genes to score with bigWigAverageOverBed for each bigwig track "
        "(ChIP-seq) files.",
    )

    parser.add_argument(
        "-T",
        "--tracks_dir",
        dest="tracks_dir",
        action="store",
        type=str,
        required=True,
        help="Path to directory with bigwig track (ChIP-seq) files.",
    )

    parser.add_argument(
        "-d",
        "--tracks",
        dest="tracks_list_filename",
        action="store",
        type=str,
        required=True,
        help='Filename with list of track IDs to be scored from directory specified by "--tracks_dir".',
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="db_prefix",
        action="store",
        type=str,
        required=True,
        help="Feather database prefix output filename.",
    )

    parser.add_argument(
        "-a",
        "--bwaob",
        dest="bigwig_average_over_bed_path",
        action="store",
        type=str,
        required=False,
        default="bigWigAverageOverBed",
        help="Path to bigWigAverageOverBed "
        "(http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed). "
        'Default: "bigWigAverageOverBed".',
    )

    parser.add_argument(
        "-t",
        "--threads",
        dest="nbr_threads",
        action="store",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use when scoring tracks. Default: 1.",
    )

    parser.add_argument(
        "-p",
        "--partial",
        dest="partial",
        nargs=2,
        metavar=("CURRENT_PART", "NBR_TOTAL_PARTS"),
        action="store",
        type=int,
        required=False,
        help="Divide the tracks list in a number of total parts (of similar size) and score only the part defined by "
        "current_part. This allows creating partial databases on machines which do not have enough RAM to score "
        "all tracks in one iteration. This will only create a partial regions/genes vs tracks scoring database "
        "({db_prefix}.part_000{current_part}_of_000{nbr_total_parts}.regions_vs_tracks.scores.feather or "
        "{db_prefix}.part_000{current_part}_of_000{nbr_total_parts}.genes_vs_tracks.scores.feather).",
    )

    parser.add_argument(
        "-g",
        "--genes",
        dest="extract_gene_id_from_region_id_regex_replace",
        action="store",
        type=str,
        required=False,
        default=None,
        help="Take top score for a gene by taking the maximum score of multiple regions for that gene. "
        "Define a regex which will remove the non-gene part of the region ID, so only the gene ID remains. "
        'Examples: "gene_id#some_number": "#[0-9]+$" or "region_id@@gene_id": "^.+@@".',
    )

    parser.add_argument(
        "-s",
        "--seed",
        dest="seed",
        action="store",
        type=int,
        required=False,
        help="Random seed used for breaking ties when creating rankings for a range of tied scores. "
        "When setting this seed to a specific value and running this script with the same input, will result in "
        "the same rankings databases as output.",
    )

    parser.add_argument(
        "-r",
        "--ssh",
        dest="ssh_command",
        action="store",
        type=str,
        required=False,
        help="If defined, run bigWigAverageOverBed over ssh by running the provided command to make the connection before "
        "running bigWigAverageOverBed itself. "
        "Example: 'ssh -o ControlMaster=auto -o ControlPath=/tmp/ssh-control-path-%%l-%%h-%%p-%%r -o ControlPersist=600 <hostname>'",
    )

    args = parser.parse_args()

    if not os.path.exists(args.bed_filename):
        print(
            f'Error: BED region IDs filename "{args.bed_filename}" does not exist.',
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.exists(args.tracks_dir):
        print(
            f'Error: Track directory "{args.tracks_dir}" does not exist.',
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.exists(args.tracks_list_filename):
        print(
            f'Error: Tracks list filename "{args.tracks_list_filename}" does not exist.',
            file=sys.stderr,
        )
        sys.exit(1)

    if os.path.dirname(args.db_prefix) and not os.path.exists(
        os.path.dirname(args.db_prefix)
    ):
        print(
            f'Error: Parent directory "{os.path.dirname(args.db_prefix)}" for Feather database prefix output filename '
            "does not exist.",
            file=sys.stderr,
        )
        sys.exit(1)

    if args.partial:
        current_part, nbr_total_parts = args.partial

        if current_part < 1 or current_part > nbr_total_parts:
            print(
                f"Error: Current part ({current_part}) should be between 1 and the number of total parts "
                f"({nbr_total_parts}).",
                file=sys.stderr,
            )
            sys.exit(1)

        # Add info about which part of the database this wil be.
        db_prefix = f"{args.db_prefix}.part_{current_part:04d}_of_{nbr_total_parts:04d}"
    else:
        db_prefix = args.db_prefix

    # Get absolute path to bigWigAverageOverBed binary and see if it can be executed.
    bigwig_average_over_bed_path = shutil.which(args.bigwig_average_over_bed_path)

    if not bigwig_average_over_bed_path:
        print(
            f'Error: bigWigAverageOverBed binary ("{args.bigwig_average_over_bed_path}") could not be found or is not executable.'
        )
        sys.exit(1)

    # Set random seed to provided input value or a random integer.
    seed = args.seed if args.seed else random.randint(0, 2**64)

    # Get all region or gene IDs from the FASTA sequence names as a RegionOrGeneIDs object.
    region_or_gene_ids = RegionOrGeneIDs.get_region_or_gene_ids_from_bed(
        bed_filename=args.bed_filename,
        extract_gene_id_from_region_id_regex_replace=args.extract_gene_id_from_region_id_regex_replace,
    )

    # Get absolute path name for BED filename so in case bigWigAverageOverBed is ran over ssh, the BED file can be found.
    bed_filename = os.path.abspath(args.bed_filename)

    # Get track ID to track bigWig file name mapping for
    # a(n optionally) filtered list of track IDs:
    #   - if partial is set
    track_id_to_filename_dict = get_track_id_to_filename_dict(
        tracks_dir=os.path.abspath(args.tracks_dir),
        tracks_list_filename=args.tracks_list_filename,
        partial=(current_part, nbr_total_parts) if args.partial else None,
    )

    # Create MotifOrTracksIDs object from plain motif IDs.
    track_ids = MotifOrTrackIDs(
        motif_or_track_ids=set(track_id_to_filename_dict),
        motifs_or_tracks_type=MotifsOrTracksType.TRACKS,
    )

    nbr_region_or_gene_ids = len(region_or_gene_ids)
    nbr_tracks = len(track_id_to_filename_dict)

    if nbr_region_or_gene_ids == 0:
        print(f"Error: No {region_or_gene_ids.type.value} provided.", file=sys.stderr)
        sys.exit(1)

    if nbr_tracks == 0:
        print("Error: No tracks provided.", file=sys.stderr)
        sys.exit(1)

    print(
        f"Initialize dataframe ({nbr_region_or_gene_ids} {region_or_gene_ids.type.value} "
        f"x {nbr_tracks} tracks) for storing track scores for each {region_or_gene_ids.type.value} per track.",
        file=sys.stderr,
    )

    ct_scores_db_tracks_vs_regions_or_genes = CisTargetDatabase.create_db(
        db_type=DatabaseTypes.from_strings(
            scores_or_rankings="scores",
            column_kind="tracks",
            row_kind=region_or_gene_ids.type.value,
        ),
        region_or_gene_ids=region_or_gene_ids,
        motif_or_track_ids=track_ids,
        order="F",
    )

    def write_track_scores_for_track_to_ct_scores_db(
        track_id_and_df_track_scores: Tuple[str, pd.DataFrame]
    ) -> None:
        if (
            "nbr_of_scored_tracks"
            not in write_track_scores_for_track_to_ct_scores_db.__dict__
        ):
            write_track_scores_for_track_to_ct_scores_db.nbr_of_scored_tracks = 0

        track_id, df_track_scores = track_id_and_df_track_scores

        start_time = time.monotonic()
        ct_scores_db_tracks_vs_regions_or_genes.update_scores_for_motif_or_track(
            motif_or_track_id=track_id,
            df_scores_for_motif_or_track=df_track_scores["track_score"],
        )
        elapsed_time = time.monotonic() - start_time

        write_track_scores_for_track_to_ct_scores_db.nbr_of_scored_tracks += 1

        print(
            f"Adding bigWigAverageOverBed track scores ({write_track_scores_for_track_to_ct_scores_db.nbr_of_scored_tracks:d} of "
            f'{nbr_tracks:d}) for track "{track_id:s}" took {elapsed_time:0.6f} seconds.',
            file=sys.stderr,
        )

    def report_error(exception: BaseException) -> None:
        print(exception, file=sys.stderr)

    with mp.Pool(processes=args.nbr_threads) as pool:
        for track_id, track_filename in track_id_to_filename_dict.items():
            # Score all regions/genes in the FASTA file for the current track and write the result in the
            # ct_scores_db_tracks_vs_regions_or_genes CisTargetDatabase object.
            pool.apply_async(
                func=run_bigwig_average_over_bed_for_track,
                args=[
                    bigwig_average_over_bed_path,
                    bed_filename,
                    track_filename,
                    track_id,
                    args.extract_gene_id_from_region_id_regex_replace,
                    "max",
                    args.ssh_command,
                ],
                callback=write_track_scores_for_track_to_ct_scores_db,
                error_callback=report_error,
            )

        # Prevents any more task from being submitted to the pool.
        pool.close()

        # Wait for worker processes to exit.
        pool.join()

    if (
        "nbr_of_scored_tracks"
        not in write_track_scores_for_track_to_ct_scores_db.__dict__
    ):
        print(
            f"Error: None of {nbr_tracks:d} tracks were scored successfully.",
            file=sys.stderr,
        )
        sys.exit(1)
    elif (
        write_track_scores_for_track_to_ct_scores_db.nbr_of_scored_tracks != nbr_tracks
    ):
        print(
            f"Error: Only {write_track_scores_for_track_to_ct_scores_db.nbr_of_scored_tracks:d} out of {nbr_tracks:d} "
            f"tracks were scored successfully.",
            file=sys.stderr,
        )
        sys.exit(1)

    print("", file=sys.stderr)

    def write_db(ct_db: CisTargetDatabase, db_prefix: str):
        """
        Write cisTarget database to a Feather file and print database location and elapsed time.

        :param ct_db: cisTarget database object.
        :param db_prefix: Feather database file prefix.
        :return:
        """
        db_filename = ct_db.create_db_filename_from_db_prefix(
            db_prefix=db_prefix, extension="feather"
        )

        print(
            f"Writing cisTarget {ct_db.db_type.row_kind} vs {ct_db.db_type.column_kind} "
            f'{ct_db.db_type.scores_or_rankings} db: "{db_filename}"'
        )

        start_time = time.monotonic()
        ct_db.write_db(
            db_prefix=db_prefix,
            version=2,
        )
        elapsed_time = time.monotonic() - start_time

        print(
            f"Writing cisTarget {ct_db.db_type.row_kind} vs {ct_db.db_type.column_kind} "
            f"{ct_db.db_type.scores_or_rankings} db took: {elapsed_time:.06f} seconds\n"
        )

    if not args.partial:
        # Write cisTarget scores database (tracks vs regions or genes) to Feather file.
        write_db(ct_db=ct_scores_db_tracks_vs_regions_or_genes, db_prefix=db_prefix)

    # Create cisTarget scores database (regions or genes vs tracks) from (tracks vs regions or genes) version.
    ct_scores_db_regions_or_genes_vs_tracks = (
        ct_scores_db_tracks_vs_regions_or_genes.transpose()
    )

    # Write cisTarget scores database (regions or genes vs tracks) to Feather file.
    write_db(ct_db=ct_scores_db_regions_or_genes_vs_tracks, db_prefix=db_prefix)

    if not args.partial:
        # Create cisTarget rankings database (tracks vs regions or genes) from cisTarget scores database filename
        # (tracks vs regions or genes).
        print(
            f"""Create rankings from "{
                ct_scores_db_tracks_vs_regions_or_genes.create_db_filename_from_db_prefix(
                    db_prefix=db_prefix,
                    extension='feather'
                )
            }" with random seed set to {seed}.""",
            file=sys.stderr,
        )

        start_time = time.monotonic()
        ct_rankings_db_tracks_vs_regions_or_genes = (
            ct_scores_db_tracks_vs_regions_or_genes.convert_scores_db_to_rankings_db(
                seed=seed
            )
        )
        elapsed_time = time.monotonic() - start_time

        print(
            f"Creating cisTarget rankings db from cisTarget scores db took: "
            f"{elapsed_time:.06f} seconds\n"
        )

        # Reclaim memory occupied by cisTarget scores databases.
        del ct_scores_db_tracks_vs_regions_or_genes
        del ct_scores_db_regions_or_genes_vs_tracks

        # Do not write cisTarget rankings database (tracks vs regions or genes) to Feather file
        # as it can take a very long time to write it (1.5 hours for 1 million regions) as the
        # rankings database numpy array is in "C" order and writing a Feather database requires
        # traversing the numpy array in column order.
        # write_db(ct_db=ct_rankings_db_tracks_vs_regions_or_genes, db_prefix=db_prefix)

        # Create cisTarget rankings database (regions or genes vs tracks) from (tracks vs regions or genes) version.
        ct_rankings_db_regions_or_genes_vs_tracks = (
            ct_rankings_db_tracks_vs_regions_or_genes.transpose()
        )

        # Write cisTarget rankings database (regions or genes vs tracks) to Feather file.
        write_db(ct_db=ct_rankings_db_regions_or_genes_vs_tracks, db_prefix=db_prefix)


if __name__ == "__main__":
    main()
