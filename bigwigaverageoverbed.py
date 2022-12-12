import io
import os
import subprocess
from typing import Dict, Optional, Tuple, Union

import numpy as np
import pandas as pd


def get_track_id_to_filename_dict(
    tracks_dir: str,
    tracks_list_filename: str,
    partial: Optional[Tuple[int, int]] = None,
) -> (Dict[str, str], Dict[str, int]):
    """
    Create track ID to bigWig track file mapping.

    :param tracks_dir:
        Directory with bigWig track files.
    :param tracks_list_filename:
        File with track IDs.
    :param partial: (current_part, nbr_total_parts)
        Divide the track list in a number of total parts (of similar size) and return only the part defined by
        current_part. This makes it easier to create partial cisTarget motifs or tracks vs regions or genes scores
        database on machines which do not have enough RAM to score all motifs in one iteration, while still being able
        to give the same list of tracks to each instance.
    :return: track_id_to_filename_dict:
        track ID to bigWig filename mapping.
    """

    track_id_to_filename_dict = dict()

    # Create track ID to bigWig track filename mapping.
    with open(tracks_list_filename, "r") as fh:
        for line in fh:
            track_id = line.rstrip()

            if track_id and not track_id.startswith("#"):
                if track_id.endswith(".bw"):
                    track_filename = os.path.join(tracks_dir, track_id)

                    # Remove ".bw" extension from track ID.
                    track_id = track_id[:-3]
                elif track_id.lower().endswith(".bigwig"):
                    track_filename = os.path.join(tracks_dir, track_id)

                    # Remove ".bigwig" extension from track ID.
                    track_id = track_id[:-7]
                else:
                    # If we got a track ID, try to find the corresponding bigWig file.
                    track_filenames = [
                        os.path.join(tracks_dir, track_id + bigwig_suffix)
                        for bigwig_suffix in [".bw", ".bigWig", "bigwig"]
                    ]

                    for track_filename in track_filenames:
                        if os.path.exists(track_filename):
                            break

                if not os.path.exists(track_filename):
                    raise OSError(
                        f'Error: bigWig track ID filename "{track_filename}" does not exist for track {track_id}.'
                    )
                else:
                    track_id_to_filename_dict[track_id] = track_filename

    # Sort track IDs.
    track_ids_sorted = sorted(track_id_to_filename_dict.keys())

    if partial:
        current_part, nbr_total_parts = partial

        if nbr_total_parts < 1:
            raise ValueError(
                f'"nbr_total_parts" ({nbr_total_parts}) of partial argument should be >= 1.'
            )

        if current_part < 1:
            raise ValueError(
                f'"current_part" ({current_part}) of partial argument should be >= 1.'
            )

        if current_part > nbr_total_parts:
            raise ValueError(
                f'"current_part" ({current_part}) of partial argument should be <= "nbr_total_parts" '
                f"({nbr_total_parts})."
            )

        # Get partial track IDs list for current requested part of the track IDs.
        # If this function is run with a different current_part number, each partial track IDs list should have
        # a similar number of tracks.
        partial_track_ids_list = [
            track_ids_sorted[i]
            for i in range(
                current_part - 1, len(track_id_to_filename_dict), nbr_total_parts
            )
        ]

        # Recreate dictionaries with track IDs sorted by track ID.
        track_id_to_filename_dict = {
            track_id: track_id_to_filename_dict[track_id]
            for track_id in partial_track_ids_list
        }
    else:
        # Recreate dictionaries with track IDs sorted by track ID.
        track_id_to_filename_dict = {
            track_id: track_id_to_filename_dict[track_id]
            for track_id in track_ids_sorted
        }

    return track_id_to_filename_dict


def run_bigwig_average_over_bed_for_track(
    bigwig_average_over_bed_path: str,
    bed_filename: str,
    bigwig_filename: str,
    track_id: str,
    extract_gene_id_from_region_id_regex_replace: Optional[str] = None,
    scoring_column: str = "max",
    ssh_command: Optional[Union[str, list]] = None,
) -> Tuple[str, pd.DataFrame]:
    """
    Score each region in the BED file with bigWigAverageOverBed over the bigWig file
    and only keep the top track score per region ID/gene ID.

    :param bigwig_average_over_bed_path: Path to bigWigAverageOverBed binary.
    :param bed_filename:        BED filename with regions to score.
    :param bigwig_filename:     BigWig filename from which to extract the score for each region in the BED filename.
    :param track_id:            Track ID.
    :param extract_gene_id_from_region_id_regex_replace:
                                Define a regex which will remove the non-gene part of the region ID of each sequence
                                name in the FASTA file, so only the gene ID remains. If set to None the whole region ID
                                will be kept instead. In case of region IDs, the best track score per region is kept.
                                In case of gene IDs, the best track score from multiple regions is kept.
    :param scoring_column:      Scoring column to use from bigWigAverageOverBed output (default: "max").
    :param ssh_command:         If defined, run bigWigAverageOverBed over ssh by running the provided command to make the
                                connection before running Cluster-Buster.
                                Example : 'ssh -o ControlMaster=auto -o ControlPath=/tmp/ssh-control-path-%l-%h-%p-%r -o ControlPersist=600 <hostname>'
    :return:                    (track_id, df_track_scores): track ID and dataframe with top track score per region/gene ID.
    """

    bigwig_average_over_bed_command = []

    if ssh_command:
        # Add SSH command to the start of the Cluster-Buster command.
        if isinstance(ssh_command, str):
            bigwig_average_over_bed_command.extend(ssh_command.split())
        elif isinstance(ssh_command, list):
            bigwig_average_over_bed_command.extend(ssh_command)

    # Construct Cluster-Buster command line.
    bigwig_average_over_bed_command.extend(
        [
            bigwig_average_over_bed_path,
            "-minMax",
            bigwig_filename,
            bed_filename,
            "/dev/stdout",
        ]
    )

    # Score each region in FASTA file with Cluster-Buster for the provided motif and get top CRM score for each region.
    try:
        pid = subprocess.Popen(
            args=bigwig_average_over_bed_command,
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
            creationflags=0,
        )
        stdout_data, stderr_data = pid.communicate()
    except OSError as msg:
        raise RuntimeError(
            "Execution error for: '"
            + " ".join(bigwig_average_over_bed_command)
            + "': "
            + str(msg)
        )

    if pid.returncode != 0:
        raise RuntimeError(
            "Error: Non-zero exit status for: '" + " ".join(bigwig_average_over_bed_command) + "'"
        )

    # Read bigWigAverageOverBed standard out as a pandas dataframe.
    df_track_scores = pd.read_csv(
        filepath_or_buffer=io.BytesIO(stdout_data),
        sep="\t",
        header=None,
        names=[
            "name",     # name field from bed, which should be unique
            "size",     # size of bed (sum of exon sizes
            "covered",  # bases within exons covered by bigWig
            "sum",      # sum of values over all bases covered
            "mean0",    # average over bases with non-covered bases counting as zeroes
            "mean",     # average over just covered bases
            "min",      # minimum coverage value
            "max",      # maximum coverage value
        ],
        index_col="name",
        usecols=["name", scoring_column],
        dtype={"name": str, scoring_column: np.float32},
        engine="c",
    )

    df_track_scores.rename(columns={scoring_column: "track_score"}, inplace=True)

    if extract_gene_id_from_region_id_regex_replace:
        # Extract gene ID from the region ID by removing the non-gene part.
        #
        # Take the top track score for each gene ID by taking the maximum track
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
        df_track_scores = (
            df_track_scores.assign(
                gene_ids=df_track_scores.index.str.replace(
                    extract_gene_id_from_region_id_regex_replace, "", regex=True
                )
            )
            .groupby("gene_ids")
            .max()
        )

    return track_id, df_track_scores
