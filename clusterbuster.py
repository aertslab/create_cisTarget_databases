import io
import os
import subprocess
from typing import Optional, Dict, Union, Tuple

import numpy as np
import pandas as pd


def get_motif_id_to_filename_and_nbr_motifs_dict(
        motifs_dir: str,
        motifs_list_filename: str,
        partial: Optional[Tuple[int, int]] = None,
        min_nbr_motifs: Optional[int] = 1,
        max_nbr_motifs: Optional[int] = None,
        motif_md5_to_motif_id_filename: Optional[str] = None,
) -> (Dict[str, str], Dict[str, int]):
    """
    Create motif ID to Cluster-Buster motif file and number of motifs per motif file mapping.

    :param motifs_dir:
        Directory with Cluster-Buster motif files (with motif MD5 name or motif ID motif files).
    :param motifs_list_filename:
        File with Cluster-Buster motif MD5 names or motif IDs.
    :param partial: (current_part, nbr_total_parts)
        Divide the motif list in a number of total parts (of similar size) and return only the part defined by
        current_part. This makes it easier to create partial cisTarget motifs or tracks vs regions or genes scores
        database on machines which do not have enough RAM to score all motifs in one iteration, while still being able
        to give the same list of motifs to each instance.
    :param min_nbr_motifs:
        Only include motif IDs which have at least min_nbr_motifs of motifs in motif ID Cluster-Buster file.
        This allows to easily score e.g clustered Cluster-Buster motif files with at least 2 motifs as the singletons
        might already be scored before.
        Default: 1.
    :param max_nbr_motifs:
        Only include motif IDs which have at maximum max_nbr_motifs of motifs in motif ID Cluster-Buster file.
        Default: None.
    :param motif_md5_to_motif_id_filename: TSV file with motif MD5 names to motif IDs mapping (optional).
    :return: (motif_id_to_filename_dict, motif_id_to_nbr_motifs_dict):
        motif ID to CLuster-Buster motif filename mapping.
        motif ID to number of motifs in motif ID file mapping.
    """

    motif_id_to_filename_dict = dict()
    motif_id_to_nbr_motifs_dict = dict()
    motif_md5_to_motif_id_dict = dict()
    motif_id_to_motif_md5_dict = dict()

    if not min_nbr_motifs:
       min_nbr_motifs = 1

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
                        raise ValueError(
                            f'Error: Could not find motif MD5 name <=> motif ID association for "{motif_md5_or_id}".'
                        )

                    # Cluster-Buster motif MD5 name filename.
                    motif_filename = os.path.join(motifs_dir, motif_md5 + '.cb')
                else:
                    # No motif_md5_to_motif_id_filename was provided, so assume Cluster-Buster motif filenames in
                    # motifs_dir have motif IDs.

                    motif_id = motif_md5_or_id
                    # Cluster-Buster motif ID filename.
                    motif_filename = os.path.join(motifs_dir, motif_id + '.cb')

                if not os.path.exists(motif_filename):
                    raise OSError(
                        f'Error: Cluster-Buster motif filename "{motif_filename}" does not exist for motif {motif_id}.'
                    )
                else:
                    # Number of motifs in current motif ID Cluster-Buster file.
                    nbr_motifs_for_motif_id = 0

                    with open(motif_filename, 'r') as motif_fh:
                        # Count number of motifs in one motif file.
                        for motif_line in motif_fh:
                            if motif_line.startswith('>'):
                                nbr_motifs_for_motif_id += 1

                    # Filter out if motif ID Cluster-Buster file does have too little or too many motifs.
                    if nbr_motifs_for_motif_id >= min_nbr_motifs:
                        if max_nbr_motifs:
                            if nbr_motifs_for_motif_id <= max_nbr_motifs:
                                motif_id_to_filename_dict[motif_id] = motif_filename
                                motif_id_to_nbr_motifs_dict[motif_id] = nbr_motifs_for_motif_id
                        else:
                            motif_id_to_filename_dict[motif_id] = motif_filename
                            motif_id_to_nbr_motifs_dict[motif_id] = nbr_motifs_for_motif_id

    # Sort motif IDs by number of motifs in motif ID CLuster-Buster file (high to low) and then by motif ID name, so
    # motif IDs with a lot of motifs appear first, so the chance of having a few motif IDs with a lot of motifs only
    # running after all other ones are finished, should be much lower.
    motif_ids_sorted_by_nbr_motifs_for_motif_id = [
        motif_id
        for motif_id, nbr_motifs_for_motif_id in sorted(
            motif_id_to_nbr_motifs_dict.items(),
            key=lambda nbr_motifs_and_motif_id: (-nbr_motifs_and_motif_id[1], nbr_motifs_and_motif_id[0])
        )
    ]

    if partial:
        current_part, nbr_total_parts = partial

        if nbr_total_parts < 1:
            raise ValueError(f'"nbr_total_parts" ({nbr_total_parts}) of partial argument should be >= 1.')

        if current_part < 1:
            raise ValueError(f'"current_part" ({current_part}) of partial argument should be >= 1.')

        if current_part > nbr_total_parts:
            raise ValueError(
                f'"current_part" ({current_part}) of partial argument should be <= "nbr_total_parts" '
                f'({nbr_total_parts}).'
            )

        # Get partial motif IDs list for current requested part of the motif IDs.
        # If this function is run with a different current_part number, each partial motif IDs list should have
        # a similar number of motifs with a lot of motifs per motif ID Cluster-Buster file (in case a clustered motif
        # collection was given).
        partial_motif_ids_list = [
            motif_ids_sorted_by_nbr_motifs_for_motif_id[i]
            for i in range(current_part - 1, len(motif_id_to_nbr_motifs_dict), nbr_total_parts)
        ]

        # Recreate dictionaries with subset (corresponding to current_part) of motif IDs sorted by number of motifs in
        # motif ID Cluster-Buster file.
        motif_id_to_filename_dict = {
            motif_id: motif_id_to_filename_dict[motif_id]
            for motif_id in partial_motif_ids_list
        }

        motif_id_to_nbr_motifs_dict = {
            motif_id: motif_id_to_nbr_motifs_dict[motif_id]
            for motif_id in partial_motif_ids_list
        }
    else:
        # Recreate dictionaries with motif IDs sorted by number of motifs in motif ID Cluster-Buster file (high to low).
        motif_id_to_filename_dict = {
            motif_id: motif_id_to_filename_dict[motif_id]
            for motif_id in motif_ids_sorted_by_nbr_motifs_for_motif_id
        }

        motif_id_to_nbr_motifs_dict = {
            motif_id: motif_id_to_nbr_motifs_dict[motif_id]
            for motif_id in motif_ids_sorted_by_nbr_motifs_for_motif_id
        }

    return motif_id_to_filename_dict, motif_id_to_nbr_motifs_dict


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
        names=['seq_name', 'crm_score'],
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
