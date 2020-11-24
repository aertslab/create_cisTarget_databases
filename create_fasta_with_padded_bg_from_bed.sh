#!/bin/bash
#
# Copyright (C) 2020 - Gert Hulselmans
#
# Purpose:
#   Create FASTA file for usage with:
#     - Cluster-Buster "-b X"
#     - create_cistarget_motif_databases.py "-b X"
#   to extend regions with X bp in both directions to use the extended part
#   only for creating the background model while only scoring the original
#   region for the motif.
#   This allows scoring small sequences better as they will be less affected
#   by their too local small background.



create_fasta_with_padded_bg_from_bed () {
    local genomic_fasta_filename="${1}";
    local chrom_sizes_filename="${2}";
    local bed_regions_filename="${3}";
    local fasta_padded_bg_filename="${4}";
    local -i bg_padding="${5}";
    local create_name="${6:-0}";

    if [ ${#@} -lt 5 ] ; then
        printf 'Usage: create_fasta_with_padded_bg_from_bed genomic_fasta_filename chrom_sizes_filename bed_regions_filename fasta_padded_bg_filename bg_padding create_name\n\n';
        printf 'Arguments:\n';
        printf '    genomic_fasta_filename:   Genomic FASTA file.\n';
        printf '    chrom_sizes_filename:     File with chromosome names and their sizes.\n';
        printf '    bed_regions_filename:     File with regions in BED format which will\n';
        printf '                              be extended in both directions with X bp to\n';
        printf '                              create padded background FASTA file.\n';
        printf '    fasta_padded_bg_filename: FASTA output file with padded background.\n';
        printf '    bg_padding:               Extend each region in the BED file with X bp\n';
        printf '                              to create FASTA padded background file.\n';
        printf '    create_name:              Create region name based on the original BED\n';
        printf '                              region coordinates and use it for naming the\n';
        printf '                              sequences in the FASTA file:\n';
        printf '                                 - yes or 1: use cooridinates to create name\n';
        printf '                                 - no or 0: use original region name\n\n';
        printf 'Purpose:\n';
        printf '    Create FASTA file for usage with Cluster-Buster "-b X" or\n'
        printf '    create_cistarget_motif_databases.py "-b X" to extend regions with X bp in\n';
        printf '    both directions to use the extended part only for creating the background\n';
        printf '    model while only scoring the original region for the motif.\n';
        printf '    This allows scoring small sequences better as they will be less affected\n';
        printf '    by their too local small background.\n\n';
        return $?;
    fi

    if [ "${create_name}" = "yes" ] ; then
        create_name=1;
    fi

    bedtools slop \
        -i "${bed_regions_filename}" \
        -g "${chrom_sizes_filename}" \
        -b "${bg_padding}" \
      | paste - <(cut -f 2-3 "${bed_regions_filename}") \
      | awk \
            -F '\t' \
            -v 'OFS=\t' \
            -v create_name="${create_name}" \
            -v bg_padding="${bg_padding}" \
            -v create_name="${create_name}" \
            '
            {
                # Column number which contains start and end coordinate from original BED file.
                original_start_idx = NF - 1;
                original_end_idx = NF;

                if (create_name == 1 || original_start_idx == 4) {
                    # Original BED file does not contain a name column,
                    # so construct a name based on the original chromosomal positions.
                    name = $1 ":" $original_start_idx "-" $original_end_idx;
                } else {
                    name = $4;
                }

                # Calculate how many Ns to add in front and end of the sequence in case
                # a chromosome boundary was crossed when extending the original regions
                # with the background padding:
                #     name add_nbr_Ns_start_seq add_nbr_Ns_end_seq
                name_with_nbr_Ns = name " " bg_padding - ($original_start_idx - $2) " " bg_padding - ($3 - $original_end_idx);

                # Print cooridinates of BED file created with "bedtools slop" and
                # region name with number of Ns to add in front and end of the sequence.
                print $1, $2, $3, name_with_nbr_Ns;
            }
            ' \
      | bedtools getfasta \
            -nameOnly \
            -fi "${genomic_fasta_filename}" \
            -bed - \
      | awk \
            -F ' ' \
            '
            # Repeat character c, n times.
            function repeat(c, n,    t) {
                if (n == 0) {
                    return "";
                }

                # Create temporary string of n spaces.
                t = sprintf("%" n "s", "");

                # Replace each space with character of interest.
                gsub(/ /, c, t);

                return t
            }

            {
                if ($1 ~ /^>/) {
                    # Extract number of Ns to add in front and end from the last 2 values from the sequence name string.
                    add_nbr_Ns_start_seq = $(NF - 1);
                    add_nbr_Ns_end_seq = $NF;

                    # Discard the last 2 columns and print the sequence name only.
                    NF = NF - 2;
                    print $0;
                } else {
                    # Add padding with Ns in front and end of the sequence.
                    print repeat("N", add_nbr_Ns_start_seq) $0 repeat("N", add_nbr_Ns_end_seq);
                }
            }
            ' \
      > "${fasta_padded_bg_filename}";
}



create_fasta_with_padded_bg_from_bed "${@}";
