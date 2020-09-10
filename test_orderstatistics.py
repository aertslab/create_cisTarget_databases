import numpy as np
import pytest

from orderstatistics import _calculate_cross_species_rank_ratio_with_order_statistics, create_cross_species_ranking_for_motif


@pytest.fixture
def motif_id_rankings_per_species():
    # Numpy array with motif rankings:
    #   - per row: ranking for motif for a certain region/gene.
    #   - per column: ranking for motif for a certain species.
    motif_id_rankings_per_species = np.array(
        [[ 3, 17, 20,  6, 13,  0, 11, 22],
         [16,  3, 17,  5, 28,  8, 28,  1],
         [ 9, 28,  0,  7, 24, 28, 23,  4],
         [21, 23, 21, 16, 18, 15,  9,  9],
         [11,  2, 15,  3, 20, 18, 12, 15],
         [12,  5, 28, 11, 26, 14, 26, 26],
         [ 7, 21,  9,  0, 23, 25, 20, 14],
         [ 4,  4, 14, 25, 25,  6, 24, 29],
         [ 5, 27, 18, 19, 14, 22, 13, 11],
         [ 0, 11, 10, 10, 16,  2,  3, 24],
         [14,  6, 25,  4,  6, 27, 29, 18],
         [19, 26,  6, 23, 17, 11, 19,  8],
         [20, 18, 27,  1,  4, 21, 25,  2],
         [24, 24, 13, 15, 21,  5,  0,  7],
         [13, 15, 23, 18,  2, 29,  7,  5],
         [ 2, 12, 24,  8,  8, 20, 16, 13],
         [26,  9,  2,  2,  7,  4,  8, 28],
         [23, 19,  3, 14, 15, 23,  1, 23],
         [27, 22,  1, 13,  1, 13,  4, 27],
         [25, 13, 16, 20, 19,  7, 22, 10],
         [10, 10,  7, 17, 10, 10, 18, 19],
         [28, 16,  8, 21,  3,  1, 14,  6],
         [17,  1, 19, 24,  0, 16,  6,  3],
         [ 8, 29,  5, 29,  5, 12, 21,  0],
         [29, 20, 22, 22,  9, 26, 17, 12],
         [ 1,  8, 26,  9, 29,  3, 15, 16],
         [18, 25, 11, 12, 11, 19, 10, 20],
         [15,  0,  4, 28, 12, 17,  2, 21],
         [ 6,  7, 29, 26, 22, 24,  5, 25],
         [22, 14, 12, 27, 27,  9, 27, 17]],
        dtype=np.int16
    )

    return motif_id_rankings_per_species


@pytest.fixture
def motif_id_rank_ratios_per_species(motif_id_rankings_per_species):
    # Create motif_id_rank_ratios_per_species by adding 1 to each ranking and dividing the
    # motif_id_rankings_per_species by the number of regions/genes.
    motif_id_rank_ratios_per_species = (
            (motif_id_rankings_per_species.astype(np.float64) + 1) / motif_id_rankings_per_species.shape[0]
    )

    return motif_id_rank_ratios_per_species


def test_calculate_cross_species_rank_ratio_with_order_statistics(motif_id_rank_ratios_per_species):
    """
    Check if _calculate_cross_species_rank_ratio_with_order_statistics returns the correct cross species combined rank
    ratio value by comparing the output with the output of OrderStatistcs.jar
    (https://github.com/aertslab/orderstatistics).
    """

    # Check if an empty array returns 1.0.
    assert _calculate_cross_species_rank_ratio_with_order_statistics(np.array([])) == 1.0

    # Check if the rank ratios for the first region/gene match.
    assert np.all(motif_id_rank_ratios_per_species[0, :] == np.array(
        [0.1333333333333333333,
         0.6,
         0.7,
         0.2333333333333333333,
         0.466666666666666666,
         0.0333333333333333333,
         0.4,
         0.766666666666666666])
    )

    # Use cross species combined rank ratios produced by OrderStatistcs.jar:
    #   https://github.com/aertslab/orderstatistics
    #
    #   - Create rank ratios TSV file from motif_id_rank_ratios_per_species:
    #       motif_id_rank_ratios_per_species_df = pd.DataFrame(
    #           motif_id_rank_ratios_per_species,
    #           index=['reg' + str(i) for i in range(1, motif_id_rank_ratios_per_species.shape{0})]
    #       )
    #
    #       motif_id_rank_ratios_per_species_df.to_csv('motif_id_rank_ratios_per_species.rr', sep='\t', header=False)
    #
    #   - Run OrderStatistcs.jar:
    #
    #       java -jar OrderStatistics.jar motif_id_rank_ratios_per_species.rr \
    #         | sort -k1,1V \
    #         > motif_id_rank_ratios_per_species.r
    #
    cross_species_rank_ratios = np.array(
        [1.471485534979422500e-02,
         5.167835571406806000e-02,
         6.632278021490627000e-02,
         1.285181832037805300e-01,
         1.904465591678098000e-02,
         4.238827173662555500e-01,
         5.788900044048133000e-02,
         1.840457054183819000e-01,
         2.213220727572009000e-01,
         4.062766817558305600e-03,
         2.229295056927296200e-01,
         2.081463199497027600e-01,
         5.077855507392157000e-02,
         3.619868131534811000e-02,
         1.254108782304526000e-01,
         5.421597502057618000e-02,
         8.512892304526771000e-03,
         5.102039135497620600e-02,
         2.598971437280905000e-02,
         1.999737075201946700e-01,
         1.811859236092062000e-02,
         4.093855263526901400e-02,
         6.014680050297216000e-03,
         4.193637887669566400e-02,
         6.253578472793826000e-01,
         6.760248552354815000e-02,
         1.274203830102119400e-01,
         2.300545926992833000e-02,
         3.078276922908097000e-01,
         4.727090892242037000e-01],
        dtype=np.float64
    )

    # Create motif_id_rank_ratios_for_one_region_or_gene numpy array which will be reused in each loop as we have to
    # copy the original data to this array as _calculate_cross_species_rank_ratio_with_order_statistics will modify this
    # array inplace.
    motif_id_rank_ratios_for_one_region_or_gene = np.zeros(
        (motif_id_rank_ratios_per_species.shape[1],),
        dtype=np.float64
    )

    for region_or_gene_id_idx in range(motif_id_rank_ratios_per_species.shape[0]):
        # Get rank ratios for the current gene region/gene.
        motif_id_rank_ratios_for_one_region_or_gene[:] = motif_id_rank_ratios_per_species[region_or_gene_id_idx, :]

        # Check if calculated cross species combined rank ratio matches with the OrderStatistics.jar version.
        assert _calculate_cross_species_rank_ratio_with_order_statistics(
            motif_id_rank_ratios_for_one_region_or_gene) == cross_species_rank_ratios[region_or_gene_id_idx]


def test_create_cross_species_ranking_for_motif(motif_id_rankings_per_species):
    """
    Test creating of cross species ranking for the input motif from the rankings for that motif for each region/gene
    for each species.
    """

    # Use cross species combined rank ratios produced by OrderStatistcs.jar
    # (https://github.com/aertslab/orderstatistics) and sort the values from low to high and assign a zero-based rank:
    #
    #  - First sort value form low to high, add zero-based rank (by appending the line number minus 1) and resort by
    #    region name, so everything is back in the original order. Column 3 contains the ranking.
    #
    #      sort -k 2,2g rank_ratios_orderstatistics.r | awk -F '\t' '{print $0 "\t" NR - 1}' | sort -k 1,1V
    #
    combined_rankings_per_region_or_gene = np.array(
        [3,
         13,
         16,
         20,
         5,
         27,
         15,
         21,
         24,
         0,
         25,
         23,
         11,
         8,
         18,
         14,
         2,
         12,
         7,
         22,
         4,
         9,
         1,
         10,
         29,
         17,
         19,
         6,
         26,
         28],
        dtype=np.int16
    )

    assert np.all(create_cross_species_ranking_for_motif(motif_id_rankings_per_species)
                  == combined_rankings_per_region_or_gene)
