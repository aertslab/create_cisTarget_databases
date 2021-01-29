import numba
import numpy as np


@numba.jit(nopython=True)
def _calculate_cross_species_rank_ratio_with_order_statistics(motif_id_rank_ratios_for_one_region_or_gene: np.ndarray) -> np.ndarray:
    """
    Calculate cross-species combined rank ratio for a region/gene from rank ratios of a certain region/gene scored for
    a certain motif in multiple species with order statistics.

    Code based on applyOrderStatistics function:
      https://github.com/aertslab/orderstatistics/blob/master/OrderStatistics.java

    Paper:
      https://www.nature.com/articles/nbt1203

    :param motif_id_rank_ratios_for_one_region_or_gene:
        Numpy array of rank ratios of a certain region/gene scored for a certain motif in multiple species.
        This array is sorted inplace, so if the original array is required afterwards, provide a copy to this function.
    :return: Cross species combined rank ratio.
    """

    # Number of species for which to calculate a cross-species combined rank ratio score.
    rank_ratios_size = motif_id_rank_ratios_for_one_region_or_gene.shape[0]

    if rank_ratios_size == 0:
        return np.float64(1.0)
    else:
        # Sort rank ratios inplace.
        motif_id_rank_ratios_for_one_region_or_gene.sort()

        w = np.zeros((rank_ratios_size + 1,), dtype=np.float64)
        w[0] = np.float64(1.0)
        w[1] = motif_id_rank_ratios_for_one_region_or_gene[rank_ratios_size - 1]

        for k in range(2, rank_ratios_size + 1):
            f = np.float64(-1.0)
            for j in range(0, k):
                f = -(f * (k - j) * motif_id_rank_ratios_for_one_region_or_gene[rank_ratios_size - k]) / (j + 1.0)
                w[k] = w[k] + (w[k - j - 1] * f)

        # Cross species combined rank ratio.
        return w[rank_ratios_size]


@numba.jit(nopython=True)
def create_cross_species_ranking_for_motif(motif_id_rankings_per_species: np.ndarray) -> np.ndarray:
    """
    Create cross-species ranking for the input motif from the rankings for that motif for each region/gene for each
    species.
    
    :param motif_id_rankings_per_species:
        Numpy array with rankings for a certain motif, organized with regions/genes as rows and species as columns.
    :return: Cross species combined ranking for motif.
    """

    assert motif_id_rankings_per_species.ndim == 2, "motif_id_rankings_per_species needs to be a 2D array."

    # Rankings for the input motif:
    #   - per row: rankings for input motif for a certain region/gene for each species.
    #   - per column: rankings for input motif for a certain species for all regions/genes.
    nbr_regions_or_genes = motif_id_rankings_per_species.shape[0]
    nbr_species = motif_id_rankings_per_species.shape[1]

    # Get dtype of input ranking (np.int16 or np.int32), so the same type can be returned for the combined ranking.
    ranking_dtype = motif_id_rankings_per_species.dtype

    # Create rank ratios from rankings (starts at zero) by adding 1 to each ranking and by dividing by the total number
    # of regions or genes.
    motif_id_rank_ratios_per_species = (motif_id_rankings_per_species.astype(np.float64) + 1) / nbr_regions_or_genes

    # Create array to store rank ratios for one region/gene for each species.
    # This array is reused for each iteration of the loop.
    motif_id_rank_ratios_for_one_region_or_gene = np.zeros((nbr_species,), dtype=np.float64)

    # Create array to store the cross-species combined rank ratios for each region or gene for the input motif.
    cross_species_rank_ratios_per_region_or_gene = np.zeros((nbr_regions_or_genes,), dtype=np.float64)

    # Calculate combined rank ratio with order statistics for each region or gene for the input motif.
    for i in range(nbr_regions_or_genes):
        # Copy rank ratios for a certain region or gene for all species to a temporary array
        # (as calculate_combined_rank_ratio_with_orderstatistics modifies the input array when sorting it inplace).
        motif_id_rank_ratios_for_one_region_or_gene[:] = motif_id_rank_ratios_per_species[i, :]

        # Calculate cross-species combined rank ratio with order statistics for the current region or gene for the
        # input motif.
        cross_species_rank_ratios_per_region_or_gene[i] = _calculate_cross_species_rank_ratio_with_order_statistics(
            motif_id_rank_ratios_for_one_region_or_gene
        )

    # Convert cross-species combined rank ratios per region or gene to a cross-species combined ranking per region or
    # gene.
    cross_species_rankings_per_region_or_gene = cross_species_rank_ratios_per_region_or_gene.argsort().argsort().astype(
        ranking_dtype
    )

    # Cross species combined ranking for input motif.
    return cross_species_rankings_per_region_or_gene
