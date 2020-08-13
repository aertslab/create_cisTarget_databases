import numpy as np
import pytest

from cistarget_db import FeaturesType, MotifsOrTracksType, FeatureIDs, MotifsOrTracksIDs, DatabaseTypes


def test_FeaturesType_from_str():
    """Check if a member of FeaturesType Enum can be made from a string."""
    assert FeaturesType.from_str('regions') == FeaturesType.REGIONS
    assert FeaturesType.from_str('REGIONS') == FeaturesType.REGIONS
    assert FeaturesType.from_str('genes') == FeaturesType.GENES
    assert FeaturesType.from_str('GENES') == FeaturesType.GENES

    with pytest.raises(ValueError, match=r'Unsupported FeaturesType "NON_EXISTING_FEATURE_TYPE".'):
        FeaturesType.from_str('NON_EXISTING_FEATURE_TYPE')


def test_MotifsOrTracksType_from_str():
    """Check if a member of MotifsOrTracksType Enum can be made from a string."""
    assert MotifsOrTracksType.from_str('motifs') == MotifsOrTracksType.MOTIFS
    assert MotifsOrTracksType.from_str('MOTIFS') == MotifsOrTracksType.MOTIFS
    assert MotifsOrTracksType.from_str('tracks') == MotifsOrTracksType.TRACKS
    assert MotifsOrTracksType.from_str('TRACKS') == MotifsOrTracksType.TRACKS

    with pytest.raises(ValueError, match=r'Unsupported MotifsOrTracksType "NON_EXISTING_MOTIFS_OR_TRACKS_TYPE".'):
        MotifsOrTracksType.from_str('NON_EXISTING_MOTIFS_OR_TRACKS_TYPE')


def test_FeatureIDs_with_regions():
    """Check if a FeatureIDs object can be constructed from a list of region IDs."""
    features_ids_instance = FeatureIDs(
        feature_ids=['reg2', 'reg1', 'reg6', 'reg2'], features_type=FeaturesType.REGIONS
    )
    assert features_ids_instance.features_type == FeaturesType.REGIONS
    assert features_ids_instance.feature_ids == ('reg1', 'reg2', 'reg6')

    assert eval(features_ids_instance.__repr__()) == features_ids_instance


def test_FeatureIDs_with_genes():
    """Check if a FeatureIDs object can be constructed from a set of gene IDs."""
    features_ids_instance = FeatureIDs(
        feature_ids={'gene2', 'gene1', 'gene6', 'gene2'}, features_type=FeaturesType.GENES
    )
    assert features_ids_instance.features_type == FeaturesType.GENES
    assert features_ids_instance.feature_ids == ('gene1', 'gene2', 'gene6')

    assert eval(features_ids_instance.__repr__()) == features_ids_instance


def test_FeatureIDs_with_features_type_str():
    """
    Check if a FeatureIDs object can be constructed from a tuple of gene IDs where features_type is given as a string.
    """
    features_ids_instance = FeatureIDs(
        feature_ids=('gene2', 'gene1', 'gene6', 'gene2'), features_type='gEnES'
    )
    assert features_ids_instance.features_type == FeaturesType.GENES
    assert features_ids_instance.feature_ids == ('gene1', 'gene2', 'gene6')

    assert eval(features_ids_instance.__repr__()) == features_ids_instance


def test_MotifsOrTracksIDs_with_motifs():
    """Check if a MotifsOrTracksIDs object can be constructed from a list of motif IDs."""
    motif_or_track_ids_instance = MotifsOrTracksIDs(
        motif_or_track_ids=['motif5', 'motif10', 'motif3', 'motif10'], motifs_or_tracks_type=MotifsOrTracksType.MOTIFS
    )
    assert motif_or_track_ids_instance.motifs_or_tracks_type == MotifsOrTracksType.MOTIFS
    assert motif_or_track_ids_instance.motif_or_track_ids == ('motif10', 'motif3', 'motif5')

    assert eval(motif_or_track_ids_instance.__repr__()) == motif_or_track_ids_instance


def test_MotifsOrTracksIDs_with_tracks():
    """Check if a MotifsOrTracksIDs object can be constructed from a set of track IDs."""
    motif_or_track_ids_instance = MotifsOrTracksIDs(
        motif_or_track_ids={'track5', 'track10', 'track3', 'track10'}, motifs_or_tracks_type=MotifsOrTracksType.TRACKS
    )
    assert motif_or_track_ids_instance.motifs_or_tracks_type == MotifsOrTracksType.TRACKS
    assert motif_or_track_ids_instance.motif_or_track_ids == ('track10', 'track3', 'track5')

    assert eval(motif_or_track_ids_instance.__repr__()) == motif_or_track_ids_instance


def test_MotifsOrTracksIDs_with_motifs_or_tracks_type_str():
    """
    Check if a MotifsOrTracksIDs object can be constructed from a tuple of track IDs,
    where motifs_or_tracks_type is given as a string.
    """
    motif_or_track_ids_instance = MotifsOrTracksIDs(
        motif_or_track_ids=('track5', 'track10', 'track3', 'track10'), motifs_or_tracks_type='tracks'
    )
    assert motif_or_track_ids_instance.motifs_or_tracks_type == MotifsOrTracksType.TRACKS
    assert motif_or_track_ids_instance.motif_or_track_ids == ('track10', 'track3', 'track5')

    assert eval(motif_or_track_ids_instance.__repr__()) == motif_or_track_ids_instance


def test_DatabaseTypes():
    """
    Check if all needed DatabaseTypes exist by constructing all combinations, check if the name of each member matches
    with the associated values and check if a member of DatabaseTypes Enum can be constructed from the string name.
    """
    for scores_or_rankings in ('scores', 'rankings'):
        for motif_or_tracks_type in MotifsOrTracksType.__members__:
            for feature_type in FeaturesType.__members__:
                database_type_name = f'{scores_or_rankings.upper()}_DB_{MotifsOrTracksType[motif_or_tracks_type].value.upper()}_VS_{FeaturesType[feature_type].value.upper()}'
                assert database_type_name in DatabaseTypes.__members__
                assert DatabaseTypes[database_type_name].value == (
                    scores_or_rankings,
                    MotifsOrTracksType[motif_or_tracks_type].value,
                    FeaturesType[feature_type].value
                )
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_str(database_type_name)
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_strings(
                    scores_or_rankings,
                    MotifsOrTracksType[motif_or_tracks_type].value,
                    FeaturesType[feature_type].value
                )
                del database_type_name

                database_type_name = f'{scores_or_rankings.upper()}_DB_{FeaturesType[feature_type].value.upper()}_VS_{MotifsOrTracksType[motif_or_tracks_type].value.upper()}'
                assert database_type_name in DatabaseTypes.__members__
                assert DatabaseTypes[database_type_name].value == (
                    scores_or_rankings,
                    FeaturesType[feature_type].value,
                    MotifsOrTracksType[motif_or_tracks_type].value
                )
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_str(database_type_name)
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_strings(
                    scores_or_rankings,
                    FeaturesType[feature_type].value,
                    MotifsOrTracksType[motif_or_tracks_type].value
                )
                del database_type_name

    with pytest.raises(ValueError, match=r'Unsupported DatabaseTypes "NON_EXISTING_DB_TYPE".'):
        DatabaseTypes.from_str('NON_EXISTING_DB_TYPE')

    with pytest.raises(ValueError, match=r'''"\('scores', 'motifs', 'unsupported'\)" could not be converted to a valid DatabaseTypes member.'''):
        DatabaseTypes.from_strings('scores', 'motifs', 'unsupported')


def test_DatabaseTypes_create_database_type_and_db_prefix_from_db_filename():
    """
    Check if a database filename (which includes the type of database in the name) can be converted to a member of
    DatabaseTypes Enum and a database prefix.
    """
    assert DatabaseTypes.create_database_type_and_db_prefix_and_extension_from_db_filename(
        db_filename='/some/path/test_db.tracks_vs_genes.scores.feather'
    ) == (DatabaseTypes.SCORES_DB_TRACKS_VS_GENES, '/some/path/test_db', 'feather')

    assert DatabaseTypes.create_database_type_and_db_prefix_and_extension_from_db_filename(
        db_filename='/some/path/test_db.with.extra.dots.tracks_vs_genes.scores.feather'
    ) == (DatabaseTypes.SCORES_DB_TRACKS_VS_GENES, '/some/path/test_db.with.extra.dots', 'feather')

    with pytest.raises(
            ValueError,
            match=r'Database filename "/some/path/test_db.tracks_vs_genes_scores.feather" does not contain 3 dots.'
    ):
        DatabaseTypes.create_database_type_and_db_prefix_and_extension_from_db_filename(
            db_filename='/some/path/test_db.tracks_vs_genes_scores.feather'
        )

    with pytest.raises(
            ValueError,
            match=r'Database filename "/some/path/test_db.tracks_versus_genes.scores.feather" does not contain "_vs_" in "tracks_versus_genes" part.'
    ):
        DatabaseTypes.create_database_type_and_db_prefix_and_extension_from_db_filename(
            db_filename='/some/path/test_db.tracks_versus_genes.scores.feather'
        )


def test_DatabaseTypes_create_db_filename():
    """
    Check if a database filename (which includes the type of database in the name) can be constructed from a member of
    DatabaseTypes Enum by providing a database prefix and extension.
    """
    assert DatabaseTypes.SCORES_DB_TRACKS_VS_GENES.create_db_filename(
        db_prefix='/some/path/test_db',
        extension='feather'
    ) == '/some/path/test_db.tracks_vs_genes.scores.feather'

    assert DatabaseTypes.RANKINGS_DB_REGIONS_VS_MOTIFS.create_db_filename(
        db_prefix='/some/path/test_db',
        extension='feather'
    ) == '/some/path/test_db.regions_vs_motifs.rankings.feather'


def test_DatabaseTypes_properties_and_get_dtype():
    """
    Check properties of member of DatabaseType Enum.
    """
    scores_db_tracks_vs_genes = DatabaseTypes.SCORES_DB_TRACKS_VS_GENES
    assert scores_db_tracks_vs_genes.is_scores_db is True
    assert scores_db_tracks_vs_genes.is_rankings_db is False
    assert scores_db_tracks_vs_genes.is_regions_db is False
    assert scores_db_tracks_vs_genes.is_genes_db is True
    assert scores_db_tracks_vs_genes.is_motifs_db is False
    assert scores_db_tracks_vs_genes.is_tracks_db is True
    assert scores_db_tracks_vs_genes.scores_or_rankings == 'scores'
    assert scores_db_tracks_vs_genes.features_type == FeaturesType.GENES
    assert scores_db_tracks_vs_genes.motifs_or_tracks_type == MotifsOrTracksType.TRACKS
    assert scores_db_tracks_vs_genes.column_kind == 'tracks'
    assert scores_db_tracks_vs_genes.row_kind == 'genes'

    # Score databases always store the data as 32-bit floats.
    assert scores_db_tracks_vs_genes.get_dtype(nbr_rows=20000) == np.float32
    assert scores_db_tracks_vs_genes.get_dtype(nbr_rows=32766) == np.float32
    assert scores_db_tracks_vs_genes.get_dtype(nbr_rows=32767) == np.float32
    assert scores_db_tracks_vs_genes.get_dtype(nbr_rows=32768) == np.float32
    assert scores_db_tracks_vs_genes.get_dtype(nbr_rows=32769) == np.float32
    assert scores_db_tracks_vs_genes.get_dtype(nbr_rows=1000000) == np.float32

    del scores_db_tracks_vs_genes

    rankings_db_region_vs_motifs = DatabaseTypes.RANKINGS_DB_REGIONS_VS_MOTIFS
    assert rankings_db_region_vs_motifs.is_scores_db is False
    assert rankings_db_region_vs_motifs.is_rankings_db is True
    assert rankings_db_region_vs_motifs.is_regions_db is True
    assert rankings_db_region_vs_motifs.is_genes_db is False
    assert rankings_db_region_vs_motifs.is_motifs_db is True
    assert rankings_db_region_vs_motifs.is_tracks_db is False
    assert rankings_db_region_vs_motifs.scores_or_rankings == 'rankings'
    assert rankings_db_region_vs_motifs.features_type == FeaturesType.REGIONS
    assert rankings_db_region_vs_motifs.motifs_or_tracks_type == MotifsOrTracksType.MOTIFS
    assert rankings_db_region_vs_motifs.column_kind == 'regions'
    assert rankings_db_region_vs_motifs.row_kind == 'motifs'

    # Rankings databases store the zero-based rankings as optimally as possible in a:
    #   - 16-bit signed integer: max value = 2^15 - 1 = 32767 ==> can store 32768 rankings.
    #   - 32-bit signed integer: max value = 2^31 - 1 = 2147483647 ==> can store 2147483648
    assert rankings_db_region_vs_motifs.get_dtype(nbr_rows=20000) == np.int16
    assert rankings_db_region_vs_motifs.get_dtype(nbr_rows=32766) == np.int16
    assert rankings_db_region_vs_motifs.get_dtype(nbr_rows=32767) == np.int16
    assert rankings_db_region_vs_motifs.get_dtype(nbr_rows=32768) == np.int16
    assert rankings_db_region_vs_motifs.get_dtype(nbr_rows=32769) == np.int32
    assert rankings_db_region_vs_motifs.get_dtype(nbr_rows=1000000) == np.int32

    del rankings_db_region_vs_motifs
