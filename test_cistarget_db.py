import numpy as np
import pandas as pd
import pytest

from cistarget_db import FeaturesType, MotifsOrTracksType, FeatureIDs, MotifOrTrackIDs, DatabaseTypes, CisTargetDatabase


def test_FeaturesType_from_str():
    """Check if a member of FeaturesType Enum can be made from a string."""
    assert FeaturesType.from_str('regions') == FeaturesType.REGIONS
    assert FeaturesType.from_str('REGIONS') == FeaturesType.REGIONS
    assert FeaturesType.from_str('genes') == FeaturesType.GENES
    assert FeaturesType.from_str('GENES') == FeaturesType.GENES

    with pytest.raises(ValueError, match=r'Unsupported FeaturesType "NON_EXISTING_FEATURES_TYPE".'):
        FeaturesType.from_str('NON_EXISTING_FEATURES_TYPE')


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
    assert features_ids_instance.type == FeaturesType.REGIONS
    assert features_ids_instance.ids == ('reg1', 'reg2', 'reg6')

    assert eval(features_ids_instance.__repr__()) == features_ids_instance
    assert len(features_ids_instance) == 3


def test_FeatureIDs_with_genes():
    """Check if a FeatureIDs object can be constructed from a set of gene IDs."""
    features_ids_instance = FeatureIDs(
        feature_ids={'gene2', 'gene1', 'gene6', 'gene2'}, features_type=FeaturesType.GENES
    )
    assert features_ids_instance.type == FeaturesType.GENES
    assert features_ids_instance.ids == ('gene1', 'gene2', 'gene6')

    assert eval(features_ids_instance.__repr__()) == features_ids_instance
    assert len(features_ids_instance) == 3


def test_FeatureIDs_with_features_type_str():
    """
    Check if a FeatureIDs object can be constructed from a tuple of gene IDs where features_type is given as a string.
    """
    features_ids_instance = FeatureIDs(
        feature_ids=('gene2', 'gene1', 'gene6', 'gene2'), features_type='gEnES'
    )
    assert features_ids_instance.type == FeaturesType.GENES
    assert features_ids_instance.ids == ('gene1', 'gene2', 'gene6')

    assert eval(features_ids_instance.__repr__()) == features_ids_instance
    assert len(features_ids_instance) == 3


def test_FeatureIDs_subset_superset():
    """Check if feature IDs of a FeatureIDs object are a subset or a superset of another FeatureIDs object."""
    features_ids_instance1 = FeatureIDs(
        feature_ids=['reg1', 'reg2', 'reg6'], features_type=FeaturesType.REGIONS
    )
    features_ids_instance2 = FeatureIDs(
        feature_ids=['reg1', 'reg2', 'reg4', 'reg6'], features_type=FeaturesType.REGIONS
    )

    assert features_ids_instance1.issubset(features_ids_instance2)
    assert features_ids_instance2.issuperset(features_ids_instance1)


def test_MotifsOrTracksIDs_with_motifs():
    """Check if a MotifOrTrackIDs object can be constructed from a list of motif IDs."""
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids=['motif5', 'motif10', 'motif3', 'motif10'], motifs_or_tracks_type=MotifsOrTracksType.MOTIFS
    )
    assert motif_or_track_ids_instance.type == MotifsOrTracksType.MOTIFS
    assert motif_or_track_ids_instance.ids == ('motif10', 'motif3', 'motif5')

    assert eval(motif_or_track_ids_instance.__repr__()) == motif_or_track_ids_instance
    assert len(motif_or_track_ids_instance) == 3


def test_MotifsOrTracksIDs_with_tracks():
    """Check if a MotifOrTrackIDs object can be constructed from a set of track IDs."""
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids={'track5', 'track10', 'track3', 'track10'}, motifs_or_tracks_type=MotifsOrTracksType.TRACKS
    )
    assert motif_or_track_ids_instance.type == MotifsOrTracksType.TRACKS
    assert motif_or_track_ids_instance.ids == ('track10', 'track3', 'track5')

    assert eval(motif_or_track_ids_instance.__repr__()) == motif_or_track_ids_instance
    assert len(motif_or_track_ids_instance) == 3


def test_MotifsOrTracksIDs_with_motifs_or_tracks_type_str():
    """
    Check if a MotifOrTrackIDs object can be constructed from a tuple of track IDs,
    where motifs_or_tracks_type is given as a string.
    """
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids=('track5', 'track10', 'track3', 'track10'), motifs_or_tracks_type='tracks'
    )
    assert motif_or_track_ids_instance.type == MotifsOrTracksType.TRACKS
    assert motif_or_track_ids_instance.ids == ('track10', 'track3', 'track5')

    assert eval(motif_or_track_ids_instance.__repr__()) == motif_or_track_ids_instance
    assert len(motif_or_track_ids_instance) == 3


def test_DatabaseTypes():
    """
    Check if all needed DatabaseTypes exist by constructing all combinations, check if the name of each member matches
    with the associated values and check if a member of DatabaseTypes Enum can be constructed from the string name.
    """
    for scores_or_rankings in ('scores', 'rankings'):
        for motif_or_tracks_type in MotifsOrTracksType.__members__:
            for features_type in FeaturesType.__members__:
                database_type_name = f'{scores_or_rankings.upper()}_DB_{MotifsOrTracksType[motif_or_tracks_type].value.upper()}_VS_{FeaturesType[features_type].value.upper()}'
                assert database_type_name in DatabaseTypes.__members__
                assert DatabaseTypes[database_type_name].value == (
                    scores_or_rankings,
                    MotifsOrTracksType[motif_or_tracks_type].value,
                    FeaturesType[features_type].value
                )
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_str(database_type_name)
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_str(
                    f'DatabaseTypes.{database_type_name}'
                )
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_strings(
                    scores_or_rankings,
                    MotifsOrTracksType[motif_or_tracks_type].value,
                    FeaturesType[features_type].value
                )
                del database_type_name

                database_type_name = f'{scores_or_rankings.upper()}_DB_{FeaturesType[features_type].value.upper()}_VS_{MotifsOrTracksType[motif_or_tracks_type].value.upper()}'
                assert database_type_name in DatabaseTypes.__members__
                assert DatabaseTypes[database_type_name].value == (
                    scores_or_rankings,
                    FeaturesType[features_type].value,
                    MotifsOrTracksType[motif_or_tracks_type].value
                )
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_str(database_type_name)
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_str(
                    f'DatabaseTypes.{database_type_name}'
                )
                assert DatabaseTypes[database_type_name] == DatabaseTypes.from_strings(
                    scores_or_rankings,
                    FeaturesType[features_type].value,
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


@pytest.fixture
def db_numpy_array_scores_db_motifs_vs_regions():
    db_numpy_array_scores_db_motifs_vs_regions = np.array(
        [[1.2, 3.0, 0.3, 5.6],
         [6.7, 3.0, 4.3, 5.6],
         [3.5, 3.0, 0.0, 0.0],
         [0.0, 3.0, 0.0, 5.6],
         [2.4, 3.0, 7.8, 1.2],
         [2.4, 3.0, 0.6, 0.0],
         [2.4, 3.0, 7.7, 0.0]],
        dtype=np.float32,
        order='C'
    )

    return db_numpy_array_scores_db_motifs_vs_regions


@pytest.fixture
def db_numpy_array_rankings_db_genes_vs_tracks():
    # Create numpy array with values which will be written to the cisTarget database dataframe.
    db_numpy_array_rankings_db_genes_vs_tracks = np.array(
        [[0, 1, 2, 3, 4, 5, 6],
         [2, 4, 3, 0, 1, 6, 5],
         [0, 6, 2, 4, 1, 3, 5],
         [2, 0, 4, 6, 5, 3, 1]],
        dtype=np.int16,
        order='F'
    )

    return db_numpy_array_rankings_db_genes_vs_tracks


def test_cistargetdatabase_basic(db_numpy_array_scores_db_motifs_vs_regions,
                                 db_numpy_array_rankings_db_genes_vs_tracks):
    # Test cisTarget SCORES_DB_MOTIFS_VS_REGIONS databases.

    features_ids_instance = FeatureIDs(
        feature_ids=['reg1', 'reg2', 'reg3', 'reg4', 'reg5', 'reg6', 'reg7'], features_type=FeaturesType.REGIONS
    )
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids=['motif1', 'motif2', 'motif3', 'motif4'], motifs_or_tracks_type=MotifsOrTracksType.MOTIFS
    )

    def check_ct_scores_db_motifs_vs_regions(ct_scores_db_motifs_vs_regions, db_numpy_array, order):
        # Check if creation of cisTarget SCORES_DB_MOTIFS_VS_REGIONS database succeeded (float32 datatype).
        assert np.all(ct_scores_db_motifs_vs_regions.df.to_numpy() == db_numpy_array)
        assert ct_scores_db_motifs_vs_regions.shape == (7, 4)
        assert ct_scores_db_motifs_vs_regions.nbr_rows == 7
        assert ct_scores_db_motifs_vs_regions.nbr_columns == 4
        assert ct_scores_db_motifs_vs_regions.dtype == np.float32
        if order == 'C':
            assert ct_scores_db_motifs_vs_regions.df.to_numpy().flags.c_contiguous is True
        elif order == 'F':
            assert ct_scores_db_motifs_vs_regions.df.to_numpy().flags.f_contiguous is True
        # Check if feature IDs and motif and track IDs are properly set.
        assert ct_scores_db_motifs_vs_regions.feature_ids == features_ids_instance
        assert ct_scores_db_motifs_vs_regions.motif_or_track_ids == motif_or_track_ids_instance
        # Columns contain motifs.
        assert ct_scores_db_motifs_vs_regions.df.columns.to_list() == list(motif_or_track_ids_instance.ids)
        # Rows contain regions.
        assert ct_scores_db_motifs_vs_regions.df.index.to_list() == list(features_ids_instance.ids)

    # Create zeroed cisTarget SCORES_DB_MOTIFS_VS_REGIONS database in C order.
    check_ct_scores_db_motifs_vs_regions(
        ct_scores_db_motifs_vs_regions=CisTargetDatabase.create_db(
            db_type=DatabaseTypes.SCORES_DB_MOTIFS_VS_REGIONS,
            feature_ids=features_ids_instance,
            motif_or_track_ids=motif_or_track_ids_instance
        ),
        db_numpy_array=np.zeros((7, 4), dtype=np.float32),
        order='C'
    )

    # Create zeroed cisTarget SCORES_DB_MOTIFS_VS_REGIONS database in Fortran order.
    check_ct_scores_db_motifs_vs_regions(
        ct_scores_db_motifs_vs_regions=CisTargetDatabase.create_db(
            db_type=DatabaseTypes.SCORES_DB_MOTIFS_VS_REGIONS,
            feature_ids=features_ids_instance,
            motif_or_track_ids=motif_or_track_ids_instance,
            order='F'
        ),
        db_numpy_array=np.zeros((7, 4), dtype=np.float32),
        order='F'
    )

    # Create cisTarget SCORES_DB_MOTIFS_VS_REGIONS database from numpy array
    # (db_numpy_array_scores_db_motifs_vs_regions is in C order).
    check_ct_scores_db_motifs_vs_regions(
        ct_scores_db_motifs_vs_regions=CisTargetDatabase.create_db(
            db_type=DatabaseTypes.SCORES_DB_MOTIFS_VS_REGIONS,
            feature_ids=features_ids_instance,
            motif_or_track_ids=motif_or_track_ids_instance,
            db_numpy_array=db_numpy_array_scores_db_motifs_vs_regions
        ),
        db_numpy_array=db_numpy_array_scores_db_motifs_vs_regions,
        order='C'
    )

    # Delete some objects so we don't accidentally reuse them in the next section.
    del features_ids_instance
    del motif_or_track_ids_instance

    # Test cisTarget RANKINGS_DB_GENES_VS_TRACKS databases.

    features_ids_instance = FeatureIDs(
        feature_ids=['gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6', 'gene7'], features_type=FeaturesType.GENES
    )
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids=['track1', 'track2', 'track3', 'track4'], motifs_or_tracks_type=MotifsOrTracksType.TRACKS
    )

    def check_ct_rankings_db_genes_vs_tracks(ct_rankings_db_genes_vs_tracks, db_numpy_array, order):
        # Check if creation of zeroed cisTarget RANKINGS_DB_GENES_VS_TRACKS database succeeded (int16 datatype).
        assert np.all(ct_rankings_db_genes_vs_tracks.df.to_numpy() == db_numpy_array)
        assert ct_rankings_db_genes_vs_tracks.shape == (4, 7)
        assert ct_rankings_db_genes_vs_tracks.nbr_rows == 4
        assert ct_rankings_db_genes_vs_tracks.nbr_columns == 7
        assert ct_rankings_db_genes_vs_tracks.dtype == np.int16
        if order == 'C':
            assert ct_rankings_db_genes_vs_tracks.df.to_numpy().flags.c_contiguous is True
        elif order == 'F':
            assert ct_rankings_db_genes_vs_tracks.df.to_numpy().flags.f_contiguous is True
        # Check if feature IDs and motif and track IDs are properly set.
        assert ct_rankings_db_genes_vs_tracks.feature_ids == features_ids_instance
        assert ct_rankings_db_genes_vs_tracks.motif_or_track_ids == motif_or_track_ids_instance
        # Columns contain genes.
        assert ct_rankings_db_genes_vs_tracks.df.columns.to_list() == list(features_ids_instance.ids)
        # Rows contain tracks.
        assert ct_rankings_db_genes_vs_tracks.df.index.to_list() == list(motif_or_track_ids_instance.ids)

    # Create zeroed cisTarget RANKINGS_DB_GENES_VS_TRACKS database in C order.
    check_ct_rankings_db_genes_vs_tracks(
        ct_rankings_db_genes_vs_tracks=CisTargetDatabase.create_db(
            db_type=DatabaseTypes.RANKINGS_DB_GENES_VS_TRACKS,
            feature_ids=features_ids_instance,
            motif_or_track_ids=motif_or_track_ids_instance
        ),
        db_numpy_array=np.zeros((4, 7), dtype=np.int16),
        order='C'
    )

    # Create zeroed cisTarget RANKINGS_DB_GENES_VS_TRACKS database in Fortran order.
    check_ct_rankings_db_genes_vs_tracks(
        ct_rankings_db_genes_vs_tracks=CisTargetDatabase.create_db(
            db_type=DatabaseTypes.RANKINGS_DB_GENES_VS_TRACKS,
            feature_ids=features_ids_instance,
            motif_or_track_ids=motif_or_track_ids_instance,
            order='F'
        ),
        db_numpy_array=np.zeros((4, 7), dtype=np.int16),
        order='F'
    )

    # Create cisTarget RANKINGS_DB_GENES_VS_TRACKS database from numpy array
    # (db_numpy_array_rankings_db_genes_vs_tracks is in Fortran order).
    check_ct_rankings_db_genes_vs_tracks(
        ct_rankings_db_genes_vs_tracks=CisTargetDatabase.create_db(
            db_type=DatabaseTypes.RANKINGS_DB_GENES_VS_TRACKS,
            feature_ids=features_ids_instance,
            motif_or_track_ids=motif_or_track_ids_instance,
            db_numpy_array=db_numpy_array_rankings_db_genes_vs_tracks
        ),
        db_numpy_array=db_numpy_array_rankings_db_genes_vs_tracks,
        order='F'
    )

    # Delete some objects so we don't accidentally reuse them in the next section.
    del features_ids_instance
    del motif_or_track_ids_instance


@pytest.fixture
def ct_scores_db_motifs_vs_regions(db_numpy_array_scores_db_motifs_vs_regions):
    # Create cisTarget SCORES_DB_MOTIFS_VS_REGIONS database.

    # Create zeroed cisTarget SCORES_DB_MOTIFS_VS_REGIONS database.
    features_ids_instance = FeatureIDs(
        feature_ids=['reg1', 'reg2', 'reg3', 'reg4', 'reg5', 'reg6', 'reg7'], features_type=FeaturesType.REGIONS
    )
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids=['motif1', 'motif2', 'motif3', 'motif4'], motifs_or_tracks_type=MotifsOrTracksType.MOTIFS
    )

    ct_scores_db_motifs_vs_regions = CisTargetDatabase.create_db(
        db_type=DatabaseTypes.SCORES_DB_MOTIFS_VS_REGIONS,
        feature_ids=features_ids_instance,
        motif_or_track_ids=motif_or_track_ids_instance,
        db_numpy_array=db_numpy_array_scores_db_motifs_vs_regions
    )

    return ct_scores_db_motifs_vs_regions


@pytest.fixture
def ct_rankings_db_genes_vs_tracks(db_numpy_array_rankings_db_genes_vs_tracks):
    # Create cisTarget RANKINGS_DB_GENES_VS_TRACKS database.

    # Create zeroed cisTarget RANKINGS_DB_GENES_VS_TRACKS database.
    features_ids_instance = FeatureIDs(
        feature_ids=['gene1', 'gene2', 'gene3', 'gene4', 'gene5', 'gene6', 'gene7'], features_type=FeaturesType.GENES
    )
    motif_or_track_ids_instance = MotifOrTrackIDs(
        motif_or_track_ids=['track1', 'track2', 'track3', 'track4'], motifs_or_tracks_type=MotifsOrTracksType.TRACKS
    )

    ct_rankings_db_genes_vs_tracks = CisTargetDatabase.create_db(
        db_type=DatabaseTypes.RANKINGS_DB_GENES_VS_TRACKS,
        feature_ids=features_ids_instance,
        motif_or_track_ids=motif_or_track_ids_instance,
        db_numpy_array=db_numpy_array_rankings_db_genes_vs_tracks
    )

    return ct_rankings_db_genes_vs_tracks


def test_cistargetdatabase_read_db_and_write_db(ct_scores_db_motifs_vs_regions, ct_rankings_db_genes_vs_tracks):
    # Generate cisTarget database name based on prefix.
    ct_scores_db_motifs_vs_regions_db_filename = ct_scores_db_motifs_vs_regions.create_db_filename_from_db_prefix(
        db_prefix='test/ct_scores_db_motifs_vs_regions'
    )

    assert ct_scores_db_motifs_vs_regions_db_filename \
           == 'test/ct_scores_db_motifs_vs_regions.motifs_vs_regions.scores.feather'

    # Write cisTarget database to Feather file.
    ct_scores_db_motifs_vs_regions_db_filename_returned = ct_scores_db_motifs_vs_regions.write_db(
        db_prefix='test/ct_scores_db_motifs_vs_regions'
    )

    # Check if the same database name is generated by create_db_filename_from_db_prefix as by write_db.
    assert ct_scores_db_motifs_vs_regions_db_filename == ct_scores_db_motifs_vs_regions_db_filename_returned

    # Read cisTarget database from Feather file.
    ct_scores_db_motifs_vs_regions_read_from_feather = CisTargetDatabase.read_db(
        db_filename=ct_scores_db_motifs_vs_regions_db_filename
    )

    # Check if the cisTarget database object read from the Feather file is the same than the one that was written
    # to the Feather file. The numpy array underlying the cisTarget database object will be in Fortran order when
    # reading it from a Feather file.
    assert ct_scores_db_motifs_vs_regions_read_from_feather.db_type == ct_scores_db_motifs_vs_regions.db_type
    assert ct_scores_db_motifs_vs_regions_read_from_feather.dtype == ct_scores_db_motifs_vs_regions.dtype
    assert ct_scores_db_motifs_vs_regions_read_from_feather.shape == ct_scores_db_motifs_vs_regions.shape
    assert ct_scores_db_motifs_vs_regions_read_from_feather.feature_ids == ct_scores_db_motifs_vs_regions.feature_ids
    assert ct_scores_db_motifs_vs_regions_read_from_feather.motif_or_track_ids == ct_scores_db_motifs_vs_regions.motif_or_track_ids
    assert np.all(ct_scores_db_motifs_vs_regions_read_from_feather.df == ct_scores_db_motifs_vs_regions.df)
    assert ct_scores_db_motifs_vs_regions_read_from_feather.df.to_numpy().flags.f_contiguous is True

    # Delete some objects so we don't accidentally reuse them in the next section.
    del ct_scores_db_motifs_vs_regions
    del ct_scores_db_motifs_vs_regions_db_filename
    del ct_scores_db_motifs_vs_regions_db_filename_returned
    del ct_scores_db_motifs_vs_regions_read_from_feather

    # Generate cisTarget database name based on prefix.
    ct_rankings_db_genes_vs_tracks_db_filename = ct_rankings_db_genes_vs_tracks.create_db_filename_from_db_prefix(
        db_prefix='test/ct_rankings_db_genes_vs_tracks'
    )

    assert ct_rankings_db_genes_vs_tracks_db_filename \
           == 'test/ct_rankings_db_genes_vs_tracks.genes_vs_tracks.rankings.feather'

    # Write cisTarget database to Feather file.
    ct_rankings_db_genes_vs_tracks_db_filename_returned = ct_rankings_db_genes_vs_tracks.write_db(
        db_prefix='test/ct_rankings_db_genes_vs_tracks'
    )

    # Check if the same database name is generated by create_db_filename_from_db_prefix as by write_db.
    assert ct_rankings_db_genes_vs_tracks_db_filename == ct_rankings_db_genes_vs_tracks_db_filename_returned

    # Read cisTarget database from Feather file.
    ct_rankings_db_genes_vs_tracks_read_from_feather = CisTargetDatabase.read_db(
        db_filename=ct_rankings_db_genes_vs_tracks_db_filename
    )

    # Check if the cisTarget database object read from the Feather file is the same than the one that was written
    # to the Feather file. The numpy array underlying the cisTarget database object will be in Fortran order when
    # reading it from a Feather file.
    assert ct_rankings_db_genes_vs_tracks_read_from_feather.db_type == ct_rankings_db_genes_vs_tracks.db_type
    assert ct_rankings_db_genes_vs_tracks_read_from_feather.dtype == ct_rankings_db_genes_vs_tracks.dtype
    assert ct_rankings_db_genes_vs_tracks_read_from_feather.shape == ct_rankings_db_genes_vs_tracks.shape
    assert ct_rankings_db_genes_vs_tracks_read_from_feather.feature_ids == ct_rankings_db_genes_vs_tracks.feature_ids
    assert ct_rankings_db_genes_vs_tracks_read_from_feather.motif_or_track_ids == ct_rankings_db_genes_vs_tracks.motif_or_track_ids
    assert np.all(ct_rankings_db_genes_vs_tracks_read_from_feather.df == ct_rankings_db_genes_vs_tracks.df)
    assert ct_rankings_db_genes_vs_tracks_read_from_feather.df.to_numpy().flags.f_contiguous is True

    # Delete some objects so we don't accidentally reuse them in the next section.
    del ct_rankings_db_genes_vs_tracks_db_filename
    del ct_rankings_db_genes_vs_tracks_read_from_feather

    # Write cisTarget database to Feather file with a custom name.
    ct_rankings_db_genes_vs_tracks_db_filename_with_custom_name_returned = ct_rankings_db_genes_vs_tracks.write_db(
        db_filename='test/ct_rankings_db_genes_vs_tracks_with_custom_name.db'
    )

    # Check if the correct database name is used by write_db.
    assert ct_rankings_db_genes_vs_tracks_db_filename_with_custom_name_returned \
           == 'test/ct_rankings_db_genes_vs_tracks_with_custom_name.db'

    # Read cisTarget database from Feather file with a custom name (database type can not be automatically retrieved).
    ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name = CisTargetDatabase.read_db(
        db_filename=ct_rankings_db_genes_vs_tracks_db_filename_with_custom_name_returned,
        db_type=DatabaseTypes.RANKINGS_DB_GENES_VS_TRACKS
    )

    # Check if the cisTarget database object read from the Feather file is the same than the one that was written
    # to the Feather file. The numpy array underlying the cisTarget database object will be in Fortran order when
    # reading it from a Feather file.
    assert ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name.db_type == ct_rankings_db_genes_vs_tracks.db_type
    assert ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name.dtype == ct_rankings_db_genes_vs_tracks.dtype
    assert ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name.shape == ct_rankings_db_genes_vs_tracks.shape
    assert ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name.feature_ids == ct_rankings_db_genes_vs_tracks.feature_ids
    assert ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name.motif_or_track_ids == ct_rankings_db_genes_vs_tracks.motif_or_track_ids
    assert np.all(ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name.df == ct_rankings_db_genes_vs_tracks.df)
    assert ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name.df.to_numpy().flags.f_contiguous is True

    # Delete some objects so we don't accidentally reuse them in the next section.
    del ct_rankings_db_genes_vs_tracks
    del ct_rankings_db_genes_vs_tracks_db_filename_with_custom_name_returned
    del ct_rankings_db_genes_vs_tracks_read_from_feather_with_custom_name


def test_cistargetdatabase_transpose(ct_scores_db_motifs_vs_regions, ct_rankings_db_genes_vs_tracks):
    # Test creating cisTarget SCORES_DB_REGIONS_VS_MOTIFS databases by transposing the cisTarget
    # SCORES_DB_MOTIFS_VS_REGIONS database with different options.

    # Create transposed numpy array from ct_scores_db_motifs_vs_regions.
    db_numpy_array_scores_db_motifs_vs_regions_transposed = ct_scores_db_motifs_vs_regions.df.to_numpy().transpose()

    # Create a cisTarget SCORES_DB_REGIONS_VS_MOTIFS database by transposing the cisTarget SCORES_DB_MOTIFS_VS_REGIONS
    # database.
    ct_scores_db_regions_vs_motifs = ct_scores_db_motifs_vs_regions.transpose()

    assert np.all(
        ct_scores_db_regions_vs_motifs.df.to_numpy() == db_numpy_array_scores_db_motifs_vs_regions_transposed
    )
    assert ct_scores_db_regions_vs_motifs.db_type == DatabaseTypes.SCORES_DB_REGIONS_VS_MOTIFS

    del ct_scores_db_regions_vs_motifs

    # Create a cisTarget SCORES_DB_REGIONS_VS_MOTIFS database in C order by transposing the cisTarget
    # SCORES_DB_MOTIFS_VS_REGIONS database.
    ct_scores_db_regions_vs_motifs_order_C = ct_scores_db_motifs_vs_regions.transpose(order='C')

    assert np.all(
        ct_scores_db_regions_vs_motifs_order_C.df.to_numpy() == db_numpy_array_scores_db_motifs_vs_regions_transposed
    )
    assert ct_scores_db_regions_vs_motifs_order_C.db_type == DatabaseTypes.SCORES_DB_REGIONS_VS_MOTIFS
    assert ct_scores_db_regions_vs_motifs_order_C.df.to_numpy().flags.c_contiguous is True

    del ct_scores_db_regions_vs_motifs_order_C

    # Create a cisTarget SCORES_DB_REGIONS_VS_MOTIFS database in Fortran order by transposing the cisTarget
    # SCORES_DB_MOTIFS_VS_REGIONS database.
    ct_scores_db_regions_vs_motifs_order_F = ct_scores_db_motifs_vs_regions.transpose(order='F')

    assert np.all(
        ct_scores_db_regions_vs_motifs_order_F.df.to_numpy() == db_numpy_array_scores_db_motifs_vs_regions_transposed
    )
    assert ct_scores_db_regions_vs_motifs_order_F.db_type == DatabaseTypes.SCORES_DB_REGIONS_VS_MOTIFS
    assert ct_scores_db_regions_vs_motifs_order_F.df.to_numpy().flags.f_contiguous is True

    del ct_scores_db_regions_vs_motifs_order_F
    del ct_scores_db_motifs_vs_regions
    del db_numpy_array_scores_db_motifs_vs_regions_transposed

    # Test creating cisTarget RANKINGS_DB_TRACKS_VS_GENES databases by transposing the cisTarget
    # RANKINGS_DB_GENES_VS_TRACK database with different options.

    # Create transposed numpy array from ct_rankings_db_genes_vs_tracks.
    db_numpy_array_rankings_db_genes_vs_tracks_transposed = ct_rankings_db_genes_vs_tracks.df.to_numpy().transpose()

    # Create a cisTarget RANKINGS_DB_TRACKS_VS_GENES database by transposing the cisTarget RANKINGS_DB_GENES_VS_TRACKS
    # database.
    ct_rankings_db_tracks_vs_genes = ct_rankings_db_genes_vs_tracks.transpose()

    assert np.all(
        ct_rankings_db_tracks_vs_genes.df.to_numpy() == db_numpy_array_rankings_db_genes_vs_tracks_transposed
    )
    assert ct_rankings_db_tracks_vs_genes.db_type == DatabaseTypes.RANKINGS_DB_TRACKS_VS_GENES

    del ct_rankings_db_tracks_vs_genes

    # Create a cisTarget RANKINGS_DB_TRACKS_VS_GENES database in C order by transposing the cisTarget
    # RANKINGS_DB_GENES_VS_TRACKS database.
    ct_rankings_db_tracks_vs_genes_order_C = ct_rankings_db_genes_vs_tracks.transpose(order='C')

    assert np.all(
        ct_rankings_db_tracks_vs_genes_order_C.df.to_numpy() == db_numpy_array_rankings_db_genes_vs_tracks_transposed
    )
    assert ct_rankings_db_tracks_vs_genes_order_C.db_type == DatabaseTypes.RANKINGS_DB_TRACKS_VS_GENES
    assert ct_rankings_db_tracks_vs_genes_order_C.df.to_numpy().flags.c_contiguous is True

    del ct_rankings_db_tracks_vs_genes_order_C

    # Create a cisTarget RANKINGS_DB_TRACKS_VS_GENES database in C order by transposing the cisTarget
    # RANKINGS_DB_GENES_VS_TRACKS database.
    ct_rankings_db_tracks_vs_genes_order_F = ct_rankings_db_genes_vs_tracks.transpose(order='F')

    assert np.all(
        ct_rankings_db_tracks_vs_genes_order_F.df.to_numpy() == db_numpy_array_rankings_db_genes_vs_tracks_transposed
    )
    assert ct_rankings_db_tracks_vs_genes_order_F.db_type == DatabaseTypes.RANKINGS_DB_TRACKS_VS_GENES
    assert ct_rankings_db_tracks_vs_genes_order_F.df.to_numpy().flags.f_contiguous is True

    del ct_rankings_db_tracks_vs_genes_order_F
    del ct_rankings_db_genes_vs_tracks


def test_cistargetdatabase_update_scores_for_motif_or_track(ct_scores_db_motifs_vs_regions):
    # Dataframe with (new) scores for certain motif.
    df_scores_for_motif_or_track = pd.DataFrame(
        np.array(
            [[2.4, 6.7],
             [4.5, 7.3],
             [6.7, 0.2]],
            dtype=np.float32
        ),
        index=['reg2', 'reg7', 'reg5'],
        columns=['some_random_column', 'crm_score']
    )

    # Update values (from "df_scores_for_motif_or_track ==> pandas dataframe with "crm_score" column) for some regions
    # for "motif3" in cisTarget SCORES_DB_MOTIFS_VS_REGIONS database.
    ct_scores_db_motifs_vs_regions.update_scores_for_motif_or_track(
        motif_or_track_id='motif3',
        df_scores_for_motif_or_track=df_scores_for_motif_or_track
    )

    # Create numpy array with updated values (in column 3 = "motif3" related values).
    ct_scores_db_motifs_vs_regions_numpy = np.array(
        [[1.2, 3.0, 0.3, 5.6],
         [6.7, 3.0, 6.7, 5.6],
         [3.5, 3.0, 0.0, 0.0],
         [0.0, 3.0, 0.0, 5.6],
         [2.4, 3.0, 0.2, 1.2],
         [2.4, 3.0, 0.6, 0.0],
         [2.4, 3.0, 7.3, 0.0]],
        dtype=np.float32
    )

    assert np.all(ct_scores_db_motifs_vs_regions.df.to_numpy() == ct_scores_db_motifs_vs_regions_numpy)

    del ct_scores_db_motifs_vs_regions_numpy

    # Update values (from df_scores_for_motif_or_track["some_random_column"] ==> pandas series) for some regions for
    # "motif3" in cisTarget SCORES_DB_MOTIFS_VS_REGIONS database.
    ct_scores_db_motifs_vs_regions.update_scores_for_motif_or_track(
        motif_or_track_id='motif3',
        df_scores_for_motif_or_track=df_scores_for_motif_or_track["some_random_column"]
    )

    # Create numpy array with updated values (in column 3 = "motif3" related values).
    ct_scores_db_motifs_vs_regions_numpy = np.array(
        [[1.2, 3.0, 0.3, 5.6],
         [6.7, 3.0, 2.4, 5.6],
         [3.5, 3.0, 0.0, 0.0],
         [0.0, 3.0, 0.0, 5.6],
         [2.4, 3.0, 6.7, 1.2],
         [2.4, 3.0, 0.6, 0.0],
         [2.4, 3.0, 4.5, 0.0]],
        dtype=np.float32
    )

    assert np.all(ct_scores_db_motifs_vs_regions.df.to_numpy() == ct_scores_db_motifs_vs_regions_numpy)

    # Create a cisTarget SCORES_DB_REGIONS_VS_MOTIFS database by transposing the cisTarget SCORES_DB_MOTIFS_VS_REGIONS
    # database.
    ct_scores_db_regions_vs_motifs = ct_scores_db_motifs_vs_regions.transpose()

    del ct_scores_db_motifs_vs_regions
    del ct_scores_db_motifs_vs_regions_numpy

    # Update values (from first column of df_scores_for_motifs_or_track ==> pandas dataframe with 1 column) for some
    # regions for "motif2" in cisTarget SCORES_DB_REGIONS_VS_MOTIFS database.
    ct_scores_db_regions_vs_motifs.update_scores_for_motif_or_track(
        motif_or_track_id='motif2',
        df_scores_for_motif_or_track=pd.DataFrame(
            np.array(
                [[2.4],
                 [4.5],
                 [6.7]],
                dtype=np.float32
            ),
            index=['reg6', 'reg1', 'reg3'],
        )
    )

    # Create numpy array with updated values (in row 2 = "motif2" related values).
    ct_scores_db_regions_vs_motifs_numpy = np.array(
        [[1.2, 6.7, 3.5, 0.0, 2.4, 2.4, 2.4],
         [4.5, 3.0, 6.7, 3.0, 3.0, 2.4, 3.0],
         [0.3, 2.4, 0.0, 0.0, 6.7, 0.6, 4.5],
         [5.6, 5.6, 0.0, 5.6, 1.2, 0.0, 0.0]],
        dtype=np.float32
    )

    assert np.all(ct_scores_db_regions_vs_motifs.df.to_numpy() == ct_scores_db_regions_vs_motifs_numpy)

    del ct_scores_db_regions_vs_motifs
    del ct_scores_db_regions_vs_motifs_numpy


def test_cistargetdatabase_convert_scores_db_to_rankings_db(ct_scores_db_motifs_vs_regions):
    # Initialize random number generator.
    rng = np.random.default_rng(seed=123456)

    # Check the shape.
    assert ct_scores_db_motifs_vs_regions.shape == (7, 4)

    # Get 4 times a permutation (once for each motif) as those same values will be used later in
    # convert_scores_db_to_rankings_db for sorting ties.
    assert np.all(rng.permutation(np.arange(7)) == np.array([2, 0, 5, 4, 6, 1, 3]))
    assert np.all(rng.permutation(np.arange(7)) == np.array([0, 3, 5, 6, 1, 2, 4]))
    assert np.all(rng.permutation(np.arange(7)) == np.array([5, 1, 3, 0, 6, 4, 2]))
    assert np.all(rng.permutation(np.arange(7)) == np.array([5, 3, 2, 6, 1, 4, 0]))

    # Convert cisTarget SCORES_DB_MOTIFS_VS_REGIONS database to cisTarget RANKINGS_DB_MOTIFS_VS_REGIONS database.
    ct_rankings_db_motifs_vs_regions = ct_scores_db_motifs_vs_regions.convert_scores_db_to_rankings_db(seed=123456)

    assert ct_rankings_db_motifs_vs_regions.db_type == DatabaseTypes.RANKINGS_DB_MOTIFS_VS_REGIONS
    assert ct_rankings_db_motifs_vs_regions.dtype == np.int16

    ct_rankings_db_motifs_vs_regions_numpy = np.array(
        [[5, 0, 4, 2],
         [0, 4, 2, 1],
         [1, 5, 6, 5],
         [6, 1, 5, 0],
         [3, 6, 0, 3],
         [2, 2, 3, 4],
         [4, 3, 1, 6]],
        dtype=np.int16
    )

    # Check if it creates this ranking when seed is set to 123456 (to resolve ties).
    assert np.all(ct_rankings_db_motifs_vs_regions.df.to_numpy() == ct_rankings_db_motifs_vs_regions_numpy)

    # Create a cisTarget SCORES_DB_REGIONS_VS_MOTIFS database by transposing the cisTarget SCORES_DB_MOTIFS_VS_REGIONS
    # database.
    ct_scores_db_regions_vs_motifs = ct_scores_db_motifs_vs_regions.transpose()
    ct_rankings_db_regions_vs_motifs_numpy = ct_rankings_db_motifs_vs_regions_numpy.transpose()

    del ct_scores_db_motifs_vs_regions
    del ct_rankings_db_motifs_vs_regions
    del ct_rankings_db_motifs_vs_regions_numpy

    # Convert cisTarget SCORES_DB_REGIONS_VS_MOTIFS database to cisTarget RANKINGS_DB_REGIONS_VS_MOTIFS database.
    ct_rankings_db_regions_vs_motifs = ct_scores_db_regions_vs_motifs.convert_scores_db_to_rankings_db(seed=123456)

    assert ct_rankings_db_regions_vs_motifs.db_type == DatabaseTypes.RANKINGS_DB_REGIONS_VS_MOTIFS
    assert ct_rankings_db_regions_vs_motifs.dtype == np.int16

    # Check if it creates this ranking when seed is set to 123456 (to resolve ties).
    assert np.all(ct_rankings_db_regions_vs_motifs.df.to_numpy() == ct_rankings_db_regions_vs_motifs_numpy)


@pytest.fixture()
def ct_rankings_db_motifs_vs_regions(ct_scores_db_motifs_vs_regions):
    # Convert cisTarget SCORES_DB_MOTIFS_VS_REGIONS database to cisTarget RANKINGS_DB_MOTIFS_VS_REGIONS database.
    ct_rankings_db_motifs_vs_regions = ct_scores_db_motifs_vs_regions.convert_scores_db_to_rankings_db(seed=123456)

    return ct_rankings_db_motifs_vs_regions


def test_cistargetdatabase_create_cross_species_rankings_db(ct_rankings_db_motifs_vs_regions):
    # Write cisTarget RANKINGS_DB_MOTIFS_VS_REGIONS database.
    ct_rankings_db_motifs_vs_regions_db_filename = ct_rankings_db_motifs_vs_regions.write_db(
        db_prefix='test/ct_rankings_db_genes_vs_tracks'
    )

    # Create cross species CisTarget rankings database from individual rankings databases (in this case the same ones).
    ct_cross_species_rankings_db_motifs_vs_regions = CisTargetDatabase.create_cross_species_rankings_db(
        (ct_rankings_db_motifs_vs_regions_db_filename,
         ct_rankings_db_motifs_vs_regions_db_filename,
         ct_rankings_db_motifs_vs_regions_db_filename,
         ct_rankings_db_motifs_vs_regions_db_filename,
         ct_rankings_db_motifs_vs_regions_db_filename,
         ct_rankings_db_motifs_vs_regions_db_filename,
         ct_rankings_db_motifs_vs_regions_db_filename)
    )

    # Check if the cross species CisTarget rankings database made from the same 7 species rankings CisTarget databases
    # give the same database.
    assert ct_cross_species_rankings_db_motifs_vs_regions.db_type == ct_rankings_db_motifs_vs_regions.db_type
    assert np.all(ct_cross_species_rankings_db_motifs_vs_regions.df.to_numpy()
                  == ct_rankings_db_motifs_vs_regions.df.to_numpy())
    assert ct_cross_species_rankings_db_motifs_vs_regions.feature_ids == ct_rankings_db_motifs_vs_regions.feature_ids
    assert ct_cross_species_rankings_db_motifs_vs_regions.motif_or_track_ids == ct_rankings_db_motifs_vs_regions.motif_or_track_ids
