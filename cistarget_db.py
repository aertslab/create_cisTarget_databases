import numpy as np
import pandas as pd

from enum import Enum, unique
from typing import Dict, List, Set, Tuple, Type, Union

from cytoolz import memoize


@unique
class FeaturesType(Enum):
    """Enum describing all possible features types."""

    REGIONS = 'regions'
    GENES = 'genes'

    @classmethod
    def from_str(cls, feature_type: str) -> 'FeaturesType':
        """
        Create FeaturesType Enum member from string.

        :param feature_type: 'regions' or 'genes'.
        :return: FeatureType Enum member.
        """

        feature_type = feature_type.upper()
        feature_type_instance = cls.__members__.get(feature_type)
        if feature_type_instance:
            return feature_type_instance
        else:
            raise ValueError(f'Unsupported FeaturesType "{feature_type}".')


@unique
class MotifsOrTracksType(Enum):
    """Enum describing all possible motif or track types."""

    MOTIFS = 'motifs'
    TRACKS = 'tracks'

    @classmethod
    def from_str(cls, motifs_or_tracks_type: str) -> 'MotifsOrTracksType':
        """
        Create MotifsOrTracksType Enum member from string.

        :param motifs_or_tracks_type: 'motifs' or 'tracks'.
        :return: MotifsOrTracksType Enum member.
        """

        motifs_or_tracks_type = motifs_or_tracks_type.upper()
        motifs_or_tracks_type_instance = cls.__members__.get(motifs_or_tracks_type)
        if motifs_or_tracks_type_instance:
            return motifs_or_tracks_type_instance
        else:
            raise ValueError(f'Unsupported MotifsOrTracksType "{motifs_or_tracks_type}".')


class FeatureIDs:
    """
    FeatureIDs class represents a unique sorted tuple of region IDs or gene IDs for constructing a Pandas dataframe
    index for a cisTarget database.
    """

    def __init__(self,
                 feature_ids: Union[List[str], Set[str], Tuple[str, ...]],
                 features_type: Union[FeaturesType, str]):
        """
        Create unique sorted tuple of region IDs or gene IDs from a list, set or tuple of region IDs or gene IDs,
        annotated with FeaturesType Enum.

        :param feature_ids: list, set or tuple of region IDs or gene IDs.
        :param features_type: FeaturesType.REGIONS ("regions") or FeaturesType.GENES ("genes").
        """

        if isinstance(features_type, str):
            features_type = FeaturesType.from_str(features_type)

        self.feature_ids = tuple(sorted(set(feature_ids)))
        self.features_type = features_type

    def __repr__(self) -> str:
        return f'FeatureIDs(\n  feature_ids={self.feature_ids},\n  features_type={self.features_type}\n)'

    def __eq__(self, other: 'FeatureIDs'):
        return self.features_type == other.features_type and self.feature_ids == other.feature_ids

    def __len__(self) -> int:
        return len(self.feature_ids)


class MotifsOrTracksIDs:
    """
    MotifsOrTracksIDs class represents a unique sorted tuple of motif IDs or track IDs for constructing a Pandas
    dataframe index for a cisTarget database.
    """

    def __init__(self,
                 motif_or_track_ids: Union[List[str], Set[str], Tuple[str, ...]],
                 motifs_or_tracks_type: Union[MotifsOrTracksType, str]):
        """
        Create unique sorted tuple of motif IDs or track IDs from a list, set or tuple of motif IDs or track IDs,
        annotated with MotifsOrTracksType Enum.

        :param motif_or_track_ids: list, set or tuple of motif IDs or track IDs.
        :param motifs_or_tracks_type: MotifsOrTracksType.MOTIFS ("motifs") or MotifsOrTracksType.TRACKS ("tracks").
        """

        if isinstance(motifs_or_tracks_type, str):
            motifs_or_tracks_type = MotifsOrTracksType.from_str(motifs_or_tracks_type)

        self.motif_or_track_ids = tuple(sorted(set(motif_or_track_ids)))
        self.motifs_or_tracks_type = motifs_or_tracks_type

    def __repr__(self) -> str:
        return f'MotifsOrTracksIDs(\n  motif_or_track_ids={self.motif_or_track_ids},\n  motifs_or_tracks_type={self.motifs_or_tracks_type}\n)'

    def __eq__(self, other: 'MotifsOrTracksIDs'):
        return (self.motifs_or_tracks_type == other.motifs_or_tracks_type and
                self.motif_or_track_ids == other.motif_or_track_ids)

    def __len__(self) -> int:
        return len(self.motif_or_track_ids)


@unique
class DatabaseTypes(Enum):
    """
    Enum describing all cisTarget database types.
    """

    def __init__(self, scores_or_rankings: str, column_kind: str, row_kind: str):
        # Type of data contained in the cisTarget database: scores or rankings.
        self._scores_or_rankings = scores_or_rankings
        # Describe type of data in the column index of the cisTarget database.
        self._column_kind = column_kind
        # Describe type of data in the row index of the cisTarget database.
        self._row_kind = row_kind

    SCORES_DB_MOTIFS_VS_REGIONS = ('scores', 'motifs', 'regions')
    SCORES_DB_MOTIFS_VS_GENES = ('scores', 'motifs', 'genes')
    SCORES_DB_TRACKS_VS_REGIONS = ('scores', 'tracks', 'regions')
    SCORES_DB_TRACKS_VS_GENES = ('scores', 'tracks', 'genes')

    SCORES_DB_REGIONS_VS_MOTIFS = ('scores', 'regions', 'motifs')
    SCORES_DB_GENES_VS_MOTIFS = ('scores', 'genes', 'motifs')
    SCORES_DB_REGIONS_VS_TRACKS = ('scores', 'regions', 'tracks')
    SCORES_DB_GENES_VS_TRACKS = ('scores', 'genes', 'tracks')

    RANKINGS_DB_MOTIFS_VS_REGIONS = ('rankings', 'motifs', 'regions')
    RANKINGS_DB_MOTIFS_VS_GENES = ('rankings', 'motifs', 'genes')
    RANKINGS_DB_TRACKS_VS_REGIONS = ('rankings', 'tracks', 'regions')
    RANKINGS_DB_TRACKS_VS_GENES = ('rankings', 'tracks', 'genes')

    RANKINGS_DB_REGIONS_VS_MOTIFS = ('rankings', 'regions', 'motifs')
    RANKINGS_DB_GENES_VS_MOTIFS = ('rankings', 'genes', 'motifs')
    RANKINGS_DB_REGIONS_VS_TRACKS = ('rankings', 'regions', 'tracks')
    RANKINGS_DB_GENES_VS_TRACKS = ('rankings', 'genes', 'tracks')

    @classmethod
    def from_str(cls, database_type: str) -> 'DatabaseTypes':
        """
        Create DatabaseTypes Enum member from string.

        :param database_type: Type of database as string (member of DatabaseTypes Enum).
        :return: DatabaseTypes Enum member.
        """

        database_type = database_type.upper()
        if database_type.startswith('DATABASETYPES.'):
            database_type = database_type[14:]

        database_type_instance = cls.__members__.get(database_type)
        if database_type_instance:
            return database_type_instance
        else:
            raise ValueError(f'Unsupported DatabaseTypes "{database_type}".')

    @classmethod
    def from_strings(cls, scores_or_rankings: str, column_kind: str, row_kind: str) -> 'DatabaseTypes':
        """
        Create DatabaseTypes Enum member.

        :param scores_or_rankings:
            Type of data contained in the cisTarget database: "scores" or "rankings".
        :param column_kind:
            Type of data in the column index of the cisTarget database: 'regions', 'genes', 'motifs' or 'tracks'.
        :param row_kind:
            Type of data in the row index of the cisTarget database: 'regions', 'genes', 'motifs' or 'tracks'.
        :return: DatabaseTypes Enum member.
        """

        database_type_tuple = (scores_or_rankings, column_kind, row_kind)

        for database_type in DatabaseTypes.__members__.values():
            if database_type.value == database_type_tuple:
                return database_type

        raise ValueError(f'"{database_type_tuple}" could not be converted to a valid DatabaseTypes member.')

    @classmethod
    def create_database_type_and_db_prefix_and_extension_from_db_filename(cls, db_filename: str) -> Tuple['DatabaseTypes', str, str]:
        """
        Create DatabaseTypes Enum member from cisTarget database filename.

        :param db_filename:
            cisTarget database filename (in a format generated by create_db_name() on a DatabaseTypes Enum member).
        :return:
            DatabaseTypes Enum member and database prefix and extension.
        """

        db_filename_dot_splitted = db_filename.split('.')

        if len(db_filename_dot_splitted) < 4:
            raise ValueError(f'Database filename "{db_filename}" does not contain 3 dots.')

        db_prefix = '.'.join(db_filename_dot_splitted[0:-3])
        extension = db_filename_dot_splitted[-1]
        scores_or_rankings = db_filename_dot_splitted[-2]
        column_kind_vs_row_kind = db_filename_dot_splitted[-3].split('_vs_')

        if len(column_kind_vs_row_kind) != 2:
            raise ValueError(
                f'Database filename "{db_filename}" does not contain "_vs_" in "{db_filename_dot_splitted[-3]}" part.'
            )

        column_kind, row_kind = column_kind_vs_row_kind

        return cls.from_strings(scores_or_rankings, column_kind, row_kind), db_prefix, extension

    def create_db_filename(self, db_prefix: str, extension: str = 'feather') -> str:
        """
        Create cisTarget database filename based on a database prefix string and extension.
        Between the database prefix string and extension the database type will be encoded.

        :param db_prefix: Database prefix.
        :param extension: Database extension.
        :return: Database filename.
        """

        return f'{db_prefix}.{self._column_kind}_vs_{self._row_kind}.{self._scores_or_rankings}.{extension}'

    @property
    def is_scores_db(self) -> bool:
        """Check if cisTarget database contains scores."""
        return 'scores' == self._scores_or_rankings

    @property
    def is_rankings_db(self) -> bool:
        """Check if cisTarget database contains rankings."""
        return 'rankings' == self._scores_or_rankings

    def _is_some_kind_of_db_by_checking_column_and_row_kind(self, some_kind: str) -> bool:
        """Check if cisTarget database has some_kind set in column_kind or row_kind."""
        return some_kind == self._column_kind or some_kind == self._row_kind

    @property
    def is_regions_db(self) -> bool:
        """Check if cisTarget database has regions in columns or rows."""
        return self._is_some_kind_of_db_by_checking_column_and_row_kind('regions')

    @property
    def is_genes_db(self) -> bool:
        """Check if cisTarget database has genes in columns or rows."""
        return self._is_some_kind_of_db_by_checking_column_and_row_kind('genes')

    @property
    def is_motifs_db(self) -> bool:
        """Check if cisTarget database has motifs in columns or rows."""
        return self._is_some_kind_of_db_by_checking_column_and_row_kind('motifs')

    @property
    def is_tracks_db(self) -> bool:
        """Check if cisTarget database has tracks in columns or rows."""
        return self._is_some_kind_of_db_by_checking_column_and_row_kind('tracks')

    @property
    def scores_or_rankings(self) -> str:
        """Return 'scores' or 'rankings' for DatabaseTypes member."""
        return self._scores_or_rankings

    @property
    def features_type(self) -> 'FeaturesType':
        """Return FeaturesType Enum member for DatabaseTypes member."""
        if self.is_regions_db:
            return FeaturesType.REGIONS
        elif self.is_genes_db:
            return FeaturesType.GENES

    @property
    def motifs_or_tracks_type(self) -> 'MotifsOrTracksType':
        """Return MotifsOrTracksType Enum member for DatabaseTypes member."""
        if self.is_motifs_db:
            return MotifsOrTracksType.MOTIFS
        elif self.is_tracks_db:
            return MotifsOrTracksType.TRACKS

    @property
    def column_kind(self) -> str:
        """Return column kind for DatabaseTypes member."""
        return self._column_kind

    @property
    def row_kind(self) -> str:
        """Return row kind for DatabaseTypes member."""
        return self._row_kind

    def get_dtype(self, nbr_rows: int) -> Type[Union[np.core.single, np.core.short, np.core.intc]]:
        """
        Get optimal dtype for storing values in cisTarget database.

        :param nbr_rows: Number of rows in the database.
        :return: dtype most suited for DatabaseTypes member.
        """

        if self.is_scores_db:
            # The precision of a 32-bit float should be good enough for storing scores in the database.
            return np.float32
        elif self.is_rankings_db:
            # Rankings databases store the zero-based rankings as optimally as possible in a:
            #   - 16-bit signed integer: max value = 2^15 - 1 = 32767 ==> can store 32768 rankings.
            #   - 32-bit signed integer: max value = 2^31 - 1 = 2147483647 ==> can store 2147483648
            if nbr_rows <= 2 ** 15:
                # Range int16: -2^15 (= -32768) to 2^15 - 1 (= 32767).
                return np.int16
            else:
                # Range int32: -2^31 (= -2147483648) to 2^31 - 1 (= 2147483647).
                return np.int32

