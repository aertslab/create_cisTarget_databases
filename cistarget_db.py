import numpy as np
import pandas as pd
import pyarrow.feather as pf

from enum import Enum, unique
from typing import List, Optional, Set, Tuple, Type, Union


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

        self.ids = tuple(sorted(set(feature_ids)))
        self.type = features_type

    def __repr__(self) -> str:
        return f'FeatureIDs(\n  feature_ids={self.ids},\n  features_type={self.type}\n)'

    def __eq__(self, other: object):
        if not isinstance(other, FeatureIDs):
            return NotImplemented

        return self.type == other.type and self.ids == other.ids

    def __len__(self) -> int:
        return len(self.ids)


class MotifOrTrackIDs:
    """
    MotifOrTrackIDs class represents a unique sorted tuple of motif IDs or track IDs for constructing a Pandas
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

        self.ids = tuple(sorted(set(motif_or_track_ids)))
        self.type = motifs_or_tracks_type

    def __repr__(self) -> str:
        return f'MotifOrTrackIDs(\n  motif_or_track_ids={self.ids},\n  motifs_or_tracks_type={self.type}\n)'

    def __eq__(self, other: object):
        if not isinstance(other, MotifOrTrackIDs):
            return NotImplemented

        return self.type == other.type and self.ids == other.ids

    def __len__(self) -> int:
        return len(self.ids)


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

        assert False, f'"features_type" is not handled for {self}'

    @property
    def motifs_or_tracks_type(self) -> 'MotifsOrTracksType':
        """Return MotifsOrTracksType Enum member for DatabaseTypes member."""
        if self.is_motifs_db:
            return MotifsOrTracksType.MOTIFS
        elif self.is_tracks_db:
            return MotifsOrTracksType.TRACKS

        assert False, f'"motifs_or_tracks_type" is not handled for {self}'

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

        assert False, f'"get_dtype" is not handled for {self}'


class CisTargetDatabase:
    """
    Class to create/update/read/write cisTarget databases.
    """

    @staticmethod
    def create_db(db_type: Union['DatabaseTypes', str],
                  feature_ids: FeatureIDs,
                  motif_or_track_ids: MotifOrTrackIDs) -> 'CisTargetDatabase':
        """
        Create zeroed cisTarget scores or rankings database of one of the DatabaseTypes types
        for FeatureIDs vs MotifOrTrackIDs or vice versa.

        :param db_type: Type of database.
        :param feature_ids: FeatureIDs object or list, set or tuple of feature IDs.
        :param motif_or_track_ids: MotifOrTrackIDs object or list, set or tuple of motif or track IDs.
        :return: CisTargetDatabase object.
        """

        if not isinstance(db_type, DatabaseTypes):
            if isinstance(db_type, str):
                # If the database type was given as a string, try to convert it to a member of DatabaseTypes Enum.
                try:
                    db_type = DatabaseTypes.from_str(database_type=db_type)
                except ValueError as e:
                    raise e
            else:
                raise ValueError('db_type must be of "DatabaseTypes" type.')

        if not isinstance(feature_ids, FeatureIDs):
            if isinstance(feature_ids, List) or isinstance(feature_ids, Set) or isinstance(feature_ids, Tuple):
                # If feature IDs were given as a list, set or tuple, convert it to a FeatureIDs object.
                feature_ids = FeatureIDs(feature_ids=feature_ids, features_type=db_type.features_type)
            else:
                raise ValueError('feature_ids must be of "FeatureIDs" type.')
        else:
            if feature_ids.type != db_type.features_type:
                raise ValueError(f'"feature_ids" type ({feature_ids.type}) is not of the same as the one defined for "db_type" ({db_type.features_type}).')

        if not isinstance(motif_or_track_ids, MotifOrTrackIDs):
            if isinstance(motif_or_track_ids, List) or isinstance(motif_or_track_ids, Set) or isinstance(motif_or_track_ids, Tuple):
                # If motif or track IDs were given as a list, set or tuple, convert it to a MotifOrTrackIDs object.
                motif_or_track_ids = MotifOrTrackIDs(
                    motif_or_track_ids=motif_or_track_ids,
                    motifs_or_tracks_type=db_type.motifs_or_tracks_type
                )
            else:
                raise ValueError('motif_or_track_ids must be of "MotifOrTrackIDs" type.')
        else:
            if motif_or_track_ids.type != db_type.motifs_or_tracks_type:
                raise ValueError(f'"motif_or_track_ids" type ({motif_or_track_ids.type}) is not of the same as the one defined for "db_type" ({db_type.motifs_or_tracks_type}).')

        # Create feature IDs index and motif or track IDs index.
        feature_ids_idx = pd.Index(
            data=feature_ids.ids,
            name=feature_ids.type.value
        )
        motif_or_track_ids_idx = pd.Index(
            data=motif_or_track_ids.ids,
            name=motif_or_track_ids.type.value
        )

        if db_type.column_kind == 'regions' or db_type.column_kind == 'genes':
            # Create zeroed dataframe for all region IDs or gene IDs vs all motif or track IDs.
            df = pd.DataFrame(
                data=np.zeros((len(motif_or_track_ids), len(feature_ids)),
                              dtype=db_type.get_dtype(nbr_rows=len(feature_ids))),
                index=motif_or_track_ids_idx,
                columns=feature_ids_idx
            )
        elif db_type.column_kind == 'motifs' or db_type.column_kind == 'tracks':
            # Create zeroed dataframe for all motif or track IDs vs all region IDs or gene IDs.
            df = pd.DataFrame(
                data=np.zeros((len(feature_ids), len(motif_or_track_ids)),
                              dtype=db_type.get_dtype(nbr_rows=len(motif_or_track_ids))),
                index=feature_ids_idx,
                columns=motif_or_track_ids_idx
            )

        return CisTargetDatabase(db_type, df)

    @staticmethod
    def read_db(db_filename: str, db_type: Optional[Union['DatabaseTypes', str]] = None) -> 'CisTargetDatabase':
        """
        Read cisTarget database from Feather file to CisTargetDatabase object.

        :param db_filename: Feather database filename.
        :param db_type: Type of database (can be automatically determined from the filename if written with write_db).
        :return: CisTargetDatabase object.
        """

        assert db_filename is not None

        # Try to extract the database type from database filename if database type was not specified.
        if not db_type:
            try:
                db_type, db_prefix, extension = \
                    DatabaseTypes.create_database_type_and_db_prefix_and_extension_from_db_filename(
                        db_filename=db_filename
                    )
            except ValueError:
                raise ValueError(
                    'cisTarget database type could not be automatically determined from the db_filename. '
                    'Specify db_type instead.'
                )
        else:
            if not isinstance(db_type, DatabaseTypes):
                if isinstance(db_type, str):
                    # If the database type was given as a string, try to convert it to a member of DatabaseTypes Enum.
                    try:
                        db_type = DatabaseTypes.from_str(database_type=db_type)
                    except ValueError as e:
                        raise e
                else:
                    raise ValueError('db_type must be of "DatabaseTypes" type.')

        # Read Feather file in dataframe.
        df = pf.read_feather(source=db_filename, columns=None)

        # Set indexes and index names for rows.
        if db_type.row_kind in df:
            df.set_index(db_type.row_kind, inplace=True)
        elif 'index' in df:
            df.set_index('index', inplace=True)
            df.index.set_names([db_type.row_kind], inplace=True)
        else:
            raise ValueError(
                f'cisTarget database file "{db_filename}" does not contain a "{db_type.row_kind}" or "index" column.'
            )

        # Set column name.
        df.columns.set_names([db_type.column_kind], inplace=True)

        return CisTargetDatabase(db_type, df)

    def __init__(self, db_type: DatabaseTypes, df: pd.DataFrame):
        self.db_type: DatabaseTypes = db_type
        self.df: pd.DataFrame = df

    @property
    def feature_ids(self) -> FeatureIDs:
        """
        Get Feature IDs present in the cisTarget database.

        :return: feature_ids
        """

        if self.db_type.column_kind == 'regions' or self.db_type.column_kind == 'genes':
            feature_ids = FeatureIDs(
                feature_ids=self.df.columns.to_list(),
                features_type=FeaturesType.from_str(feature_type=self.db_type.column_kind)
            )
        elif self.db_type.row_kind == 'regions' or self.db_type.row_kind == 'genes':
            feature_ids = FeatureIDs(
                feature_ids=self.df.index.to_list(),
                features_type=FeaturesType.from_str(feature_type=self.db_type.row_kind)
            )

        return feature_ids

    @property
    def motif_or_track_ids(self) -> MotifOrTrackIDs:
        """
        Get MotifOrTrackIDs present in the cisTarget database.

        :return: motif_or_track_ids
        """

        if self.db_type.column_kind == 'motifs' or self.db_type.column_kind == 'tracks':
            motif_or_track_ids = MotifOrTrackIDs(
                motif_or_track_ids=self.df.columns.to_list(),
                motifs_or_tracks_type=MotifsOrTracksType.from_str(motifs_or_tracks_type=self.db_type.column_kind)
            )
        elif self.db_type.row_kind == 'motifs' or self.db_type.row_kind == 'tracks':
            motif_or_track_ids = MotifOrTrackIDs(
                motif_or_track_ids=self.df.index.to_list(),
                motifs_or_tracks_type=MotifsOrTracksType.from_str(motifs_or_tracks_type=self.db_type.row_kind)
            )

        return motif_or_track_ids

    @property
    def dtype(self) -> Union[np.float32, np.int16, np.int32]:
        """
        Get dtype of scores or rankings stored in the cisTarget database.

        :return: numpy dtype
        """

        return self.df.iloc[0, 0].dtype.type

    @property
    def shape(self) -> Tuple[int, int]:
        """
        Get shape of cisTarget database Dataframe.

        :return: shape
        """

        return self.df.shape

    @property
    def nbr_rows(self) -> int:
        """
        Get number of rows in cisTarget database.

        :return: nbr_rows
        """

        return self.df.shape[0]

    @property
    def nbr_columns(self) -> int:
        """
        Get number of columns in cisTarget database.

        :return: nbr_columns
        """

        return self.df.shape[1]

    def write_db(self, db_prefix: Optional[str] = None, db_filename: Optional[str] = None) -> str:
        """
        Write cisTarget database to Feather file.
        If db_prefix is used, the database type will be encoded in the Feather database filename.
        If db_filename is used instead, that exact filename will be used.

        :param db_prefix: Database prefix path.
        :param db_filename: Full database filename.
        :return: db_filename: cisTarget database filename (constructed from db_prefix).
        """

        assert isinstance(db_prefix, str) and db_filename is None or db_prefix is None and isinstance(db_filename, str)

        if db_prefix:
            db_filename = self.create_db_filename_from_db_prefix(db_prefix=db_prefix, extension='feather')

        assert isinstance(db_filename, str)

        # Temporarily add the index column with the name of the row kind to the dataframe,
        # so row names of the dataframe are written to the Feather file.
        self.df[self.db_type.row_kind] = self.df.index.to_series()

        # Write cisTarget database in Feather v2 format with lz4 compression.
        pf.write_feather(df=self.df, dest=db_filename, compression='lz4', version=2)

        # Delete index column with row names from the dataframe.
        del self.df[self.db_type.row_kind]

        # Add names of row and column indexes again.
        self.df.index.set_names([self.db_type.row_kind], inplace=True)
        self.df.columns.set_names([self.db_type.column_kind], inplace=True)

        return db_filename

    def create_db_filename_from_db_prefix(self, db_prefix: str, extension: str = 'feather') -> str:
        """
        Create cisTarget database filename based on a database prefix string and extension.
        Between the database prefix string and extension the database type will be encoded.

        :param db_prefix: Database prefix.
        :param extension: Database extension.
        :return: Database filename.
        """

        return self.db_type.create_db_filename(
            db_prefix=db_prefix,
            extension=extension
        )

    def transpose(self, copy: bool = False):
        """
        Transpose cisTarget database (switch rows and columns).

        For example transpose a cisTarget database of DatabaseTypes.SCORES_DB_MOTIFS_VS_REGIONS
        to DatabaseTypes.SCORES_DB_REGIONS_VS_MOTIFS.

        :param copy: Copy data or not.
        :return: Transposed CisTargetDatabase object.
        """

        # Transpose dataframe (switch rows and columns) and change database type.
        return CisTargetDatabase(
            db_type=DatabaseTypes.from_strings(
                scores_or_rankings=self.db_type.scores_or_rankings,
                column_kind=self.db_type.row_kind,
                row_kind=self.db_type.column_kind
            ),
            df=self.df.transpose(copy=copy)
        )

    def update_scores_for_motif_or_track(self, motif_or_track_id: str, df_scores_for_motif_or_track: pd.DataFrame):
        """
        Update CRM or track scores for 1 motif or track ID for specific region or gene IDs in cisTarget scores database.

        :param motif_or_track_id:
            motif or track ID for which scores are provided.
        :param df_scores_for_motif_or_track:
            Dataframe with region IDs or gene IDs as index
            and motif or track scores in the first column (and only column) or in 'crm_score' or 'track_score' column.
        :return:
        """

        assert self.db_type.is_scores_db, 'cisTarget database must be a scores database.'

        # Write CRM or track scores for motif ID or track ID to cisTarget scores dataframe.
        if self.db_type.column_kind == 'regions' or self.db_type.column_kind == 'genes':
            assert motif_or_track_id in self.df.index, \
                f'"{motif_or_track_id}" not found in row of CisTargetDatabase dataframe.'

            # Column names contain regions or genes.
            if df_scores_for_motif_or_track.shape[1] == 1:
                self.df.loc[
                    motif_or_track_id,
                    df_scores_for_motif_or_track.index,
                ] = df_scores_for_motif_or_track.iloc[:, 0]
            else:
                score_name = 'crm_score' if self.db_type.is_motifs_db else 'track_score'
                self.df.loc[
                    motif_or_track_id,
                    df_scores_for_motif_or_track.index,
                ] = df_scores_for_motif_or_track[score_name]
        elif self.db_type.column_kind == 'motifs' or self.db_type.column_kind == 'tracks':
            assert motif_or_track_id in self.df.columns, \
                f'"{motif_or_track_id}" not found in column of CisTargetDatabase dataframe.'

            # Row names contain regions or genes.
            if df_scores_for_motif_or_track.shape[1] == 1:
                self.df.loc[
                    df_scores_for_motif_or_track.index,
                    motif_or_track_id
                ] = df_scores_for_motif_or_track.iloc[:, 0]
            else:
                score_name = 'crm_score' if self.db_type.is_motifs_db else 'track_score'
                self.df.loc[
                    df_scores_for_motif_or_track.index,
                    motif_or_track_id
                ] = df_scores_for_motif_or_track[score_name]

    def convert_scores_db_to_rankings_db(self, rand_seed: Optional[int] = None) -> 'CisTargetDatabase':
        """
        Convert scores cisTarget database to rankings database.

        :param rand_seed: Set seed for random number generator if generating the same rankings database
                          from the same scores database is required.
                          If set to None, unpredictable entropy will be pulled from the OS.
        :return: rankings cisTarget database
        """

        assert self.db_type.is_scores_db, 'cisTarget database must be a scores database.'

        # Initialize random number generator, so same rankings database can be generated if the same seed is set.
        rng = np.random.default_rng(seed=rand_seed)

        def rank_scores_and_assign_random_ranking_in_range_for_ties(scores_with_ties_for_motif_or_track_numpy):
            # Create random permutation so tied scores will have a different ranking each time.
            random_permutations_to_break_ties_numpy = rng.permutation(
                scores_with_ties_for_motif_or_track_numpy.shape[0]
            )

            # Rank scores for each region/gene for the current motif/track and break ties:
            #   - Get scores for each region/gene for a certain motif/track and multiply by -1 (scores >= 0) so sorting
            #     later will result in ranking the highest score first (descending):
            #
            #       (-scores_with_ties_for_motif_or_track_numpy)
            #
            #   - Access the negated scores in a random order:
            #
            #       (-scores_with_ties_for_motif_or_track_numpy)[random_permutations_to_break_ties_numpy]
            #
            #     so when sorting it, scores for regions/genes at the start of the array do not get artificially better
            #     rankings than regions/genes more at the bottom of the array (as argsort works on a first come, first
            #     served basis).
            #
            #   - Sort the negated scores (accessed in a random order) and return an array with indices that would sort
            #     those scores from high to low (first position in the returned array contains the index to the value in
            #     scores_with_ties_for_motif_or_track_numpy with the highest score):
            #
            #       (-scores_with_ties_for_motif_or_track_numpy)[random_permutations_to_break_ties_numpy].argsort()
            #
            #   - Undo the random order access of the array created in the previous step, so the indices that would sort
            #     scores_with_ties_for_motif_or_track_numpy from high scores to low scores correspond again with the
            #     input array:
            #
            #       random_permutations_to_break_ties_numpy[
            #           (-scores_with_ties_for_motif_or_track_numpy)[random_permutations_to_break_ties_numpy].argsort()
            #       ]
            #
            #   - Finally convert the array (previous step) which contains indices which would sort
            #     scores_with_ties_for_motif_or_track_numpy from high scores to low scores and which would break tied
            #     scores in a fair (random) way to a ranking (and store the result in a int16 or int32 numpy array):
            #
            #       ... .argsort().astype(rankings_db_dtype)
            #
            ranking_with_broken_ties_for_motif_or_track_numpy = random_permutations_to_break_ties_numpy[
                (-scores_with_ties_for_motif_or_track_numpy)[random_permutations_to_break_ties_numpy].argsort()
            ].argsort().astype(rankings_db_dtype)

            return ranking_with_broken_ties_for_motif_or_track_numpy

        # Create zeroed rankings database.
        rankings_db = CisTargetDatabase.create_db(
            db_type=DatabaseTypes.from_strings(
                scores_or_rankings='rankings',
                column_kind=self.db_type.column_kind,
                row_kind=self.db_type.row_kind
            ),
            feature_ids=self.feature_ids,
            motif_or_track_ids=self.motif_or_track_ids
        )

        # Get dtype of rankings cisTarget database as this is used to return a numpy array of the correct type in
        # rank_scores_and_assign_random_ranking_in_range_for_ties function.
        rankings_db_dtype = rankings_db.dtype

        # Rank all scores per motif/track and assign a random ranking in range for regions/genes with the same score.
        if self.db_type.column_kind == 'regions' or self.db_type.column_kind == 'genes':
            for row_idx in range(self.nbr_rows):
                rankings_db.df.iloc[row_idx, :] = rank_scores_and_assign_random_ranking_in_range_for_ties(
                    self.df.iloc[row_idx, :].to_numpy()
                )
        elif self.db_type.column_kind == 'motifs' or self.db_type.column_kind == 'tracks':
            for column_idx in range(self.nbr_columns):
                rankings_db.df.iloc[:, column_idx] = rank_scores_and_assign_random_ranking_in_range_for_ties(
                    self.df.iloc[:, column_idx].to_numpy()
                )

        return rankings_db
