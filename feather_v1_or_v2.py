import numpy as np
import pyarrow.dataset as ds

from typing import List, Optional

# To get column names from a Feather v1 file without needing to load the whole Feather file first,
# get FlatBuffer Feather v1 schema and compile it, so we can read all column names without loading all data.
#
#   - Download FlatBuffer Feather v1 schema.
#         wget -O feather_v1.fbs https://github.com/apache/arrow/raw/master/cpp/src/arrow/ipc/feather.fbs
#
#   - Modify namespace.
#         sed -i -e 's/^namespace arrow.ipc.feather.fbs/namespace feather_v1_fbs/' feather_v1.fbs
#
#   - Create python code for FlatBuffer Feather v1 schema.
#         flatc --python feather_v1.fbs
#
# This will not be needed anymore if:
#   https://issues.apache.org/jira/projects/ARROW/issues/ARROW-10344
# gets resolved.
import feather_v1_fbs.CTable as feather_v1_fbs


def is_feather_v1_or_v2(feather_file: str) -> Optional[int]:
    """
    Check if the passed filename is a Feather v1 or v2 file.

    :param feather_file: Feather v1 or v2 filename.
    :return: 1 (for Feather version 1), 2 (for Feather version 2) or None
    """

    with open(feather_file, 'rb') as fh_feather:
        # Read first 6 and last 6 bytes to see if we have a Feather v2 file.
        fh_feather.seek(0, 0)
        feather_v2_magic_bytes_header = fh_feather.read(6)
        fh_feather.seek(-6, 2)
        feather_v2_magic_bytes_footer = fh_feather.read(6)

        if feather_v2_magic_bytes_header == feather_v2_magic_bytes_footer == b'ARROW1':
            return 2

        # Read first 4 and last 4 bytes to see if we have a Feather v1 file.
        feather_v1_magic_bytes_header = feather_v2_magic_bytes_header[0:4]
        feather_v1_magic_bytes_footer = feather_v2_magic_bytes_footer[2:]

        if feather_v1_magic_bytes_header == feather_v1_magic_bytes_footer == b'FEA1':
            return 1

    return None


def get_all_column_names_from_feather(feather_file: str) -> List:
    """
    Get all column names from a Feather v1 or v2 filename.

    :param feather_file: Feather v1 or v2 filename.
    :return: list of column names.
    """

    feather_v1_or_v2 = is_feather_v1_or_v2(feather_file)

    if feather_v1_or_v2 == 1:
        with open(feather_file, 'rb') as fh_feather:
            fh_feather.seek(-8, 2)

            # Get Feather v1 metadata length.
            metadata_length = np.frombuffer(fh_feather.read(4), dtype=np.int32)[0]

            # Read Feather v1 metadata.
            fh_feather.seek(- (metadata_length + 8), 2)
            feather_metadata = feather_v1_fbs.CTable.GetRootAs(bytearray(fh_feather.read(metadata_length)), 0)

            num_columns = feather_metadata.ColumnsLength()

            column_names = [
                feather_metadata.Columns(column_idx).Name().decode('utf-8')
                for column_idx in range(0, num_columns)
            ]
    elif feather_v1_or_v2 == 2:
        feather_v2_dataset = ds.dataset(feather_file, format="feather")
        column_names = feather_v2_dataset.schema.names
    else:
        raise ValueError(f'"{feather_file}" is not a Feather v1 or v2 file.')

    return column_names
