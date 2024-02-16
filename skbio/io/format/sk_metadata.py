from skbio.io import create_format
from skbio.metadata._metadata import Metadata, MetadataColumn
from skbio.metadata.io import MetadataFileError, MetadataReader
from skbio.metadata.base import (
    SUPPORTED_COLUMN_TYPES,
    FORMATTED_ID_HEADERS,
    is_id_header,
)
from skbio.metadata._util import find_duplicates
from skbio.metadata.missing import (
    DEFAULT_MISSING,
    BUILTIN_MISSING,
    series_encode_missing,
)

import csv
import itertools
import pandas as pd


sk_metadata = create_format("sk_metadata")


@sk_metadata.reader(Metadata)
def _sk_metadata_read(fh):
    return MetadataReader(fh).read(Metadata)


@sk_metadata.writer(Metadata, monkey_patch=False)
def _sk_metadata_write(obj, fh):
    MetadataWriter(obj).write(fh)
