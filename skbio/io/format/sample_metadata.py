from skbio.io import create_format
from skbio.metadata._metadata import SampleMetadata, MetadataColumn
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


sample_metadata = create_format("sample_metadata")


@sample_metadata.sniffer()
def _sample_metadata_sniffer(fh):
    # Strategy:
    # Check if first word is in the file is in the list
    # of allowed metadata words
    try:
        tsv_reader = csv.reader(fh, dialect="excel-tab", strict=True)
        for line in tsv_reader:
            if line[0] in [
                "id",
                "sampleid",
                "sample id",
                "sample-id",
                "featureid",
                "feature id",
                "feature-id",
            ]:
                return True, {}
            else:
                return False, {}
    except Exception:
        return False, {}


@sample_metadata.reader(SampleMetadata)
def _sk_metadata_read(fh):
    return MetadataReader(fh).read(SampleMetadata)


@sample_metadata.writer(SampleMetadata, monkey_patch=False)
def _sample_metadata_write(obj, fh):
    MetadataWriter(obj).write(fh)
