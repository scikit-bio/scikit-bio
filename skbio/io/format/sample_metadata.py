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
import re


sample_metadata = create_format("sample_metadata")


@sample_metadata.sniffer()
def _sample_metadata_sniffer(fh):
    # Strategy:
    # Check if first word is in the file is in the list
    # of allowed metadata words
    try:
        tsv_reader = csv.reader(fh, dialect="excel-tab", strict=True)
        possible_ids = ["id", "sampleid", "sample-id", "featureid", "feature-id"]
        possible_spaced_ids = ["sample id", "feature id"]

        # We need to find the actual header row
        # so we loop until we find the first row that isn't empty or a comment
        for header in tsv_reader:
            # Skip empty rows
            if len(header) == 0:
                continue

            # Skip rows whose first non-whitespace character is a #
            # since they are comments.
            match = re.search(r"\S", header[0])
            if not match or match.group() == "#":
                continue

            if any([x == header[0].strip() for x in possible_ids]):
                return True, {}

            if len(header) > 1 and any(
                [
                    x == " ".join(h.strip() for h in header[:2])
                    for x in possible_spaced_ids
                ]
            ):
                return True, {}

            # if the first non-empty non-comment row doesn't have a valid id as
            # first entry we conclude that this is not a metadata file.
            return False, {}

    # if we run into errors with the csv file we assume its not a metadata file
    except csv.Error:
        return False, {}


@sample_metadata.reader(SampleMetadata)
def _sk_metadata_read(fh):
    return MetadataReader(fh).read(SampleMetadata)


@sample_metadata.writer(SampleMetadata, monkey_patch=False)
def _sample_metadata_write(obj, fh):
    MetadataWriter(obj).write(fh)
