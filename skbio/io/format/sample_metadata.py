"""
Sample Metadata object ported over from qiime2
===============================================

.. currentmodule:: skbio.io.format.sample_metadata

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.metadata.SampleMetadata`                              |
+------+------+---------------------------------------------------------------+
"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import csv
import re

from skbio.io import create_format
from skbio.metadata._metadata import SampleMetadata
from skbio.metadata.io import MetadataReader, MetadataWriter


sample_metadata = create_format("sample_metadata")


@sample_metadata.sniffer()
def _sample_metadata_sniffer(fh):
    # Strategy:
    # Check if first word is in the file is in the list
    # of allowed metadata words
    try:
        tsv_reader = csv.reader(fh, dialect="excel-tab", strict=True)
        # sample id and feature id are not separated when reading the tsv
        # since they are not tab-separated.
        possible_ids = [
            "id",
            "sampleid",
            "sample-id",
            "featureid",
            "feature-id",
            "sample id",
            "feature id",
        ]

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

            if any(
                [x.casefold() == header[0].strip().casefold() for x in possible_ids]
            ):
                return True, {}

            # if the first non-empty non-comment row doesn't have a valid id as
            # first entry we conclude that this is not a metadata file.
            return False, {}

        # In case the file is empty and has no files that non-empty non-comment
        # we return a negative result.
        return False, {}

    # if we run into errors with the csv file we assume its not a metadata file
    except csv.Error:
        return False, {}


@sample_metadata.reader(SampleMetadata)
def _sample_metadata_read(fh):
    return MetadataReader(fh).read(SampleMetadata)


@sample_metadata.writer(SampleMetadata)
def _sample_metadata_write(obj, fh):
    MetadataWriter(obj).write(fh)
