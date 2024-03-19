"""Base for the metadata module."""
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

SUPPORTED_COLUMN_TYPES = {"categorical", "numeric"}

SUPPORTED_ID_HEADERS = {
    "case_insensitive": {
        "id",
        "sampleid",
        "sample id",
        "sample-id",
        "featureid",
        "feature id",
        "feature-id",
    },
    # For backwards-compatibility with existing formats.
    "exact_match": {
        # QIIME 1 mapping files. "#Sample ID" was never supported, but
        # we're including it here for symmetry with the other supported
        # headers that allow a space between words.
        "#SampleID",
        "#Sample ID",
        # biom-format: observation metadata and "classic" (TSV) OTU tables.
        "#OTUID",
        "#OTU ID",
        # Qiita sample/prep information files.
        "sample_name",
    },
}

FORMATTED_ID_HEADERS = "Case-insensitive: %s\n\nCase-sensitive: %s" % (
    ", ".join(repr(e) for e in sorted(SUPPORTED_ID_HEADERS["case_insensitive"])),
    ", ".join(repr(e) for e in sorted(SUPPORTED_ID_HEADERS["exact_match"])),
)


def is_id_header(name):
    """Determine if a name is a valid ID column header.

    This function may be used to determine if a value in a metadata file is a
    valid ID column header, or if a pandas ``Index.name`` matches the ID header
    requirements. The "ID header" corresponds to the ``Metadata.id_header``
    and ``MetadataColumn.id_header`` properties.

    Parameters
    ----------
    name : string or None
        Name to check against ID header requirements.

    Returns
    -------
    bool
        ``True`` if `name` is a valid ID column header, ``False`` otherwise.

    """
    return name and (
        name in SUPPORTED_ID_HEADERS["exact_match"]
        or name.lower() in SUPPORTED_ID_HEADERS["case_insensitive"]
    )
