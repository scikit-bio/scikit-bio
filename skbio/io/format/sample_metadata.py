"""Sample Metadata object ported over from qiime2.

===============================================

.. currentmodule:: skbio.io.format.sample_metadata

This implements the Sample_Metadata format which is identical to the
Metadata format implemented in qiime2.
(see: https://docs.qiime2.org/2024.2/tutorials/metadata/)

An example sample_metadata file:

.. code-block:: none

    id	col1	col2	col3
    #q2:types	categorical	categorical	categorical
    id1	1	a	foo
    id2	2	b	bar
    id3	3	c	42


Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.metadata.SampleMetadata`                           |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------

Metadata Formatting Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

QIIME 2 metadata is most commonly stored in a TSV (i.e. tab-separated values)
file. These files typically have a .tsv or .txt file extension, though it
doesn't matter to QIIME 2 what file extension is used. TSV files are simple
text files used to store tabular data, and the format is supported by many
types of software, such as editing, importing, and exporting from spreadsheet
programs and databases. Thus, it's usually straightforward to manipulate
QIIME 2 metadata using the software of your choosing. If in doubt, we recommend
using a spreadsheet program such as Microsoft Excel or Google Sheets to edit
and export your metadata files.

The following sections describe formatting requirements for QIIME 2 metadata
files, and how to validate your metadata files. Since there is no universal
standard for TSV files, it is important to adhere to these requirements and
understand how QIIME 2 will interpret the file's contents to get the most out
of your (meta)data!

Metadata Validation
^^^^^^^^^^^^^^^^^^^

Sample and feature metadata files stored in Google Sheets can be validated
using Keemei. Select Add-ons > Keemei > Validate QIIME 2 metadata file to
validate metadata stored in Google Sheets.

QIIME 2 will also automatically validate a metadata file anytime it is used by
the software. However, using Keemei to validate your metadata is recommended
because a report of all validation errors and warnings will be presented each
time Keemei is run. Loading your metadata in QIIME 2 will typically present
only a single error at a time, which can make identifying and resolving
validation issues cumbersome, especially if there are many issues with the
metadata.


Leading and trailing whitespace characters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If any cell in the metadata contains leading or trailing whitespace characters
(e.g. spaces, tabs), those characters will be ignored when the file is loaded.
Thus, leading and trailing whitespace characters are not significant, so cells
containing the values 'gut' and '  gut  ' are equivalent. This rule is applied
before any other rules described below.


Comments and Empty Rows
^^^^^^^^^^^^^^^^^^^^^^^

Rows whose first cell begins with the pound sign (#) are interpreted as
comments and may appear anywhere in the file. Comment rows are ignored by
QIIME 2 and are for informational purposes only. Inline comments are not
supported.

Empty rows (e.g. blank lines or rows consisting solely of empty cells) may
appear anywhere in the file and are ignored.

Identifier Column
^^^^^^^^^^^^^^^^^

The first column in the metadata file is the identifier (ID) column. This
column defines the sample or feature IDs associated with your study. It is not
recommended to mix sample and feature IDs in a single metadata file; keep
sample and feature metadata stored in separate files.

The ID column name (i.e. ID header) must be one of the following values. The
values listed below may not be used to name other IDs or columns in the file.

Case-insensitive:

- id

- sampleid

- sample id

- sample-id

- fetureid

- feature id

- feature-id

Case-sensitive (these are mostly for backwards-compatibility with QIIME 1,
biom-format, and Qiita files):

- #SampleID

- #Sample ID

- #OTUID

- #OTU ID

- sample_name

The following rules apply to IDs:

- IDs may consist of any Unicode characters, with the exception that IDs must
  notstart with the pound sign (#), as those rows would be interpreted as comments
  and ignored. See the section Recommendations for Identifiers for
  recommendations on choosing identifiers in your study.

- IDs cannot be empty (i.e. they must consist of at least one character).

- IDs must be unique (exact string matching is performed to detect duplicates).

- At least one ID must be present in the file.

- IDs cannot use any of the reserved ID column names listed above.

Recommendations for Identifiers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Our goal with QIIME 2 is to support arbitrary Unicode characters in all cells
of metadata files. However, given that QIIME 2 plugins and interfaces can be
developed by anyone, we can’t make a guarantee that arbitrary Unicode
characters will work with all plugins and interfaces. We can therefore make
recommendations to users about characters that should be safe to use in
identifiers, and we are preparing resources for plugin and interface developers
to help them make their software as robust as possible. As developer resources
become available, we will announce them in the Developer Discussion category
on the QIIME 2 Forum.

Sample and feature identifiers with problematic characters tend to cause the
most issues for our users. Based on our experiences with QIIME 1, QIIME 2, and
other bioinformatics and command line tools, we can recommend the following
attributes for identifiers:

- Identifiers should be 36 characters long or less.

- Identifiers should contain only ASCII alphanumeric characters
  (i.e. in the range of [a-z], [A-Z], or [0-9]), the period (.) character, or
  the dash (-) character.

An important point to remember is that sometimes values in your sample metadata
can become identifiers. For example, taxonomy annotations can become feature
identifiers following qiime taxa collapse, and sample or feature metadata
values can become identifiers after applying qiime feature-table group.
If you plan to apply these or similar methods where metadata values can become
identifiers, you will be less likely to encounter problems if the values adhere
to these identifier recommendations as well.

To help users become aware of these recommendations, the Keemei metadata
validator will warn users about identifiers that don’t meet the above
recommendations.

Users may be interested in the cual-id software for assistance with creating
sample identifiers. The cual-id paper also provides some discussion on how to
design identifiers.

Metadata Columns
^^^^^^^^^^^^^^^^

The ID column is the first column in the metadata file, and can optionally be
followed by additional columns defining metadata associated with each sample or
feature ID. Metadata files are not required to have additional metadata
columns, so a file containing only an ID column is a valid QIIME 2 metadata
file.

The following rules apply to column names:

- May consist of any Unicode characters.

- Cannot be empty (i.e. column names must consist of at least one character).

- Must be unique (exact string matching is performed to detect duplicates).

- Column names cannot use any of the reserved ID column names described in the
  section Identifier Column.

The following rules apply to column values:

- May consist of any Unicode characters.

- Empty cells represent missing data. Other values such as NA are not
  interpreted as missing data; only the empty cell is recognized as “missing”.
  Note that cells consisting solely of whitespace characters are also
  interpreted as missing data because leading and trailing whitespace
  characters are always ignored, effectively making the cell empty.

Column Types
^^^^^^^^^^^^

QIIME 2 currently supports categorical and numeric metadata columns. By
default, QIIME 2 will attempt to infer the type of each metadata column: if the
column consists only of numbers or missing data, the column is inferred to be
numeric. Otherwise, if the column contains any non-numeric values, the column
is inferred to be categorical. Missing data (i.e. empty cells) are supported in
categorical columns as well as numeric columns.

QIIME 2 supports an optional comment directive to allow users to explicitly
state a column's type, avoiding the column type inference described above.
This can be useful if there is a column that appears to be numeric, but should
actually be treated as categorical metadata (e.g. a Subject column where
subjects are labeled 1, 2, 3, etc). Explicitly declaring a column's type also
makes your metadata file more descriptive because the intended column type is
included with the metadata, instead of relying on software to infer the type
(which isn't always transparent).

You can use an optional comment directive to declare column types in your
metadata file, either manually or through the q2cli developer tools.

For manual specifications within your metadata file(s), the comment directive
must appear directly below the header. The row's first cell must be #q2:types
or #sk:types to indicate the row is a comment directive. Subsequent cells may
contain the values categorical or numeric (both case-insensitive).
The empty cell is also supported if you do not wish to assign a type to a
column (the type will be inferred in that case). Thus, it is easy to include
this comment directive without having to declare types for every column in
your metadata.


Number Formatting
^^^^^^^^^^^^^^^^^

If a column is to be interpreted as a numeric metadata column (either through
column type inference or by using the #q2:types comment directive), numbers in
the column must be formatted following these rules:

- Use the decimal number system: ASCII characters [0-9], . for an optional
  decimal point, and + and - for positive and negative signs, respectively.

    - Examples: 123, 123.45, 0123.40, -0.000123, +1.23

- Scientific notation may be used with E-notation; both e and E are supported.

    - Examples: 1e9, 1.23E-4, -1.2e-08, +4.5E+6

- Only up to 15 digits total (including before and after the decimal point) are
  supported to stay within the 64-bit floating point specification. Numbers
  exceeding 15 total digits are unsupported and will result in undefined
  behavior.

- Common representations of not a number (e.g. NaN, nan) or infinity
  (e.g. inf, -Infinity) are not supported. Use an empty cell for missing data
  (e.g. instead of NaN). Infinity is not supported at this time in QIIME 2
  metadata files.

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


# newline="" is because otherwise csv.writer will write blank lines between rows
# in Windows. See: https://stackoverflow.com/questions/3348460/
sample_metadata = create_format("sample_metadata", newline="")


@sample_metadata.sniffer()
def _sample_metadata_sniffer(fh):
    # Strategy:
    # Check if first word in the file is in the list
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
            "sample_name",  # Technically this should be case-sensitive
        ]
        possible_ids_w_leading_comment_char = [
            "#SampleID",
            "#Sample ID",
            "#OTUID",
            "#OTU ID",
        ]

        # We need to find the actual header row
        # so we loop until we find the first row that isn't empty or a comment
        for header in tsv_reader:
            # Skip empty rows
            if len(header) == 0:
                continue

            match = re.search(r"\S+", header[0])

            # Check if first word is a columnID that starts with #
            if match and match.group() in possible_ids_w_leading_comment_char:
                return True, {}

            # Skip rows whose first non-whitespace character is a #
            # since they are comments. skips empty rows too.
            if not match or match.group()[0] == "#":
                continue

            if any(
                [x.casefold() == header[0].strip().casefold() for x in possible_ids]
            ):
                return True, {}

            # if the first non-empty non-comment row doesn't have a valid id as
            # first entry we conclude that this is not a metadata file.
            return False, {}

        # In case the file is empty and has no rows that are non-empty non-comment
        # we return a negative result.
        return False, {}

    # if we run into errors with the csv file we assume its not a metadata file
    except csv.Error:
        return False, {}


@sample_metadata.reader(SampleMetadata)
def _sample_metadata_read(fh, **kwargs):
    return MetadataReader(fh).read(SampleMetadata, **kwargs)


@sample_metadata.writer(SampleMetadata)
def _sample_metadata_write(obj, fh):
    MetadataWriter(obj).write(fh)
