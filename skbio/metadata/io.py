"""Contains io functionality for the Metadata module."""
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import csv
import itertools
import os.path
import re

import numpy as np
import pandas as pd

from skbio.io._fileobject import SaneTextIOWrapper
from skbio.util import find_duplicates
from .missing import DEFAULT_MISSING, BUILTIN_MISSING, series_encode_missing
from .base import SUPPORTED_COLUMN_TYPES, FORMATTED_ID_HEADERS, is_id_header
from ..metadata._metadata import SampleMetadata, MetadataColumn


class MetadataFileError(Exception):
    """Exception for errors with Metadata files."""

    _suffix = (
        "There may be more errors present in the metadata file. To get a full "
        "report, sample/feature metadata files can be validated with Keemei: "
        "https://keemei.qiime2.org\n\nFind details on QIIME 2 metadata "
        "requirements here: https://docs.qiime2.org/"
    )

    def __init__(self, message, include_suffix=True):
        """Initialize the MetadataFileError."""
        # LH NOTE/TODO: in Qiime2 this linked to the specific Qiime2 release.
        # However since this is not Qiime2 It did break and I removed this

        if include_suffix:
            message = message + "\n\n" + self._suffix
        super().__init__(message)


class MetadataReader:
    """Reader for Metadata files."""

    def __init__(self, filepath_or_filehandle):
        """Initialize the Reader for Metadata files."""
        # check if the filepath_filehandle is a path... if it is check if it
        # points to a file
        # TODO: Refine this check to be more specific
        if isinstance(filepath_or_filehandle, str):
            self._file_is_filehandle = False
            if not os.path.isfile(filepath_or_filehandle):
                raise MetadataFileError(
                    "Metadata file path doesn't exist, or the path points to "
                    "something other than a file. Please check that the path "
                    "exists, has read permissions, and points to a regular file "
                    "(not a directory): %s" % filepath_or_filehandle
                )
        else:
            self._file_is_filehandle = True

        self._filepath = filepath_or_filehandle

        # Used by `read()` to store an iterator yielding rows with
        # leading/trailing whitespace stripped from their cells (this is a
        # preprocessing step that should happen with *every* row). The iterator
        # protocol is the only guaranteed API on this object.
        self._reader = None

    def read(
        self,
        into,
        column_types=None,
        column_missing_schemes=None,
        default_missing_scheme=DEFAULT_MISSING,
    ):
        """Return a Metadata object read from the given file."""
        if column_types is None:
            column_types = {}

        try:
            # choose the appropriate context manager depending
            # on if a filehandle has been passed.
            if self._file_is_filehandle:
                cm = self._filepath
            else:
                # Newline settings based on recommendation from csv docs:
                #     https://docs.python.org/3/library/csv.html#id3

                # Ignore BOM on read (but do not write BOM)
                cm = open(self._filepath, "r", newline="", encoding="utf-8-sig")

            with cm as fh:
                tsv_reader = csv.reader(fh, dialect="excel-tab", strict=True)
                self._reader = (self._strip_cell_whitespace(row) for row in tsv_reader)
                header = self._read_header()
                directives = self._read_directives(header)
                ids, data = self._read_data(header)
        except UnicodeDecodeError as e:
            if "0xff in position 0" in str(e) or "0xfe in position 0" in str(e):
                raise MetadataFileError(
                    "Metadata file must be encoded as UTF-8 or ASCII, found "
                    "UTF-16. If this file is from Microsoft Excel, save "
                    "as a plain text file, not 'UTF-16 Unicode'"
                )

            raise MetadataFileError(
                "Metadata file must be encoded as UTF-8 or ASCII. The "
                "following error occurred when decoding the file:\n\n%s" % e
            )
        finally:
            self._reader = None

        index = pd.Index(ids, name=header[0], dtype=object)
        df = pd.DataFrame(data, columns=header[1:], index=index, dtype=object)

        # TODO: move these checks over to Metadata.__init__() so that you can
        # pass column_types with an untyped dataframe. This would require a bit
        # of a refactor and doesn't buy a whole lot at the moment, hence the
        # TODO.
        for name, type in column_types.items():
            if name not in df.columns:
                raise MetadataFileError(
                    "Column name %r specified in `column_types` is not a "
                    "column in the metadata file." % name
                )
            if type not in SUPPORTED_COLUMN_TYPES:
                fmt_column_types = ", ".join(
                    repr(e) for e in sorted(SUPPORTED_COLUMN_TYPES)
                )
                raise MetadataFileError(
                    "Column name %r specified in `column_types` has an "
                    "unrecognized column type %r. Supported column types: %s"
                    % (name, type, fmt_column_types)
                )

        resolved_column_types = directives.get("types", {})
        resolved_column_types.update(column_types)

        if column_missing_schemes is None:
            column_missing_schemes = {}

        resolved_missing = {c: default_missing_scheme for c in df.columns}
        resolved_missing.update(directives.get("missing", {}))
        resolved_missing.update(column_missing_schemes)

        try:
            # Cast each column to the appropriate dtype based on column type.
            df = df.apply(
                self._cast_column,
                axis="index",
                column_types=resolved_column_types,
                missing_schemes=resolved_missing,
            )
        except MetadataFileError as e:
            # HACK: If an exception is raised within `DataFrame.apply`, pandas
            # adds an extra tuple element to `e.args`, making the original
            # error message difficult to read because a tuple is repr'd instead
            # of a string. To work around this, we catch and reraise a
            # MetadataFileError with the original error message. We use
            # `include_suffix=False` to avoid adding another suffix to the
            # error message we're reraising.
            msg = e.args[0]
            raise MetadataFileError(msg, include_suffix=False)

        try:
            return into(
                df,
                column_missing_schemes=resolved_missing,
                default_missing_scheme=default_missing_scheme,
            )
        except Exception as e:
            raise MetadataFileError(
                "There was an issue with loading the metadata file:\n\n%s" % e
            )

    def _read_header(self):
        header = None
        for row in self._reader:
            if self._is_header(row):
                header = row
                break
            elif self._is_comment(row):
                continue
            elif self._is_empty(row):
                continue
            elif self._is_directive(row):
                raise MetadataFileError(
                    "Found directive %r while searching for header. "
                    "Directives may only appear immediately after the header." % row[0]
                )
            else:
                raise MetadataFileError(
                    "Found unrecognized ID column name %r while searching for "
                    "header. The first column name in the header defines the "
                    "ID column, and must be one of these values:\n\n%s\n\n"
                    "NOTE: Metadata files must contain tab-separated values."
                    % (row[0], FORMATTED_ID_HEADERS)
                )

        if header is None:
            raise MetadataFileError(
                "Failed to locate header. The metadata file may be empty, or "
                "consists only of comments or empty rows."
            )

        # Trim trailing empty cells from header.
        data_extent = None
        for idx, cell in enumerate(header):
            if cell != "":
                data_extent = idx
        header = header[: data_extent + 1]

        # Basic validation to 1) fail early before processing entire file; and
        # 2) make some basic guarantees about the header for things in this
        # class that use the header as part of reading the file.
        column_names = set(header)
        if "" in column_names:
            raise MetadataFileError(
                "Found at least one column without a name in the header. Each "
                "column must be named."
            )
        elif len(header) != len(column_names):
            duplicates = find_duplicates(header)
            raise MetadataFileError(
                "Column names must be unique. The following column names are "
                "duplicated: %s" % (", ".join(repr(e) for e in sorted(duplicates)))
            )

        # Skip the first element of the header because we know it is a valid ID
        # header. The other column names are validated to ensure they *aren't*
        # valid ID headers.
        for column_name in header[1:]:
            if is_id_header(column_name):
                raise MetadataFileError(
                    "Metadata column name %r conflicts with a name reserved "
                    "for the ID column header. Reserved ID column headers:"
                    "\n\n%s" % (column_name, FORMATTED_ID_HEADERS)
                )

        return header

    def _read_directives(self, header):
        directives = {}
        for row in self._reader:
            directive_kind = None

            if not self._is_directive(row):
                self._reader = itertools.chain([row], self._reader)
                break

            if self._is_column_types_directive(row):
                directive_kind = "types"
            elif self._is_missing_directive(row):
                directive_kind = "missing"
            else:
                raise MetadataFileError(
                    "Unrecognized directive %r. Only the #sk:types, #q2:types"
                    " and #sk:missing, #q2:missing directives are supported at this"
                    " time." % row[0]
                )

            if directive_kind in directives:
                raise MetadataFileError(
                    "Found duplicate directive %r. Each directive may "
                    "only be specified a single time." % row[0]
                )

            row = self._match_header_len(row, header)

            collected = {name: arg for name, arg in zip(header[1:], row[1:]) if arg}

            directives[directive_kind] = collected

        if "types" in directives:
            column_types = directives["types"]
            for column_name, column_type in column_types.items():
                type_nocase = column_type.lower()
                if type_nocase in SUPPORTED_COLUMN_TYPES:
                    column_types[column_name] = type_nocase
                else:
                    fmt_column_types = ", ".join(
                        repr(e) for e in sorted(SUPPORTED_COLUMN_TYPES)
                    )
                    raise MetadataFileError(
                        "Column %r has an unrecognized column type %r "
                        "specified in its #sk:types or #q2:types directive. "
                        "Supported column types (case-insensitive): %s"
                        % (column_name, column_type, fmt_column_types)
                    )

        if "missing" in directives:
            for column_name, column_missing in directives["missing"].items():
                if column_missing not in BUILTIN_MISSING:
                    raise MetadataFileError(
                        "Column %r has an unrecognized missing value scheme %r"
                        " specified in its #sk:missing or #q2:missing directive."
                        " Supported missing value schemes (case-sensitive): %s"
                        % (column_name, column_missing, list(BUILTIN_MISSING))
                    )

        return directives

    def _read_data(self, header):
        ids = []
        data = []
        for row in self._reader:
            if self._is_comment(row):
                continue
            elif self._is_empty(row):
                continue
            elif self._is_directive(row):
                raise MetadataFileError(
                    "Found directive %r outside of the directives section of "
                    "the file. Directives may only appear immediately after "
                    "the header." % row[0]
                )
            elif self._is_header(row):
                raise MetadataFileError(
                    "Metadata ID %r conflicts with a name reserved for the ID "
                    "column header. Reserved ID column headers:\n\n%s"
                    % (row[0], FORMATTED_ID_HEADERS)
                )

            row = self._match_header_len(row, header)
            ids.append(row[0])
            data.append(row[1:])
        return ids, data

    def _strip_cell_whitespace(self, row):
        return [cell.strip() for cell in row]

    def _match_header_len(self, row, header):
        row_len = len(row)
        header_len = len(header)

        if row_len < header_len:
            # Pad row with empty cells to match header length.
            row = row + [""] * (header_len - row_len)
        elif row_len > header_len:
            trailing_row = row[header_len:]
            if not self._is_empty(trailing_row):
                raise MetadataFileError(
                    "Metadata row contains more cells than are declared by "
                    "the header. The row has %d cells, while the header "
                    "declares %d cells." % (row_len, header_len)
                )
            row = row[:header_len]
        return row

    def _is_empty(self, row):
        # `all` returns True for an empty iterable, so this check works for a
        # row of zero elements (corresponds to a blank line in the file).
        return all((cell == "" for cell in row))

    def _is_comment(self, row):
        return (
            len(row) > 0
            and row[0].startswith("#")
            and not self._is_directive(row)
            and not self._is_header(row)
        )

    def _is_header(self, row):
        if len(row) == 0:
            return False
        return is_id_header(row[0])

    def _is_directive(self, row):
        return len(row) > 0 and row[0].startswith(("#sk:", "#q2:"))

    def _is_column_types_directive(self, row):
        return len(row) > 0 and (row[0].split(" ")[0] in ["#sk:types", "#q2:types"])

    def _is_missing_directive(self, row):
        return len(row) > 0 and (row[0].split(" ")[0] in ["#sk:missing", "#q2:missing"])

    def _cast_column(self, series, column_types, missing_schemes):
        if series.name in missing_schemes:
            scheme = missing_schemes[series.name]
            series = series_encode_missing(series, scheme)
        if series.name in column_types:
            if column_types[series.name] == "numeric":
                return self._to_numeric(series)
            else:  # 'categorical'
                return self._to_categorical(series)
        else:
            # Infer type
            try:
                return self._to_numeric(series)
            except MetadataFileError:
                return self._to_categorical(series)

    def _to_categorical(self, series):
        # Replace empty strings with `None` to force the series to remain
        # dtype=object (this only matters if the series consists solely of
        # missing data). Replacing with np.nan and casting to dtype=object
        # won't retain the correct dtype in the resulting dataframe
        # (`DataFrame.apply` seems to force series consisting solely of np.nan
        # to dtype=float64, even if dtype=object is specified.
        #
        # To replace a value with `None`, the following invocation of
        # `Series.replace` must be used because `None` is a sentinel:
        #     https://stackoverflow.com/a/17097397/3776794
        return series.replace([""], [None])

    def _to_numeric(self, series):
        series = series.replace("", np.nan)
        is_numeric = series.apply(self._is_numeric)
        if is_numeric.all():
            return pd.to_numeric(series, errors="raise")
        else:
            non_numerics = series[~is_numeric].unique()
            raise MetadataFileError(
                "Cannot convert metadata column %r to numeric. The following "
                "values could not be interpreted as numeric: %s"
                % (series.name, ", ".join(repr(e) for e in sorted(non_numerics)))
            )

    def _is_numeric(self, value):
        return isinstance(value, float) or len(_numeric_regex.findall(value)) == 1


class MetadataWriter:
    """Writer for Metadata."""

    def __init__(self, metadata):
        """Initialize Writer for Metadata."""
        self._metadata = metadata

    def write(self, filepath_or_filehandle):
        """Write metadata object to passed file or filehandle."""
        if isinstance(filepath_or_filehandle, str):
            # Newline settings based on recommendation from csv docs:
            # https://docs.python.org/3/library/csv.html#id3
            # Do NOT write a BOM, hence utf-8 not utf-8-sig
            cm = open(filepath_or_filehandle, "w", newline="", encoding="utf-8")
        else:
            cm = filepath_or_filehandle

        with cm as fh:
            tsv_writer = csv.writer(fh, dialect="excel-tab", strict=True)

            md = self._metadata
            header = [md.id_header]
            # NOTE/TODO: The Metadata files written with this method
            # will always have the directives of type #sk:
            # even if a metadata file with directives of type #q2:
            # has been read. This can be changed in the future
            # however we could also decide to just stick with the sk: types.
            types_directive = ["#sk:types"]
            missing_directive = ["#sk:missing"]

            if isinstance(md, SampleMetadata):
                for name, props in md.columns.items():
                    header.append(name)
                    types_directive.append(props.type)
                    missing_directive.append(props.missing_scheme)
            elif isinstance(md, MetadataColumn):
                header.append(md.name)
                types_directive.append(md.type)
                missing_directive.append(md.missing_scheme)
            else:
                raise NotImplementedError

            tsv_writer.writerow(header)
            tsv_writer.writerow(types_directive)
            if self._non_default_missing(missing_directive):
                tsv_writer.writerow(missing_directive)

            df = md.to_dataframe(encode_missing=True)
            df.fillna("", inplace=True)
            # since `applymap` is going to be deprecated soon
            # and `map` may not work on older versions of pandas
            try:
                mapper_ = df.map
            except AttributeError:
                mapper_ = df.applymap
            df = mapper_(self._format)

            tsv_writer.writerows(df.itertuples(index=True))

    def _non_default_missing(self, missing_directive):
        missing = missing_directive[1:]
        result = False
        for m in missing:
            if m != DEFAULT_MISSING:
                result = True
                break

        return result

    def _format(self, value):
        if isinstance(value, str):
            return value
        elif isinstance(value, float):
            # Use fixed precision or scientific notation as necessary (both are
            # roundtrippable in the metadata file format), with up to 15 digits
            # *total* precision (i.e. before and after the decimal point),
            # rounding if necessary. Trailing zeros or decimal points will not
            # be included in the formatted string (e.g. 42.0 will be formatted
            # as "42"). A precision of 15 digits is used because that is within
            # the 64-bit floating point spec (things get weird after that).
            #
            # Using repr() and str() each have their own predefined precision
            # which varies across Python versions. Using the string formatting
            # presentation types (e.g. %g, %f) without specifying a precision
            # will usually default to 6 digits past the decimal point, which
            # seems a little low.
            #
            # References:
            #
            # - https://stackoverflow.com/a/2440786/3776794
            # - https://stackoverflow.com/a/2440708/3776794
            # - https://docs.python.org/3/library/string.html#
            #       format-specification-mini-language
            # - https://stackoverflow.com/a/20586479/3776794
            # - https://drj11.wordpress.com/2007/07/03/python-poor-printing-
            #       of-floating-point/
            return "{0:.15g}".format(value)
        else:
            raise NotImplementedError


# Credit: https://stackoverflow.com/a/4703508/3776794
_numeric_pattern = r"""
    ^[-+]? # optional sign
    (?:
        (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
        |
        (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
    )
    # followed by optional exponent part if desired
    (?: [Ee] [+-]? \d+ ) ?$
"""

_numeric_regex = re.compile(_numeric_pattern, re.VERBOSE)
