"""Defines Handling of different Missing classes for Metadata module."""
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import numpy as np

from ._enan import make_nan_with_payload as _make_nan_with_payload
from ._enan import get_payload_from_nan as _get_payload_from_nan


def _encode_terms(namespace):
    enum = _MISSING_ENUMS[namespace]
    namespace = _NAMESPACE_LOOKUP.index(namespace)

    def encode(x):
        if not isinstance(x, str):
            return x
        try:
            code = enum.index(x)
        except ValueError:
            return x
        return _make_nan_with_payload(code, namespace=namespace)

    return encode


def _handle_insdc_missing(series):
    return series.apply(_encode_terms("INSDC:missing"))


def _handle_blank(series):
    return series


def _handle_no_missing(series):
    if series.isna().any():
        raise ValueError(
            "Missing values are not allowed in series/column"
            " (name=%r) when using scheme 'no-missing'." % series.name
        )
    return series


BUILTIN_MISSING = {
    "INSDC:missing": _handle_insdc_missing,
    "blank": _handle_blank,
    "no-missing": _handle_no_missing,
}
_MISSING_ENUMS = {
    "INSDC:missing": (
        "not applicable",
        "missing",
        "not collected",
        "not provided",
        "restricted access",
    )
}

# list index reflects the nan namespace, the "blank"/"no-missing" enums don't
# apply here, since they aren't actually encoded in the NaNs
_NAMESPACE_LOOKUP = ["INSDC:missing"]
DEFAULT_MISSING = "blank"


def series_encode_missing(series: pd.Series, enumeration: str) -> pd.Series:
    """Return encoded Missing values."""
    if not isinstance(enumeration, str):
        TypeError("Wrong type for `enumeration`, expected string")
    try:
        encoder = BUILTIN_MISSING[enumeration]
    except KeyError:
        raise ValueError(
            "Unknown enumeration: %r, (available: %r)"
            % (enumeration, list(BUILTIN_MISSING.keys()))
        )

    new = encoder(series)
    if series.dtype == object and new.isna().all():
        # return to categorical of all missing values
        return new.astype(object)
    return new


def series_extract_missing(series: pd.Series) -> pd.Series:
    """Return extracted Missing types from passed Series."""

    def _decode(x):
        if np.issubdtype(type(x), np.floating) and np.isnan(x):
            code, namespace = _get_payload_from_nan(x)
            if namespace is None:
                return x
            elif namespace == 255:
                raise ValueError("Custom enumerations are not yet supported")
            else:
                try:
                    enum = _MISSING_ENUMS[_NAMESPACE_LOOKUP[namespace]]
                except (IndexError, KeyError):
                    return x

            try:
                return enum[code]
            except IndexError:
                return x

        return x

    missing = series[series.isna()]
    missing = missing.apply(_decode)
    return missing.astype(object)
