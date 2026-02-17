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

# from ._enan import make_nan_with_payload as _make_nan_with_payload
# from ._enan import get_payload_from_nan as _get_payload_from_nan


# def _encode_terms(namespace):
#     enum = _MISSING_ENUMS[namespace]
#     namespace = _NAMESPACE_LOOKUP.index(namespace)

#     def encode(x):
#         if not isinstance(x, str):
#             return x
#         try:
#             code = enum.index(x)
#         except ValueError:
#             return x
#         return _make_nan_with_payload(code, namespace=namespace)

#     return encode


# def _handle_insdc_missing(series):
#     return series.apply(_encode_terms("INSDC:missing"))


# def _handle_blank(series):
#     return series


# def _handle_no_missing(series):
#     if series.isna().any():
#         raise ValueError(
#             "Missing values are not allowed in series/column"
#             " (name=%r) when using scheme 'no-missing'." % series.name
#         )
#     return series


# BUILTIN_MISSING = {
#     "INSDC:missing": _handle_insdc_missing,
#     "blank": _handle_blank,
#     "no-missing": _handle_no_missing,
# }
# _MISSING_ENUMS = {
#     "INSDC:missing": (
#         "not applicable",
#         "missing",
#         "not collected",
#         "not provided",
#         "restricted access",
#     )
# }

# # list index reflects the nan namespace, the "blank"/"no-missing" enums don't
# # apply here, since they aren't actually encoded in the NaNs
# _NAMESPACE_LOOKUP = ["INSDC:missing"]
DEFAULT_MISSING = "blank"

SCHEMES = {
    "INSDC:missing": (
        "not applicable",
        "missing",
        "not collected",
        "not provided",
        "restricted access",
    ),
    "no-missing": None,
    "blank": None,
}


def series_encode_missing(
    series: pd.Series, scheme: str
) -> tuple[pd.Series, pd.Series]:
    """This function preserves information about why a value is missing"""
    if scheme not in SCHEMES:
        raise KeyError("unrecognized NaN scheme.")
    mask = pd.Series(None, index=series.index, dtype="object")
    if scheme == "blank":
        return (series, mask)
    elif scheme == "no-missing":
        if series.isna().any():
            raise ValueError(
                "Missing values are not allowed in series/column"
                " (name=%r) when using scheme 'no-missing'." % series.name
            )
        return (series, mask)
    elif scheme == "INSDC:missing":
        new_series = series.copy()
        for idx, val in series.items():
            if val in SCHEMES[scheme]:
                new_series[idx] = np.nan
                mask[idx] = val
        return (new_series, mask)


def series_extract_missing(series: pd.Series, mask: pd.Series) -> pd.Series:
    """Get missing value reasons from NaNs"""
    # encoded_idx = series.isna() & mask.notna()
    # return mask[encoded_idx]
    res = series[series.isna()].copy()
    res[mask.notna()] = mask[mask.notna()]
    return res
