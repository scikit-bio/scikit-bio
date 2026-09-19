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
        raise ValueError(f"Unrecognized missing value scheme: {scheme!r}")
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
    # Get indices where series has NaN
    na_idx = series.isna()
    mask_values = mask[na_idx]
    has_encoded_terms = mask_values.notna().any()

    if has_encoded_terms:
        # Start with the original series values (preserves None vs nan distinction)
        res = series[na_idx].copy()
        # Ensure object dtype to hold both NaN/None and string values from the mask
        if res.dtype != object:
            res = res.astype(object)
        # Fill in the mask values where they exist
        res[mask_values.notna()] = mask_values[mask_values.notna()]
    else:
        # No encoded terms, just return the NaN values with original dtype
        res = series[na_idx].copy()

    return res
