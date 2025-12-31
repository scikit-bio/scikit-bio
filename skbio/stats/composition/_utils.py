# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import inspect

import numpy as np
import pandas as pd


def _check_metadata(metadata, matrix, samples=None):
    """Format metadata for differential abundance analysis.

    Parameters
    ----------
    metadata : dataframe_like
        Metadata table.
    matrix : ndarray of shape (n_samples, n_features)
        Data matrix.
    samples : array_like of shape (n_samples,), optional
        Sample IDs.

    Returns
    -------
    pd.DataFrame
        Validated metadata table.

    Notes
    -----
    This function resembles `_check_grouping`.

    """
    if not isinstance(metadata, pd.DataFrame):
        try:
            metadata = pd.DataFrame(metadata)
        except Exception:
            raise TypeError(
                "Metadata must be a pandas DataFrame, or a data structure that can be "
                "converted into a pandas DataFrame, such as a NumPy structured or rec "
                "array, or a dictionary."
            )

    # match lengths
    if samples is None or isinstance(metadata.index, pd.RangeIndex):
        if matrix.shape[0] != metadata.shape[0]:
            raise ValueError("Sample counts in table and metadata are not consistent.")

    # match sample IDs
    else:
        if not metadata.index.equals(pd.Index(samples)):
            try:
                metadata = metadata.loc[samples]
            except KeyError:
                raise ValueError(
                    "Metadata contains sample IDs that are absent in the table."
                )

    if metadata.isnull().values.any():
        raise ValueError("Cannot handle missing values in metadata.")

    return metadata


def _check_sig_test(test, n_groups=None):
    """Validate significance test.

    Parameters
    ----------
    test : str or callable
        Statistical testing function or its name.
    n_groups : int, optional
        Number of sample groups.

    Returns
    -------
    callable
        Statistical testing function.

    """
    if isinstance(test, str):
        import scipy.stats

        try:
            func = getattr(scipy.stats, test)
        except AttributeError:
            raise ValueError(f'Function "{test}" does not exist under scipy.stats.')
    else:
        if not callable(test):
            raise TypeError("`sig_test` must be a function or a string.")
        func = test
        test = test.__name__

    if n_groups is not None:
        sig = inspect.signature(func)
        param = next(iter(sig.parameters.values()))

        is_multi = param.kind is inspect.Parameter.VAR_POSITIONAL
        if n_groups > 2 and not is_multi:
            raise ValueError(
                f'"{test}" is a two-way statistical test whereas {n_groups} sample '
                "groups were provided."
            )

    return func


def _check_p_adjust(name):
    r"""Construct a p-value correction function based on the method name.

    Parameters
    ----------
    name : str
        The name of the p-value correction method. This should match one of the
        method names available in :func:`statsmodels.stats.multitest.multipletests`.

    Returns
    -------
    callable, optional
        Function to correct p-values.

    """
    if name is None:
        return

    from statsmodels.stats.multitest import multipletests as sm_multipletests

    name_ = name.lower()

    # Original options are kept for backwards compatibility
    if name_ in ("holm", "holm-bonferroni"):
        name_ = "holm"
    if name_ in ("bh", "fdr_bh", "benjamini-hochberg"):
        name_ = "fdr_bh"

    def func(pvals):
        r"""Correct p-values for multiple testing problems.

        Parameters
        ----------
        pvals : ndarray of shape (n_tests,)
            Original p-values.

        Returns
        -------
        qvals : ndarray of shape (n_tests,)
            Corrected p-values.

        """
        try:
            res = sm_multipletests(pvals, alpha=0.05, method=name_)
        except ValueError as e:
            if "method not recognized" in str(e):
                raise ValueError(f'"{name}" is not an available FDR correction method.')
            else:
                raise ValueError(
                    f"Cannot perform FDR correction using the {name} method."
                )
        else:
            return res[1]

    return func


def _check_grouping(grouping, matrix, samples=None):
    """Format grouping for differential abundance analysis.

    Parameters
    ----------
    grouping : 1-D array_like
        Vector indicating the assignment of samples to groups. For example,
        these could be strings or integers denoting which group a sample
        belongs to.
    matrix : ndarray of shape (n_samples, n_features)
        Data matrix.
    samples : array_like of shape (n_samples,), optional
        Sample IDs.

    Returns
    -------
    groups : ndarray of (n_groups,)
        Class names.
    labels : ndarray of (n_samples,)
        Class indices by sample.

    Notes
    -----
    If `grouping` is indexed and `samples` is provided, `grouping` will be filtered and
    reordered to match `samples`. Otherwise, `grouping` and `matrix` must have the same
    length, with the assumption that samples are in the same order.

    """
    # match sample IDs
    if samples is not None and isinstance(grouping, pd.Series):
        try:
            grouping = grouping.loc[samples]
        except KeyError:
            raise ValueError(
                "`table` contains sample IDs that are absent in `grouping`."
            )
        else:
            grouping = grouping.to_numpy()

    # match lengths
    else:
        grouping = np.asarray(grouping)

        if grouping.ndim != 1:
            raise ValueError("`grouping` must be convertible to a 1-D vector.")

        if matrix.shape[0] != grouping.shape[0]:
            raise ValueError(
                "Sample counts in `table` and `grouping` are not consistent."
            )

    # The following code achieves what `pd.isnull` does with NumPy.
    null_errmsg = "Cannot handle missing values in `grouping`."
    if np.isdtype(grouping.dtype, "numeric"):
        if np.isnan(grouping).any():
            raise ValueError(null_errmsg)
    else:
        if (grouping != grouping).any() or np.equal(grouping, None).any():
            raise ValueError(null_errmsg)

    return np.unique(grouping, return_inverse=True)


def _check_trt_ref_groups(treatment, reference, groups, labels):
    """Extract treatment and reference group indices.

    Parameters
    ----------
    treatment : str, int or None
        Treatment group label.
    reference : str, int or None
        Reference group label.
    groups : ndarray of (n_groups,)
        Class names.
    labels : ndarray of (n_samples,)
        Class indices by sample.

    Returns
    -------
    trt_idx : ndarray of (n_samples_in_treatment,)
        Sample indices in the treatment group.
    ref_idx : ndarray of (n_samples_in_reference,)
        Sample indices in the reference group.

    Raises
    ------
    ValueError
        If treatment or reference group is not found.
    ValueError
        If treatment and reference groups are the same.
    ValueError
        If there are less than two groups.

    """
    if len(groups) < 2:
        raise ValueError("There must be at least two groups in grouping.")

    if treatment is not None:
        try:
            (trt_i,) = np.flatnonzero(groups == treatment)
        except ValueError:
            raise ValueError(f"Treatment group {treatment} is not found in grouping.")
    else:
        trt_i = 0
    trt_idx = np.flatnonzero(labels == trt_i)

    if reference is not None:
        try:
            (ref_i,) = np.flatnonzero(groups == reference)
        except ValueError:
            raise ValueError(f"Reference group {reference} is not found in grouping.")
        if trt_i == ref_i:
            raise ValueError("Treatment and reference groups must not be identical.")
        ref_idx = np.flatnonzero(labels == ref_i)
    else:
        ref_idx = np.flatnonzero(labels != trt_i)

    return trt_idx, ref_idx


def _type_cast_to_float(df):
    """Attempt to cast all of the values in dataframe to float.

    This will try to type cast all of the series within the dataframe into floats. If a
    column cannot be type casted, it will be kept as is.

    Parameters
    ----------
    df : pd.DataFrame

    Returns
    -------
    pd.DataFrame

    """
    df = df.copy()
    for col in df.select_dtypes(exclude=["float64"]).columns:
        try:
            df[col] = df[col].astype("float64")
        except (ValueError, TypeError):
            continue
    return df
