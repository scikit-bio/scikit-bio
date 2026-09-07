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


def _adjust_pvalues(pval, method="bh", *, axis=0, n_tests=None, out=None):
    """Perform multiple testing correction of p-values.

    Parameters
    ----------
    pval : array_like of float
        Real p-values in [0, 1], optionally containing NaNs.
    method : str or None, optional
        Method to correct p-values. Options are: Bonferroni ("bonf" or "bonferroni"),
        Holm-Boniferroni ("holm" or "holm-bonferroni"), Benjamini-Hochberg ("bh" or
        "benjamini-hochberg") (default), and Benjamini-Yekutieli ("by" or
        "benjamini-yekutieli"), or any method supported by statsmodels' `multipletests`
        function. Case-insensitive. If None, no correction will be performed.
    axis : int or None, optional
        Axis along which correction will be performed. Each vector on this axis is
        considered as an independent family of p-values. Default is 0. If None, the
        entire array is treated as one family.
    n_tests : int or float, optional
        Effective number of hypotheses in each family, including unobserved tests.
        Must be at least the non-NaN count of every family. If None, use each
        family's non-NaN count. Integer sizes correspond to padding with p-values
        of one. Fractional sizes are supported for Bonferroni, Holm, BH and BY,
        following R's `p.adjust`. Names dispatched to statsmodels require an
        integer size. Ignored when no correction is requested.
    out : ndarray of float, optional
        Location to store the result. Must have the same shape and data type as `pval`.
        Can be `pval` itself for in-place correction. If not provided, a new array will
        be allocated.

    Returns
    -------
    ndarray of float
        Corrected p-values, preserving the shape, floating data type and NaN values of
        the input. Returns `out` if supplied.

    Notes
    -----
    This function corrects p-values independently along an axis within an N-dimensional
    array. The expected use cases involve large feature-by-covariate matrices, in which
    p-values of each covariate are to be corrected independently.

    The algorithmic core is still iterative over individual families of p-values, rather
    than full-array vectorization, which is memory inefficient. Certain calculations are
    performed prior to iteration for re-use. Overall, this function is more efficient
    than statsmodels' `multipletests`.

    NaN p-values are omitted from the calculation. This behavior matches R's `p.adjust`,
    whereas `multipletests` has inconsistent behavior in some methods. This is important
    because ill-conditioned features and covariates often produce non-estimable model
    parameters and p-values. When falling back to `multipletests`, only non-NaN p-values
    are passed to it, ensuring consistent behavior.

    With explicit `n_tests`, the native methods account for unobserved tests
    without padding. BY uses H_floor(n) = digamma(floor(n) + 1) + Euler's
    constant, matching R while avoiding a rank array of length `n_tests`.
    Fallback methods pad with ones.

    """
    # As a future optimization path, batch calculation may be more efficient than
    # per-family iteration. However, because the number of covariates is usually small
    # compared with the number of features, this optimization may not be necessary.

    # TODO: This function casts integer factors into float64 to prevent overflow when
    # the input has a low-precision floating dtype (e.g., float16). This approach is
    # safe with NumPy and CPU, but it needs revision when adopting GPU or Array API
    # backends.

    # TODO: array-like and array API adoption
    pval = np.asarray(pval)
    if not np.issubdtype(pval.dtype, np.floating):
        raise TypeError("`pval` must have a floating-point data type.")
    dtype = pval.dtype

    # TODO: There should be a generic helper validating provided `out`.
    if out is None:
        out = np.empty_like(pval, dtype=dtype)

    # No correction, or an empty testing family; just return input.
    if method is None or n_tests == 0:
        if out is not pval:
            out[...] = pval
        return out

    # Empty input
    if not pval.size:
        return out

    size = pval.size if axis is None else pval.shape[axis]
    total = size if n_tests is None else n_tests

    # Determine built-in method, or fallback to statsmodels
    key = method.lower()
    bonf = key in ("bonf", "bonferroni")
    holm = key in ("holm", "holm-bonferroni")
    bh = key in ("bh", "benjamini-hochberg")
    by = key in ("by", "benjamini-yekutieli")
    if holm:
        factors = total - np.arange(size, dtype=float)
    elif bh or by:
        rank = np.arange(1, size + 1, dtype=float)
        factors = rank / total
        if by:
            if n_tests is None:
                harmonic = np.cumsum(1 / rank)
            else:
                from scipy.special import digamma

                harmonic = digamma(np.floor(n_tests) + 1) + np.euler_gamma
    elif not bonf:
        if n_tests is not None:
            if n_tests != int(n_tests):
                raise ValueError(
                    f"Fractional `n_tests` is not supported for method {method!r}."
                )
            n_tests = int(n_tests)
        func = _sm_p_adjust(method)

    # Determine p-value family axis
    if axis is None:
        slices = ((pval.ravel(), out.flat),)
    else:
        pval_ = np.moveaxis(pval, axis, -1)
        qval_ = np.moveaxis(out, axis, -1)
        slices = ((pval_[idx], qval_[idx]) for idx in np.ndindex(pval_.shape[:-1]))

    for col, dest in slices:
        # Mask NaN values
        valid = ~np.isnan(col)
        if missing := not valid.all():
            values = col[valid]  # Gather before clearing the output to allow out=pval.
            dest[:] = np.nan
            if not values.size:
                continue
            result = values
        else:
            values, result = col, dest

        # Per-method calculation
        n = values.size
        total = n if n_tests is None else n_tests
        if bonf:
            result[:] = np.minimum(values * np.float64(total), 1)
        elif holm or bh or by:
            order = np.argsort(values)
            if holm:
                scale = factors[-n:] if n_tests is None else factors[:n]
                adjusted = values[order] * scale
                np.maximum.accumulate(adjusted, out=adjusted)
                np.minimum(adjusted, 1, out=adjusted)
            else:
                scale = (
                    factors[:n] if n_tests is not None or n == size else rank[:n] / n
                )
                adjusted = values[order] / scale
                np.minimum.accumulate(adjusted[::-1], out=adjusted[::-1])
                if by:
                    adjusted *= harmonic[n - 1] if n_tests is None else harmonic
                if by or n_tests is not None:
                    np.minimum(adjusted, 1, out=adjusted)
            result[order] = adjusted
        else:
            if total > n:
                values = np.pad(values, (0, total - n), constant_values=1.0)
            result[:] = func(values)[:n]

        if missing:
            dest[valid] = result
    return out


def _sm_p_adjust(name):
    r"""Import a p-value correction method from statsmodels.

    Parameters
    ----------
    name : str
        The name of the p-value correction method. This should match one of the
        method names supported by statsmodels' `multipletests`.

    Returns
    -------
    callable, optional
        Function to correct p-values.

    """
    if name is None:
        return
    method = name.lower()

    # TODO: Make statsmodels an optional dependency
    from statsmodels.stats.multitest import multipletests as sm_multipletests

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
            res = sm_multipletests(pvals, alpha=0.05, method=method)
        except ValueError as e:
            if "method not recognized" in str(e):
                raise ValueError(
                    f"'{name}' is not an available multiple testing correction method "
                    "supported by scikit-bio or statsmodels."
                )
            else:  # pragma: no cover
                raise ValueError(
                    f"Cannot perform multiple testing correction using the {name} "
                    "method."
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


def _build_dmatrix(formula, metadata, dtype=None):
    """Build a design matrix from sample metadata using Patsy.

    Metadata values are passed to Patsy unchanged. Therefore, Patsy infers whether
    factors are numerical or categorical from their existing data types and the
    formula specification.

    Parameters
    ----------
    formula : str or generic Formula object
        Formula defining the model.
    metadata : pd.DataFrame
        Sample metadata containing factors referenced in `formula`.
    dtype : dtype_like, optional
        Data type of the resulting design matrix. If omitted, use Patsy's
        default output data type.

    Returns
    -------
    patsy.DesignMatrix
        Design matrix generated from the formula and metadata.

    Notes
    -----
    TODO: patsy has been superseded by formulaic. Consider a replacement.

    """
    if dtype is None:
        from patsy import dmatrix

        return dmatrix(formula, metadata, eval_env=1, return_type="matrix")

    from patsy import build_design_matrices, incr_dbuilder

    def data_iter_maker():
        return iter([metadata])

    design_info = incr_dbuilder(formula, data_iter_maker, eval_env=1)
    return build_design_matrices(
        [design_info], metadata, return_type="matrix", dtype=dtype
    )[0]
