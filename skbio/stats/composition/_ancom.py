# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd

from skbio.util._decorator import params_aliased
from skbio.table._tabular import _ingest_table
from ._base import _check_composition
from ._utils import _check_grouping, _check_sig_test, _check_p_adjust


@params_aliased(
    [
        ("p_adjust", "multiple_comparisons_correction", "0.6.0", True),
        ("sig_test", "significance_test", "0.7.0", True),
    ]
)
def ancom(
    table,
    grouping,
    alpha=0.05,
    tau=0.02,
    theta=0.1,
    p_adjust="holm",
    sig_test="f_oneway",
    percentiles=None,
):
    r"""Perform differential abundance test using ANCOM.

    Analysis of composition of microbiomes (ANCOM) is done by calculating
    pairwise log ratios between all features and performing a significance
    test to determine if there is a significant difference in feature ratios
    with respect to the variable of interest.

    In an experiment with only two treatments, this tests the following
    hypothesis for feature :math:`i`:

    .. math::

        H_{0i}: \mathbb{E}[\ln(u_i^{(1)})] = \mathbb{E}[\ln(u_i^{(2)})]

    where :math:`u_i^{(1)}` is the mean abundance for feature :math:`i` in the
    first group and :math:`u_i^{(2)}` is the mean abundance for feature
    :math:`i` in the second group.

    .. versionchanged:: 0.7.0
        Computational efficiency significantly improved.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Matrix of strictly positive values (i.e. counts or proportions). See
        :ref:`supported formats <table_like>`.

        .. note::
            If the table contains zero values, one should add a pseudocount or apply
            :func:`multi_replace` to convert all values into positive numbers.

    grouping : pd.Series or 1-D array_like
        Vector indicating the assignment of samples to groups. These could be strings
        or integers denoting which group a sample belongs to. If it is a pandas Series
        and the table contains sample IDs, its index will be filtered and reordered to
        match the sample IDs. Otherwise, it must be the same length as the samples in
        the table.
    alpha : float, optional
        Significance level for each of the statistical tests. This can can be
        anywhere between 0 and 1 exclusive.
    tau : float, optional
        A constant used to determine an appropriate cutoff. A value close to
        zero indicates a conservative cutoff. This can can be anywhere between
        0 and 1 exclusive.
    theta : float, optional
        Lower bound for the proportion for the *W*-statistic. If all *W*-
        statistics are lower than theta, then no features will be detected to
        be significantly different. This can can be anywhere between 0 and 1
        exclusive.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are
        Holm-Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-Hochberg
        ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported by
        statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.
    sig_test : str or callable, optional
        A function to test for significance between classes. It must be able to
        accept at least two vectors of floats and returns a test statistic and
        a *p*-value. Functions under ``scipy.stats`` can be directly specified
        by name. The default is one-way ANOVA ("f_oneway").

        .. versionchanged:: 0.7.0
            Test function must accept 2-D arrays as input, perform batch testing, and
            return 1-D arrays. SciPy functions have this capability. Custom functions
            may need modification.

        .. versionchanged:: 0.6.0
            Accepts test names in addition to functions.

    percentiles : iterable of floats, optional
        Percentile abundances to return for each feature in each group. By
        default, will return the minimum, 25th percentile, median, 75th
        percentile, and maximum abundances for each feature in each group.

    Returns
    -------
    pd.DataFrame
        A table of features, their *W*-statistics and whether the null
        hypothesis is rejected.

        - ``W``: *W*-statistic, or the number of features that the current
          feature is tested to be significantly different against.

        - ``Signif``: Whether the feature is significantly differentially
          abundant across groups (``True``) or not (``False``).

        .. versionchanged:: 0.7.0
            Renamed ``Reject null hypothesis`` as ``Signif``.

    pd.DataFrame
        A table of features and their percentile abundances in each group. If
        ``percentiles`` is empty, this will be an empty ``pd.DataFrame``. The
        rows in this object will be features, and the columns will be a
        multi-index where the first index is the percentile, and the second
        index is the group.

    See Also
    --------
    ancombc
    multi_replace
    scipy.stats.ttest_ind
    scipy.stats.f_oneway
    scipy.stats.wilcoxon
    scipy.stats.kruskal

    Notes
    -----
    The developers of ANCOM recommend the following significance tests ([1]_,
    Supplementary File 1, top of page 11):

    - If there are two groups, use the standard parametric *t*-test
      (``ttest_ind``) or the non-parametric Mann-Whitney rank test
      (``mannwhitneyu``).

    - For paired samples, use the parametric paired *t*-test (``ttest_rel``) or
      the non-parametric Wilcoxon signed-rank test (``wilcoxon``).

    - If there are more than two groups, use the parametric one-way ANOVA
      (``f_oneway``) or the non-parametric Kruskal-Wallis test (``kruskal``).

    - If there are multiple measurements obtained from the individuals, use a
      Friedman test (``friedmanchisquare``).

    Because one-way ANOVA is equivalent to the standard *t*-test when the
    number of groups is two, we default to ``f_oneway`` here, which can be used
    when there are two or more groups.

    Users should refer to the documentation of these tests in SciPy to
    understand the assumptions made by each test.

    This method cannot handle any zero counts as input, since the logarithm
    of zero cannot be computed.  While this is an unsolved problem, many
    studies, including [1]_, have shown promising results by adding
    pseudocounts to all values in the matrix. In [1]_, a pseudocount of 0.001
    was used, though the authors note that a pseudocount of 1.0 may also be
    useful. Zero counts can also be addressed using the ``multi_replace`` method.

    References
    ----------
    .. [1] Mandal et al. "Analysis of composition of microbiomes: a novel
       method for studying microbial composition", Microbial Ecology in Health
       & Disease, (2015), 26.

    Examples
    --------
    >>> from skbio.stats.composition import ancom
    >>> import pandas as pd

    Let's load in a DataFrame with six samples and seven features (e.g., these
    may be bacterial taxa):

    >>> table = pd.DataFrame([[12, 11, 10, 10, 10, 10, 10],
    ...                       [9,  11, 12, 10, 10, 10, 10],
    ...                       [1,  11, 10, 11, 10, 5,  9],
    ...                       [22, 21, 9,  10, 10, 10, 10],
    ...                       [20, 22, 10, 10, 13, 10, 10],
    ...                       [23, 21, 14, 10, 10, 10, 10]],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'],
    ...                      columns=['b1', 'b2', 'b3', 'b4', 'b5', 'b6',
    ...                               'b7'])

    Then create a grouping vector. In this example, there is a treatment group
    and a placebo group.

    >>> grouping = pd.Series(['treatment', 'treatment', 'treatment',
    ...                       'placebo', 'placebo', 'placebo'],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'])

    Now run ``ancom`` to determine if there are any features that are
    significantly different in abundance between the treatment and the placebo
    groups. The first DataFrame that is returned contains the ANCOM test
    results, and the second contains the percentile abundance data for each
    feature in each group.

    >>> ancom_df, percentile_df = ancom(table, grouping)
    >>> ancom_df['W'] # doctest: +ELLIPSIS
    b1    0
    b2    4
    b3    0
    b4    1
    b5    1
    b6    0
    b7    1
    Name: W, dtype: ...

    The *W*-statistic is the number of features that a single feature is tested
    to be significantly different against. In this scenario, ``b2`` was
    detected to have significantly different abundances compared to four of the
    other features. To summarize the results from the *W*-statistic, let's take
    a look at the results from the hypothesis test. The ``Signif`` column in the
    table indicates whether the null hypothesis was rejected, and that a feature
    was therefore observed to be differentially abundant across the groups.

    >>> ancom_df['Signif']
    b1    False
    b2     True
    b3    False
    b4    False
    b5    False
    b6    False
    b7    False
    Name: Signif, dtype: bool

    From this we can conclude that only ``b2`` was significantly different in
    abundance between the treatment and the placebo. We still don't know, for
    example, in which group ``b2`` was more abundant. We therefore may next be
    interested in comparing the abundance of ``b2`` across the two groups.
    We can do that using the second DataFrame that was returned. Here we
    compare the median (50th percentile) abundance of ``b2`` in the treatment
    and placebo groups:

    >>> percentile_df[50.0].loc['b2']
    Group
    placebo      21.0
    treatment    11.0
    Name: b2, dtype: float64

    We can also look at a full five-number summary for ``b2`` in the treatment
    and placebo groups:

    >>> percentile_df.loc['b2'] # doctest: +NORMALIZE_WHITESPACE
    Percentile  Group
    0.0         placebo      21.0
    25.0        placebo      21.0
    50.0        placebo      21.0
    75.0        placebo      21.5
    100.0       placebo      22.0
    0.0         treatment    11.0
    25.0        treatment    11.0
    50.0        treatment    11.0
    75.0        treatment    11.0
    100.0       treatment    11.0
    Name: b2, dtype: float64

    Taken together, these data tell us that ``b2`` is present in significantly
    higher abundance in the placebo group samples than in the treatment group
    samples.

    """
    matrix, samples, features = _ingest_table(table)

    groups, labels = _check_grouping(grouping, matrix, samples)

    _check_composition(np, matrix, nozero=True)

    # validate parameters
    if not 0 < alpha < 1:
        raise ValueError("`alpha`=%f is not within 0 and 1." % alpha)

    if not 0 < tau < 1:
        raise ValueError("`tau`=%f is not within 0 and 1." % tau)

    if not 0 < theta < 1:
        raise ValueError("`theta`=%f is not within 0 and 1." % theta)

    # validate percentiles
    if percentiles is None:
        percentiles = np.arange(0, 125, 25.0)
    else:
        if not isinstance(percentiles, np.ndarray):
            percentiles = np.fromiter(percentiles, dtype=float)
        if (percentiles < 0.0).any() or (percentiles > 100.0).any():
            raise ValueError("Percentiles must be in the range [0, 100].")
        n_pcts = len(percentiles)
        percentiles = np.unique(percentiles)
        if percentiles.size != n_pcts:
            raise ValueError("Percentile values must be unique.")

    n_groups = len(groups)
    if n_groups == len(labels):
        raise ValueError(
            "All values in `grouping` are unique. This method cannot "
            "operate on a grouping vector with only unique values (e.g., "
            "there are no 'within' variance because each group of samples "
            "contains only a single sample)."
        )
    elif n_groups == 1:
        raise ValueError(
            "All values the `grouping` are the same. This method cannot "
            "operate on a grouping vector with only a single group of samples"
            "(e.g., there are no 'between' variance because there is only a "
            "single group)."
        )

    # validate significance test
    if sig_test is None:
        sig_test = "f_oneway"
    test_f = _check_sig_test(sig_test, n_groups)

    # compare log ratios
    pval_mat = _log_compare(matrix, labels, n_groups, test_f)

    # correct for multiple testing problem
    if p_adjust is not None:
        func = _check_p_adjust(p_adjust)
        pval_mat = np.apply_along_axis(func, 1, pval_mat)

    np.fill_diagonal(pval_mat, 1)

    # calculate W-statistics
    n_feats = matrix.shape[1]
    W = (pval_mat < alpha).sum(axis=1)
    c_start = W.max() / n_feats
    if c_start < theta:
        reject = np.zeros_like(W, dtype=bool)
    else:
        # Select appropriate cutoff
        cutoff = c_start - np.linspace(0.05, 0.25, 5)
        prop_cut = (W[:, None] > n_feats * cutoff).mean(axis=0)
        dels = np.abs(prop_cut - np.roll(prop_cut, -1))
        dels[-1] = 0

        if (dels[0] < tau) and (dels[1] < tau) and (dels[2] < tau):
            nu = cutoff[1]
        elif (dels[0] >= tau) and (dels[1] < tau) and (dels[2] < tau):
            nu = cutoff[2]
        elif (dels[1] >= tau) and (dels[2] < tau) and (dels[3] < tau):
            nu = cutoff[3]
        else:
            nu = cutoff[4]
        reject = W >= nu * n_feats

    ancom_df = pd.DataFrame(
        {
            "W": pd.Series(W, index=features),
            "Signif": pd.Series(reject, index=features),
        }
    )

    # calculate percentiles
    if percentiles.size == 0:
        return ancom_df, pd.DataFrame()
    data = []
    columns = []
    for i, group in enumerate(groups):
        feat_dists = matrix[labels == i]
        for percentile in percentiles:
            columns.append((percentile, group))
            data.append(np.percentile(feat_dists, percentile, axis=0))
    columns = pd.MultiIndex.from_tuples(columns, names=["Percentile", "Group"])
    percentile_df = pd.DataFrame(np.asarray(data).T, columns=columns, index=features)
    return ancom_df, percentile_df


def _log_compare(matrix, labels, n, test):
    """Compare pairwise log ratios between sample groups.

    Calculate pairwise log ratios between all features and perform a statistical test
    to determine if there is a significant difference in feature ratios with respect
    to the variable of interest.

    Parameters
    ----------
    matrix : ndarray of shape (n_samples, n_features)
        Data matrix.
    labels : ndarray of shape (n_samples,)
        Group indices (0-indexed, consecutive).
    n : int
        Number of groups.
    test : callable
        Statistical test to run.

    Returns
    -------
    ndarray of shape (n_features, n_features)
        p-value matrix.

    """
    # note: `n` can be simply computed with `labels.max()`. It is supplied instead to
    # save compute.

    # log-transform data
    log_mat = np.log(matrix)

    # divide data by sample group
    grouped = [log_mat[labels == i] for i in range(n)]

    # determine all pairs of feature indices
    m = matrix.shape[1]
    ii, jj = np.triu_indices(m, k=1)

    # calculate all log ratios (pairwise difference of log values)
    log_ratios = [x[:, ii] - x[:, jj] for x in grouped]

    # run statistical test on the 2-D arrays in a vectorized manner
    _, pvals = test(*log_ratios)

    # populate p-value matrix
    pval_mat = np.empty((m, m))
    pval_mat[ii, jj] = pval_mat[jj, ii] = pvals
    np.fill_diagonal(pval_mat, 0)

    return pval_mat
