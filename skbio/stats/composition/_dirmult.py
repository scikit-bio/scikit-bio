# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn, catch_warnings, simplefilter

import numpy as np
import pandas as pd

from skbio.util import get_rng
from skbio.table._tabular import _ingest_table
from ._base import _check_composition, clr
from ._utils import (
    _check_metadata,
    _check_grouping,
    _check_trt_ref_groups,
    _check_p_adjust,
    _type_cast_to_float,
)


def _dirmult_draw(matrix, rng):
    """Resample data from a Dirichlet-multinomial posterior distribution.

    See Also
    --------
    numpy.random.Generator.gamma
    numpy.random.Generator.dirichlet

    Notes
    -----
    This function uses a Gamma distribution to replace a Dirichlet distribution. The
    result is precisely identical to (and reproducible given the same seed):

    .. code-block:: python
       return clr(np.apply_along_axis(rng.dirichlet, axis=1, arr=matrix))

    A Dirichlet distribution is essentially a standard Gamma distribution normalized
    by row sums. Meanwhile, CLR is independent of scale, therefore the normalization
    step can be omitted.

    `gamma` can vectorize to a 2-D array whereas `dirichlet` cannot.

    """
    return clr(rng.gamma(shape=matrix, scale=1.0, size=matrix.shape), validate=False)


def dirmult_ttest(
    table,
    grouping,
    treatment=None,
    reference=None,
    pseudocount=0.5,
    draws=128,
    p_adjust="holm",
    seed=None,
):
    r"""Perform *t*-test using Dirichlet-multinomial distribution.

    The Dirichlet-multinomial distribution is a compound distribution that
    combines a Dirichlet distribution over the probabilities of a multinomial
    distribution. This distribution is used to model the distribution of
    species abundances in a community.

    To perform the *t*-test, we first fit a Dirichlet-multinomial distribution
    for each sample, and then we compute the fold change and *p*-value for each
    feature. The fold change is computed as the difference between the
    samples of the two groups. *t*-tests are then performed on the posterior
    samples, drawn from each Dirichlet-multinomial distribution. The
    log-fold changes as well as their credible intervals, the *p*-values and
    the multiple comparison corrected *p*-values are reported.

    This process mirrors the approach performed by the R package "ALDEx2" [1]_.

    Additionally, this function excludes hits with a 95% confidence interval of
    fold-change crossing zero during any draw. This step further reduces false
    positive hits, especially among low-abundance features.

    .. versionchanged:: 0.7.0
        Computational efficiency significantly improved.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing count or proportional abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    grouping : pd.Series or 1-D array_like
        Vector indicating the assignment of samples to groups. These could be strings
        or integers denoting which group a sample belongs to. If it is a pandas Series
        and the table contains sample IDs, its index will be filtered and reordered to
        match the sample IDs. Otherwise, it must be the same length as the samples in
        the table.
    treatment : str, optional
        Name of the treatment group. The *t*-test is computed between the ``treatment``
        group and the ``reference`` group specified in the ``grouping`` vector. If
        omitted, the first group in the sorted order of all group names will be the
        treatment group.
    reference : str, optional
        Name of the reference group. See above. If omitted, all groups other than the
        treatment group will be combined as the reference group.

        .. versionchanged:: 0.7.0
            ``treatment`` and ``reference`` are now optional.

    pseudocount : float, optional
        A non-zero value added to the input counts to ensure that all of the estimated
        abundances are strictly greater than zero. Default is 0.5.
    draws : int, optional
        The number of draws from the Dirichilet-multinomial posterior distribution
        More draws provide higher uncertainty surrounding the estimated
        log-fold changes and *p*-values.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are
        Holm-Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-Hochberg
        ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported by
        statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance for drawing from the
        Dirichlet distribution. See :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    pd.DataFrame
        A table of features, their log-fold changes and other relevant statistics.

        - ``T-statistic``: *t*-statistic of Welch's *t*-test. The reported value is the
          average across all of the posterior draws.

        - ``Log2(FC)``: Expected log2-fold change of abundance from the reference group
          to the treatment group. The value is expressed in the center log ratio (see
          :func:`clr`) transformed coordinates. The reported value is the average of
          all of the log2-fold changes computed from each of the posterior draws.

        - ``CI(2.5)``: 2.5% quantile of the log2-fold change. The reported value is the
          minimum of all of the 2.5% quantiles computed from each of the posterior
          draws.

        - ``CI(97.5)``: 97.5% quantile of the log2-fold change. The reported value is
          the maximum of all of the 97.5% quantiles computed from each of the posterior
          draws.

        - ``pvalue``: *p*-value of Welch's *t*-test. The reported value is the average
          of all of the *p*-values computed from each of the posterior draws.

        - ``qvalue``: Corrected *p*-value of Welch's *t*-test for multiple comparisons.
          The reported value is the average of all of the *q*-values computed from each
          of the posterior draws.

        - ``Signif``: Whether feature is significantly differentially abundant between
          the treatment and reference groups. A feature marked as "True" suffice: 1)
          The *q*-value must be less than or equal to the significance level (0.05). 2)
          The confidence interval (CI(2.5)..CI(97.5)) must not overlap with zero.

        .. versionchanged:: 0.7.0
            ``df`` (degrees of freedom) was removed from the report, as this metric is
            inconsistent across draws.

        .. versionchanged:: 0.7.0
            Renamed ``T statistic`` as ``T-statistic``.
            Renamed ``Reject null hypothesis`` as ``Signif``.

    See Also
    --------
    dirmult_ttest
    scipy.stats.ttest_ind
    statsmodels.stats.weightstats.CompareMeans

    Notes
    -----
    The confidence intervals are computed using the mininum 2.5% and maximum
    97.5% bounds computed across all of the posterior draws.

    The reference frame here is the geometric mean. Extracting absolute log
    fold changes from this test assumes that the average feature abundance
    between the ``treatment`` and the ``reference`` groups are the same. If this
    assumption is violated, then the log-fold changes will be biased, and the
    *p*-values will not be reliable. However, the bias is the same across each
    feature, as a result the ordering of the log-fold changes can still be useful.

    One benefit of using the Dirichlet-multinomial distribution is that the
    statistical power increases with regards to the abundance magnitude. More counts
    per sample will shrink the size of the confidence intervals, and can result in
    lower *p*-values.

    References
    ----------
    .. [1] Fernandes, A. D., Reid, J. N., Macklaim, J. M., McMurrough, T. A., Edgell,
       D. R., & Gloor, G. B. (2014). Unifying the analysis of high-throughput
       sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and
       selective growth experiments by compositional data analysis. Microbiome, 2,
       1-13.

    Examples
    --------
    >>> import pandas as pd
    >>> from skbio.stats.composition import dirmult_ttest
    >>> table = pd.DataFrame([[20,  110, 100, 101, 100, 103, 104],
    ...                       [33,  110, 120, 100, 101, 100, 102],
    ...                       [12,  110, 100, 110, 100, 50,  90],
    ...                       [202, 201, 9,  10, 10, 11, 11],
    ...                       [200, 202, 10, 10, 13, 10, 10],
    ...                       [203, 201, 14, 10, 10, 13, 12]],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'],
    ...                      columns=['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7'])
    >>> grouping = pd.Series(['treatment', 'treatment', 'treatment',
    ...                       'placebo', 'placebo', 'placebo'],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'])
    >>> result = dirmult_ttest(table, grouping, 'treatment', 'placebo', seed=0)
    >>> result
        T-statistic  Log2(FC)   CI(2.5)  CI(97.5)    pvalue    qvalue  Signif
    b1   -17.178600 -4.991987 -7.884498 -2.293463  0.003355  0.020131    True
    b2   -16.873187 -2.533729 -3.594590 -1.462339  0.001064  0.007446    True
    b3     6.942727  1.627677 -1.048219  4.750792  0.021130  0.068310   False
    b4     6.522786  1.707221 -0.467481  4.164998  0.013123  0.065613   False
    b5     6.654142  1.528243 -1.036910  3.978387  0.019360  0.068310   False
    b6     3.839520  1.182343 -0.702656  3.556061  0.045376  0.068310   False
    b7     7.600734  1.480232 -0.601277  4.043888  0.017077  0.068310   False

    """
    from statsmodels.stats.weightstats import CompareMeans

    rng = get_rng(seed)

    matrix, samples, features = _ingest_table(table)
    _check_composition(np, matrix)

    # handle zero values
    if pseudocount:
        matrix = matrix + pseudocount

    # get sample indices of treatment and reference groups
    groups, labels = _check_grouping(grouping, matrix, samples)
    trt_idx, ref_idx = _check_trt_ref_groups(treatment, reference, groups, labels)

    cm_params = dict(alternative="two-sided", usevar="unequal")

    # initiate results
    m = matrix.shape[1]
    delta = np.zeros(m)  # inter-group difference
    tstat = np.zeros(m)  # t-test statistic
    pval = np.zeros(m)  # t-test p-value
    lower = np.full(m, np.inf)  # 2.5% percentile of distribution
    upper = np.full(m, -np.inf)  # 97.5% percentile of distribution

    for i in range(draws):
        # Resample data in a Dirichlet-multinomial distribution.
        dir_mat = _dirmult_draw(matrix, rng)

        # Stratify data by group (treatment vs. reference).
        trt_mat = dir_mat[trt_idx]
        ref_mat = dir_mat[ref_idx]

        # Calculate the difference between the two means.
        delta += trt_mat.mean(axis=0) - ref_mat.mean(axis=0)

        # Create a CompareMeans object for statistical testing.
        # Welch's t-test is also available in SciPy's `ttest_ind` (with `equal_var=
        # False`). The current code uses statsmodels' `CompareMeans` instead because
        # it additionally returns confidence intervals.
        cm = CompareMeans.from_data(trt_mat, ref_mat)

        # Perform Welch's t-test to assess the significance of difference.
        tstat_, pval_, _ = cm.ttest_ind(value=0, **cm_params)
        tstat += tstat_
        pval += pval_

        # Calculate confidence intervals.
        # The final lower and upper bounds are the minimum and maximum of all lower
        # and upper bounds seen during sampling, respectively.
        lower_, upper_ = cm.tconfint_diff(alpha=0.05, **cm_params)
        np.minimum(lower, lower_, out=lower)
        np.maximum(upper, upper_, out=upper)

    # Normalize metrics to averages over all replicates.
    delta /= draws
    tstat /= draws
    pval /= draws

    # Correct p-values for multiple comparison.
    if p_adjust is not None:
        qval = _check_p_adjust(p_adjust)(pval)
    else:
        qval = pval
    reject = qval <= 0.05

    # Test if confidence interval includes 0.
    # A significant result (i.e., reject null hypothesis) must simultaneously suffice:
    # 1) q-value <= significance level. 2) confidence interval doesn't include 0.
    # This test is in addition to the original ALDEx2 method. It helps to reduce false
    # positive discoveries of low abundance.
    outer = ((lower > 0) & (upper > 0)) | ((lower < 0) & (upper < 0))
    reject &= outer

    # Convert all log fold changes to base 2.
    log2_ = np.log(2)
    delta /= log2_
    upper /= log2_
    lower /= log2_

    # construct report
    res = pd.DataFrame.from_dict(
        {
            "T-statistic": tstat,
            "Log2(FC)": delta,
            "CI(2.5)": lower,
            "CI(97.5)": upper,
            "pvalue": pval,
            "qvalue": qval,
            "Signif": reject,
        }
    )
    if features is not None:
        res.index = features

    return res


def dirmult_lme(
    table,
    metadata,
    formula,
    grouping,
    pseudocount=0.5,
    draws=128,
    p_adjust="holm",
    seed=None,
    re_formula=None,
    vc_formula=None,
    model_kwargs={},
    fit_method=None,
    fit_converge=False,
    fit_warnings=False,
    fit_kwargs={},
):
    r"""Fit a Dirichlet-multinomial linear mixed effects model.

    .. versionadded:: 0.7.0

    The Dirichlet-multinomial distribution is a compound distribution that
    combines a Dirichlet distribution over the probabilities of a multinomial
    distribution. This distribution is used to model the distribution of
    species abundances in a community.

    To fit the linear mixed effects model we first fit a Dirichlet-multinomial
    distribution for each sample, and then we compute the fold change and
    *p*-value for each feature. The fold change is computed as the slopes
    from the resulting model. Statistical tests are then performed on the posterior
    samples, drawn from each Dirichlet-multinomial distribution. The
    log-fold changes as well as their credible intervals, the *p*-values and
    the multiple comparison corrected *p*-values are reported.

    This function uses the :class:`~statsmodels.regression.mixed_linear_model.MixedLM`
    class from statsmodels.

    .. note::
        Because the analysis iteratively runs many numeric optimizations, it can take
        longer than usual to finish. Please allow extra time for completion.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing count or proportional abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    metadata : pd.DataFrame or 2-D array_like
        The metadata for the model. Rows correspond to samples and columns correspond
        to covariates in the model. Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    formula : str or generic Formula object
        The formula defining the model. Refer to `Patsy's documentation
        <https://patsy.readthedocs.io/en/latest/formulas.html>`_ on how to specify
        a formula.
    grouping : str, pd.Series or 1-D array_like
        A vector or a metadata column name indicating the assignment of samples to
        groups. Samples are independent between groups during model fitting.
    pseudocount : float, optional
        A non-zero value added to the input counts to ensure that all of the
        estimated abundances are strictly greater than zero. Default is 0.5.
    draws : int, optional
        Number of draws from the Dirichlet-multinomial posterior distribution.
        Default is 128.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are
        Holm-Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-Hochberg
        ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported by
        statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance for drawing from the
        Dirichlet distribution. See :func:`details <skbio.util.get_rng>`.
    re_formula : str, optional
        Random coefficient formula. See :meth:`MixedLM.from_formula
        <statsmodels.regression.mixed_linear_model.MixedLM.from_formula>` for details.
    vc_formula : str, optional
        Variance component formula. See ``MixedLM.from_formula`` for details.
    model_kwargs : dict, optional
        Additional keyword arguments to pass to ``MixedLM``.
    fit_method : str or list of str, optional
        Optimization method for model fitting. Can be a single method name, or a list
        of method names to be tried sequentially. See `statsmodels.optimization
        <https://www.statsmodels.org/stable/optimization.html>`_
        for available methods. If None, a default list of methods will be tried.
    fit_converge : bool, optional
        If True, model fittings that were completed but did not converge will be
        excluded from the calculation of final statistics. Default is False.
    fit_warnings : bool, optional
        Issue warnings if any during the model fitting process. Default is False.
        Warnings are usually issued when the optimization methods do not converge,
        which is common in the analysis. Default is False.
    fit_kwargs : dict, optional
        Additional keyword arguments to pass to :meth:`MixedLM.fit
        <statsmodels.regression.mixed_linear_model.MixedLM.fit>`.

    Returns
    -------
    pd.DataFrame
        A table of features and covariates, their log-fold changes and other relevant
        statistics.

        - ``FeatureID``: Feature identifier, i.e., dependent variable.

        - ``Covariate``: Covariate name, i.e., independent variable.

        - ``Reps``: Number of Dirichlet-multinomial posterior draws that supported the
          reported statistics, i.e., the number of successful model fittings on this
          feature. Max: ``draws`` (if none failed). Min: 0 (in which case all
          statistics are NaN).

        - ``Log2(FC)``: Expected log2-fold change of abundance from the reference
          category to the covariate category defined in the formula. The value is
          expressed in the center log ratio (see :func:`clr`) transformed coordinates.
          The reported value is the average of all of the log2-fold changes computed
          from each of the posterior draws.

        - ``CI(2.5)``: 2.5% quantile of the log2-fold change. The reported value is the
          minimum of all of the 2.5% quantiles computed from each of the posterior
          draws.

        - ``CI(97.5)``: 97.5% quantile of the log2-fold change. The reported value is
          the maximum of all of the 97.5% quantiles computed from each of the posterior
          draws.

        - ``pvalue``: *p*-value of the linear mixed effects model. The reported value
          is the average of all of the *p*-values computed from each of the posterior
          draws.

        - ``qvalue``: Corrected *p*-value of the linear mixed effects model for multiple
          comparisons. The reported value is the average of all of the *q*-values
          computed from each of the posterior draws.

        - ``Signif``: Whether the covariate category is significantly differentially
          abundant from the reference category. A feature-covariate pair marked as
          "True" suffice: 1) The *q*-value must be less than or equal to the
          significance level (0.05). 2) The confidence interval (CI(2.5)..CI(97.5))
          must not overlap with zero.

    See Also
    --------
    dirmult_ttest
    statsmodels.formula.api.mixedlm
    statsmodels.regression.mixed_linear_model.MixedLM

    Examples
    --------
    >>> import pandas as pd
    >>> from skbio.stats.composition import dirmult_lme
    >>> table = pd.DataFrame(
    ...     [[1.00000053, 6.09924644],
    ...      [0.99999843, 7.0000045],
    ...      [1.09999884, 8.08474053],
    ...      [1.09999758, 1.10000349],
    ...      [0.99999902, 2.00000027],
    ...      [1.09999862, 2.99998318],
    ...      [1.00000084, 2.10001257],
    ...      [0.9999991, 3.09998418],
    ...      [0.99999899, 3.9999742],
    ...      [1.10000124, 5.0001796],
    ...      [1.00000053, 6.09924644],
    ...      [1.10000173, 6.99693644]],
    ...     index=['u1', 'u2', 'u3', 'x1', 'x2', 'x3', 'y1', 'y2', 'y3', 'z1',
    ...            'z2', 'z3'],
    ...     columns=['Y1', 'Y2'])
    >>> metadata = pd.DataFrame(
    ...     {'patient': [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
    ...      'treatment': [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
    ...      'time': [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],},
    ...     index=['x1', 'x2', 'x3', 'y1', 'y2', 'y3', 'z1', 'z2', 'z3', 'u1',
    ...            'u2', 'u3'])
    >>> result = dirmult_lme(table, metadata, formula='time + treatment',
    ...                      grouping='patient', seed=0, p_adjust='sidak')
    >>> result
      FeatureID  Covariate  Reps  Log2(FC)   CI(2.5)  CI(97.5)    pvalue  \
    0        Y1       time   128 -0.210769 -1.532255  1.122148  0.403737
    1        Y1  treatment   128 -0.744061 -3.401978  1.581917  0.252057
    2        Y2       time   128  0.210769 -1.122148  1.532255  0.403737
    3        Y2  treatment   128  0.744061 -1.581917  3.401978  0.252057
    <BLANKLINE>
         qvalue  Signif
    0  0.644470   False
    1  0.440581   False
    2  0.644470   False
    3  0.440581   False

    """
    from patsy import dmatrix
    from scipy.optimize import OptimizeWarning
    from statsmodels.regression.mixed_linear_model import MixedLM, VCSpec
    from statsmodels.tools.sm_exceptions import ConvergenceWarning

    rng = get_rng(seed)

    matrix, samples, features = _ingest_table(table)

    _check_composition(np, matrix)

    n_feats = matrix.shape[1]
    if n_feats < 2:
        raise ValueError("Table must have at least two features.")
    if features is None:
        features = np.arange(n_feats)

    # validate metadata
    metadata = _check_metadata(metadata, matrix, samples)

    # cast metadata to numbers where applicable
    metadata = _type_cast_to_float(metadata)

    # Instead of directly calling `MixedLM.from_formula` on merged table + metadata,
    # the following code converts metadata into a design matrix based on the formula
    # (as well as re_formula and vc_formula, if applicable), and calls `MixedLM`.
    # This is because the design matrix (independent variable) is always the same
    # whereas the table (dependent variable) is resampled in every replicate. Fixing
    # the design matrix can save conversion overheads.

    # Create a design matrix based on metadata and formula.
    dmat = dmatrix(formula, metadata, return_type="matrix")

    # Obtain the list of covariates by selecting the relevant columns
    covars = dmat.design_info.column_names

    # Remove intercept since it is not a covariate, and is included by default.
    # Then determine the range of rows to be extracted from the model fitting result.
    # (See also `result.model.k_fe`, number of fixed effects.)
    if covars[0] == "Intercept":
        covars = covars[1:]
        n_covars = len(covars)
        covar_range = slice(1, n_covars + 1)
    else:
        n_covars = len(covars)
        covar_range = slice(0, n_covars)

    exog_mat = np.asarray(dmat)

    # parse grouping
    if isinstance(grouping, str):
        try:
            grouping = metadata[grouping].to_numpy()
        except KeyError:
            raise ValueError("Grouping is not a column in the metadata.")
        uniq, grouping = np.unique(grouping, return_inverse=True)
    else:
        uniq, grouping = _check_grouping(grouping, matrix, samples)
    n_groups = len(uniq)

    # random effects matrix
    if re_formula is not None:
        exog_re = np.asarray(dmatrix(re_formula, metadata, return_type="matrix"))
    else:
        exog_re = None

    # variance component matrices
    # see: https://www.statsmodels.org/v0.12.2/examples/notebooks/generated/
    # variance_components.html
    if vc_formula is not None:
        metas = [metadata.iloc[grouping == x] for x in range(n_groups)]
        names, cols, mats = [], [], []
        for name, formula in vc_formula.items():
            names.append(name)
            dmats = [dmatrix(formula, x, return_type="matrix") for x in metas]
            cols.append([x.design_info.column_names for x in dmats])
            mats.append([np.asarray(x) for x in dmats])
        exog_vc = VCSpec(names, cols, mats)
    else:
        exog_vc = None

    # handle zero values
    if pseudocount:
        matrix = matrix + pseudocount

    # initiate results
    shape = (n_feats, n_covars)
    coef = np.zeros(shape)  # coefficient (fold change)
    pval = np.zeros(shape)  # p-value
    lower = np.full(shape, np.inf)  # 2.5% CI
    upper = np.full(shape, -np.inf)  # 97.5% CI

    # number of replicates (draws) LME fitting is successful for each feature
    fitted = np.zeros(n_feats, dtype=int)

    fit_fail_msg = "LME fit failed for feature {} in replicate {}, outputting NaNs."

    with catch_warnings():
        if not fit_warnings:
            simplefilter("ignore", UserWarning)
            simplefilter("ignore", ConvergenceWarning)
            simplefilter("ignore", OptimizeWarning)
            # This is temporary because statsmodels calls scipy in a deprecated way as
            # of v0.14.4.
            simplefilter("ignore", DeprecationWarning)

        for i in range(draws):
            # Resample data in a Dirichlet-multinomial distribution.
            dir_mat = _dirmult_draw(matrix, rng)

            # Fit a linear mixed effects (LME) model for each feature.
            for j in range(n_feats):
                model = MixedLM(
                    dir_mat[:, j],
                    exog_mat,
                    grouping,
                    exog_re=exog_re,
                    exog_vc=exog_vc,
                    **model_kwargs,
                )

                # model fitting (computationally expensive)
                try:
                    result = model.fit(method=fit_method, **fit_kwargs)

                # There are many ways model fitting may fail. Examples are LinAlgError,
                # RuntimeError, OverflowError, and ZeroDivisionError. If any error
                # occurs, the function will still proceed but the current run will be
                # discarded.
                except Exception:
                    warn(fit_fail_msg.format(features[j], i), UserWarning)
                    continue

                # It is common that model fitting successfully finished (no error) but
                # the optimizer did not converge, making the calculated statistics less
                # reliable. The `fit_converge` flag can discard these runs.
                if fit_converge and not result.converged:
                    warn(fit_fail_msg.format(features[j], i), UserWarning)
                    continue

                # update results
                coef[j] += result.params[covar_range]
                pval[j] += result.pvalues[covar_range]

                # calculate confidence interval and update results
                ci = result.conf_int()
                np.minimum(lower[j], ci[covar_range, 0], out=lower[j])
                np.maximum(upper[j], ci[covar_range, 1], out=upper[j])

                fitted[j] += 1

    # deal with fitting failures
    all_fail_msg = "LME fit failed for {} features in all replicates, reporting NaNs."
    mask = fitted > 0
    n_failed = n_feats - mask.sum()
    # all succeeded
    if n_failed == 0:
        mask = slice(None)
    # some failed
    elif n_failed < n_feats:
        warn(all_fail_msg.format(n_failed), UserWarning)
        for x in (coef, pval, lower, upper):
            x[~mask] = np.nan
    # all failed
    else:
        raise ValueError("LME fit failed for all features in all replicates.")

    # normalize metrics to averages over all successful replicates
    fitted_ = fitted.reshape(-1, 1)[mask]
    for x in (coef, pval):
        x[mask] /= fitted_

    # convert all log fold changes to base 2
    log2_ = np.log(2)
    for x in (coef, lower, upper):
        x[mask] /= log2_

    # correct p-values for multiple comparison
    # (only valid replicates are included)
    if p_adjust is not None:
        func = _check_p_adjust(p_adjust)
        qval = np.full(shape, np.nan)
        qval[mask] = np.apply_along_axis(func, 1, pval[mask])
    else:
        qval = pval

    # get significant results (q-value <= 0.05 and CI doesn't cross 0)
    # see `dirmult_ttest`
    reject = np.full(shape, np.nan)
    ii = np.where(mask)[0] if n_failed else np.arange(n_feats)
    for i in ii:
        outer = ((lower[i] > 0) & (upper[i] > 0)) | ((lower[i] < 0) & (upper[i] < 0))
        reject[i] = (qval[i] <= 0.05) & outer

    # construct report
    res = pd.DataFrame.from_dict(
        {
            "FeatureID": [x for x in features for _ in range(n_covars)],
            "Covariate": list(covars) * n_feats,
            "Reps": np.repeat(fitted, n_covars),
            "Log2(FC)": coef.ravel(),
            "CI(2.5)": lower.ravel(),
            "CI(97.5)": upper.ravel(),
            "pvalue": pval.ravel(),
            "qvalue": qval.ravel(),
        }
    )
    # pandas' nullable boolean type
    res["Signif"] = pd.Series(reject.ravel(), dtype="boolean")
    return res
