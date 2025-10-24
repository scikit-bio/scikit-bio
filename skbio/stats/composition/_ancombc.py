# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# This implementation of ANCOM-BC is based on an analysis of the source code
# from the R package ANCOMBC:
# - https://github.com/FrederickHuangLin/ANCOMBC
#
# Which is licensed under Artistic-2.0:
# - https://www.bioconductor.org/packages/release/bioc/html/ANCOMBC.html
#
# We thank Dr. Huang Lin (@FrederickHuangLin) for his helpful advice.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from scipy.stats import norm, chi2
from scipy.optimize import minimize
from patsy import dmatrix

from skbio.table._tabular import _ingest_table
from ._base import _check_composition
from ._utils import _check_metadata, _check_p_adjust, _type_cast_to_float


def ancombc(
    table,
    metadata,
    formula,
    grouping=None,
    max_iter=100,
    tol=1e-5,
    alpha=0.05,
    p_adjust="holm",
):
    r"""Perform differential abundance test using ANCOM-BC.

    Analysis of compositions of microbiomes with bias correction (ANCOM-BC) [1]_ is a
    differential abundance testing method featuring the estimation and correction for
    the bias of differential sampling fractions.

    .. versionadded:: 0.7.1

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
    grouping : str, optional
        A metadata column name of interests for global test, which is used to identify
        the features that are differentially abundant between at least two groups across
        three or more groups in that attribute. The group must be one of the factors in
        the formula. Default is None to skip global test.
    max_iter : int, optional
        Maximum number of iterations for the bias estimation process. Default is 100.
    tol : float, optional
        Absolute convergence tolerance for the bias estimation process. Default is
        1e-5.
    alpha : float, optional
        Significance level for the statistical tests. Must be in the range of (0, 1).
        Default is 0.05.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are Holm-
        Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-
        Hochberg ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported
        by statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.

    Returns
    -------
    res : pd.DataFrame
        A table of features and covariates, their log-fold changes and other relevant
        statistics of the primary ANCOM-BC analysis.

        - ``FeatureID``: Feature identifier, i.e., dependent variable.

        - ``Covariate``: Covariate name, i.e., independent variable.

        - ``Log2(FC)``: Expected log2-fold change of abundance from the reference
          category to the covariate category defined in the formula. The value is
          expressed in the center log ratio (see :func:`clr`) transformed coordinates.

        - ``SE``: Standard error of the estimated Log2(FC).

        - ``W``: *W*-statistic, or the number of features that the current feature is
          tested to be significantly different against.

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

    res_global : pd.DataFrame
        A table of features and statistics of ANCOM-BC global test.

        - ``FeatureID``: Feature identifier, i.e., dependent variable.

        - ``W``: *W*-statistic of Chi-squared statistic that quantifies the overall
          evidence against null hypothesis (the mean abundance of the feature is the
          same across all groups).

        - ``pvalue``: *p*-value of the *W*-statistic in global test.

        - ``qvalue``: Corrected *p*-value in global test.

        - ``Signif``: Whether at least one group mean is different from other
          in grouping variable.

    See Also
    --------
    ancom
    struc_zero
    multi_replace

    Notes
    -----
    The input data table for ANCOM-BC must contain only positive numbers. One needs to
    remove zero values by, e.g., adding a pseudocount of 1.0.
    The categorical data in metadata is sorted alphabetically by default. The reference
    level for each attribute need to be changed manually.

    References
    ----------
    .. [1] Lin, H. and Peddada, S.D., 2020. Analysis of compositions of microbiomes
       with bias correction. Nature communications, 11(1), p.3514.

    Examples
    --------
    >>> from skbio.stats.composition import ancombc
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

    >>> metadata = pd.DataFrame(
    ...     {'grouping': ['treatment', 'treatment', 'treatment',
    ...                   'placebo', 'placebo', 'placebo']},
    ...     index=['s1', 's2', 's3', 's4', 's5', 's6'])

    Now run ``ancombc`` to determine if there are any features that are significantly
    different in abundance between the treatment and the placebo groups.

    >>> result = ancombc(table + 1, metadata, 'grouping')
    >>> result.round(3)
       FeatureID              Covariate  Log2(FC)     SE       W  pvalue  qvalue  \
    0         b1              Intercept     0.719  0.033  21.472   0.000   0.000
    1         b1  grouping[T.treatment]    -1.182  0.385  -3.070   0.002   0.013
    2         b2              Intercept     0.706  0.018  39.548   0.000   0.000
    3         b2  grouping[T.treatment]    -0.537  0.097  -5.562   0.000   0.000
    4         b3              Intercept     0.069  0.087   0.802   0.422   1.000
    5         b3  grouping[T.treatment]     0.068  0.120   0.566   0.571   1.000
    6         b4              Intercept    -0.002  0.015  -0.142   0.887   1.000
    7         b4  grouping[T.treatment]     0.113  0.120   0.946   0.344   1.000
    8         b5              Intercept     0.078  0.065   1.206   0.228   1.000
    9         b5  grouping[T.treatment]     0.004  0.115   0.032   0.974   1.000
    10        b6              Intercept    -0.002  0.015  -0.142   0.887   1.000
    11        b6  grouping[T.treatment]    -0.118  0.072  -1.641   0.101   0.504
    12        b7              Intercept    -0.002  0.015  -0.142   0.887   1.000
    13        b7  grouping[T.treatment]     0.052  0.071   0.741   0.459   1.000
    <BLANKLINE>
        Signif
    0     True
    1     True
    2     True
    3     True
    4    False
    5    False
    6    False
    7    False
    8    False
    9    False
    10   False
    11   False
    12   False
    13   False

    The following example shows how to perform global test using ANCOM-BC to identify
    features that are differentially abundant in more than one time point.

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

    The first DataFrame contains the ANCOM-BC test results that identifies
    differentially abundant feature in treatment as well as in time. The second
    DataFrame contains the ANCOM-BC global test results that indicates whether the
    features are differentially abundant across multiple time points.

    >>> result = ancombc(table+1, metadata, formula="time+treatment", grouping="time")
    >>> result[0][['FeatureID', 'Covariate', 'Signif']]
          FeatureID  Covariate  Signif
        0        Y1  Intercept    True
        1        Y1       time    True
        2        Y1  treatment    True
        3        Y2  Intercept    True
        4        Y2       time    True
        5        Y2  treatment    True
    >>> result[1].round(3)
          FeatureID       W  pvalue  qvalue  Signif
        0        Y1  11.625   0.001   0.002    True
        1        Y2  12.069   0.001   0.002    True

    """
    # Note: A pseudocount should have been added to the table by the user prior to
    # calling this function.
    matrix, samples, features = _ingest_table(table)
    _check_composition(np, matrix, nozero=True)
    n_feats = matrix.shape[1]

    # Log-transform count matrix.
    matrix = np.log(matrix)

    # Validate metadata and cast to numbers where applicable.
    metadata = _check_metadata(metadata, matrix, samples)
    metadata = _type_cast_to_float(metadata)

    # Create a design matrix based on metadata and formula.
    dmat = dmatrix(formula, metadata)

    # Obtain a list of covariates by selecting the relevant columns.
    covars = dmat.design_info.column_names
    n_covars = len(covars)

    # validate parameters
    if not 0 < alpha < 1:
        raise ValueError("`alpha`=%f is not within 0 and 1." % alpha)

    # Estimate initial model parameters.
    var_hat, beta, _, vcov_hat = _estimate_params(matrix, dmat)

    # Estimate and correct for sampling bias via expectation-maximization (EM).
    bias = np.empty((n_covars, 3))
    for i in range(n_covars):
        bias[i] = _estimate_bias_em(beta[i], var_hat[:, i], tol=tol, max_iter=max_iter)
    delta_em = bias[:, 0]

    # Correct coefficients (logFC) according to estimated bias.
    beta_hat = beta.T - delta_em

    # Calculate statistics while performing multiple testing correction.
    se_hat, W, pval, qval = _calc_statistics(beta_hat, var_hat, method=p_adjust)

    # Identify significantly differentially abundance feature-covariate pairs.
    reject = qval <= alpha

    # Output the primary results.
    res = pd.DataFrame.from_dict(
        {
            "FeatureID": [x for x in features for _ in range(n_covars)],
            "Covariate": list(covars) * n_feats,
            "Log2(FC)": beta_hat.ravel(),
            "SE": se_hat.ravel(),
            "W": W.ravel(),
            "pvalue": pval.ravel(),
            "qvalue": qval.ravel(),
        }
    )
    # Pandas' nullable boolean type
    res["Signif"] = pd.Series(reject.ravel(), dtype="boolean")

    # Output global test results
    if grouping is not None and len(metadata[grouping].unique()) >= 3:
        res_global = _global_test(dmat, grouping, beta_hat, vcov_hat, alpha, p_adjust)
        res_global_df = pd.DataFrame.from_dict(
            {
                "FeatureID": features,
                "W": res_global[0],
                "pvalue": res_global[1],
                "qvalue": res_global[2],
                "Signif": res_global[3],
            }
        )
        return res, res_global_df
    return res  # , None


def _estimate_params(data, dmat):
    """Estimate initial model parameters.

    Perform initial estimation of model parameters (coefficients, variances and
    mean residuals) based on the observed data prior to bias correction.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Data table. Must be zero-handled and log-transformed.
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.

    Returns
    -------
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances of regression coefficients.
    beta : ndarray of shape (n_features, n_covariates)
        Estimated coefficients (log-fold changes before correction).
    theta : ndarray of shape (n_samples,)
        Per-sample mean residuals of estimated data.
    beta_covmat : ndarray of shape (n_features, n_covariates, n_covariates)
        Estimated covariance matrices.

    """
    # The original R code performs iterative maximum likelihood estimation to calculate
    # the coefficients, and calls MASS:ginv to calculate the pseudo inverse Gram matrix
    # using the Moore-Penrose method. We noticed that the former can be replaced with
    # ordinary least squares by NumPy's `lstsq`; and the latter matches the method by
    # NumPy's `pinv`. We further noticed that both functions perform singular value
    # decomposition (SVD). Therefore, the following code only performs SVD once, and
    # uses the intermediates to calculate coefficients and inverse Gram matrix.

    # Perform thin SVD and set small singular values to zero. The process and threshold
    # (1e-15) are consistent with the underlying algorithm of `pinv`.
    # Note: the Gram matrix may be singular, if there are colinear covariates in the
    # design matrix. In such cases, a direct `inv` will raise an error, and `pinv` is
    # the robust choice.
    U, S, Vh = np.linalg.svd(dmat, full_matrices=False)
    S_inv = np.where(S > 1e-15 * np.max(S), 1.0 / S, 0.0)

    # Regression coefficients
    V = Vh.T
    dmat_inv = (V * S_inv) @ U.T
    beta = dmat_inv @ data

    # Inverse Gram matrix
    gmat_inv = (V * S_inv**2) @ Vh

    # Per-sample mean residuals (theta)
    diff = data - dmat @ beta
    theta = np.mean(diff, axis=1, keepdims=True)

    # Centered residuals
    eps = diff - theta

    # Calculate the covariance matrix of the coefficients. The estimated variances are
    # the diagonal of the covariance matrix.
    # Note: The original R code uses nested `for` loops over samples and features. The
    # current implementation is fully vectorized.
    # Note: The original R code patches NaN with 0.1 when calculating covariances. This
    # is not needed in the current implementation, as the input data are guaranteed to
    # contain only finite real numbers.
    intm_mat = np.einsum("ji,jp,jq->ipq", eps**2, dmat, dmat)
    beta_covmat = (gmat_inv @ intm_mat) @ gmat_inv
    var_hat = np.diagonal(beta_covmat, axis1=1, axis2=2)

    # Note: Residuals are needed for ANCOM-BC2.
    return var_hat, beta, theta.reshape(-1), beta_covmat


def _estimate_bias_em(beta, var_hat, tol=1e-5, max_iter=100):
    """Estimate sampling bias through an expectation-maximization (EM) algorithm.

    This function models the observed coefficients (log-fold changes) for a given
    covariate as a Gaussian mixture distribution with three components: 0) null,
    1) negative, 2) positive. It aims to estimate a global bias term (delta) that
    affects all features.

    Parameters
    ----------
    beta : ndarray of shape (n_features,)
        Estimated coefficients (log-fold change before correction).
    var_hat : ndarray of shape (n_features,)
        Estimated variances of regression coefficients.
    tol : float, optional
        Absolute tolerance of EM iteration. Default is 1e-5.
    max_iter : int
        Max number of iteration of EM iteration. Default is 100.

    Returns
    -------
    delta_em : float
        EM estimator of bias.
    delta_wls : float
        WLS estimator of bias.
    var_delta : float
        Estimated variances of bias.

    """
    # This function involves careful memory optimization to avoid creating temporary
    # arrays during iteration. Technically, the arrays could have been pre-allocated
    # in the outer function `ancombc` and re-used across covariates, which could
    # further enhance memory efficiency. It is left as-is for modularity.

    # The original R code has `na.rm = TRUE` in many commands. This is not necessary
    # in the current implementation, because the pre-correction coefficients (beta)
    # is guaranteed to not contain NaN values.

    # Mask NaN values (deemed unnecessary; left here for future examination).
    # beta = beta[~np.isnan(beta)]

    # There might be a chance that numeric optimization produces (near) zero weights
    # (pi) or variances (kappa), which can cause numerical stability issues in the
    # EM process. To safe-guard, one may use a small number `eps` as the floor of
    # those parameters. The original R code doesn't have this mechanism. Therefore, it
    # is currently disabled.
    # eps = 1e-12

    # Initial model parameters
    pi0, pi1, pi2 = 0.75, 0.125, 0.125  # weights of components (pi)
    delta, l1, l2, kappa1, kappa2 = _init_bias_params(beta)
    params = np.array([pi0, pi1, pi2, delta, l1, l2, kappa1, kappa2])
    updated = np.empty(8)

    # Pre-allocate memory for intermediates. Each array has three rows, representing
    # the three components (0, 1, 2), and columns representing individual features.
    # Let delta_i = mu_i1_hat - mu_i2_hat
    # The distribution of delta_i is modeled by Gaussian mixture:
    # f(delta_i) = pi0 * phi((delta_i - delta) / nu_i0) +
    #              pi1 * phi((delta_i - (delta + l1)) / nu_i1) +
    #              pi2 * phi((delta_i - (delta + l2)) / nu_i2)
    # where phi is the normal density function,
    # (delta + l1) and (delta + l2) are means for delta_i | C1 and delta_i | C2
    # nu_i0, nu_i1, and nu_i2 are variances of delta_i | C0, delta_i | C1, and
    # delta_i | C2 respectively.
    # We assume nu_i1 = nu_i0 + kappa1 and nu_i2 = nu_i0 + kappa2 for computational
    # simplicity.
    n_feats = beta.shape[0]
    shape = (3, n_feats)
    nu_inv = np.empty(shape)  # inverse of variances
    stdevs = np.empty(shape)  # standard deviations
    ratios = np.empty(shape)  # coefficients / variances

    # Mean coefficients
    means = np.empty(3)

    # Posterior probabilities of feature-component assignments (EM's responsibilities)
    resp = np.empty(shape)

    # Just a 2-row array to store random data
    intm = np.empty((2, n_feats))

    # Initialize intermediates. The 1st row is constant, representing pre-correction
    # estimates, whereas the 2nd and 3rd rows are to be modified during iteration.
    np.reciprocal(var_hat, out=nu_inv[0])
    np.sqrt(var_hat, out=stdevs[0])
    np.divide(beta, var_hat, out=ratios[0])

    # Objective function for numeric optimization of variance estimation
    # Note: `norm.logpdf` doesn't have an `out` parameter. To furhter optimize this,
    # one needs to manually implement the under-the-hood algorithm.
    def func(x, loc, resp):
        log_pdf = norm.logpdf(beta, loc=loc, scale=(var_hat + x) ** 0.5)
        return -np.dot(resp, log_pdf)

    # Optimizer arguments. The Nelder-Mead simplex algorithm is used, which is
    # consistent with the original R code.
    # Note: The Nelder–Mead method doesn't actually enforce bounds during optimization.
    # It merely clips at the bounds. For further consideration.
    args = dict(method="Nelder-Mead", bounds=[(0, None)])

    # Expectation-maximization (E-M) iterations
    loss, epoch = np.inf, 0
    while loss > tol and epoch < max_iter:
        # Update intermediates (2nd and 3rd rows only)
        np.add(var_hat, params[6:8, None], out=intm)  # kappa1, kappa2
        np.reciprocal(intm, out=nu_inv[1:])
        np.sqrt(intm, out=stdevs[1:])
        np.subtract(beta, params[4:6, None], out=ratios[1:])  # means (l)
        ratios[1:] *= nu_inv[1:]

        ### E-step ###
        # Mean coefficients
        delta = means[0] = params[3]  # global bias (delta)
        np.add(params[4:6], delta, means[1:])

        # Posterior probabilities = mean probability density functions weighted by
        # component fractions
        # Note: `norm.pdf` doesn't have an `out` parameter. To further optimize this,
        # one needs to manually implement the under-the-hood algorithm.
        # p_r,i = (pi_r * phi(delta_i - (delta + l_r) / nu_ir)) /
        #         sum_r(pi_r * phi((delta_i - (dleta + l_r)) / nu_ir)),
        # where r = 0, 1, 2; i = 1, ..., n_features
        resp[:] = norm.pdf(beta, means[:, None], stdevs)
        resp *= params[:3, None]  # weights (pi)
        resp /= np.sum(resp, axis=0, keepdims=True)

        ### M-step ###
        # Weights of components (pi)
        # pi_r_new = mean(pi_r * pdf_r / (pi0 * pdf0 + pi1 * pdf1 + pi2 * pdf2)),
        # where r = 0, 1, 2
        np.mean(resp, axis=1, out=updated[:3])

        # Avoid zero weights.
        # updated[:3] = np.maximum(updated[:3], eps)
        # updated[:3] /= updated[:3].sum()

        # Gaussian mixture modeling of global bias (delta)
        # The following code produces the same result as:
        #   updated[3] = np.sum(resp * ratios) / np.sum(resp * nu_inv)
        # But it avoids creating intermediate arrays.
        # delta_new = sum(r_0i * beta / nu0 +
        #                 r_1i * (beta - l1) / (nu0 + kappa1) +
        #                 r_2i * (beta - l2) / (nu0 + kappa2)) /
        #             sum(r0i / nu0 + r1i / (nu0 + kappa1) + r2i / (nu0 + kappa2))
        updated[3] = np.vdot(resp, ratios) / np.vdot(resp, nu_inv)

        # Negative and positive components relative to delta (l)
        np.multiply(resp[1:], nu_inv[1:], out=intm)
        denom = np.sum(intm, axis=1)
        np.subtract(beta, delta, out=resp[0])  # reuse as intermediate
        intm *= resp[0]
        numer = np.sum(intm, axis=1)
        l1, l2 = numer / denom
        # l1_new = min(sum(r1i * (beta - delta) / (nu0 + kappa1)) /
        #              sum(r1i / (nu0 + kappa1)), 0)
        updated[4] = np.minimum(l1, 0)
        # l2_new = min(sum(r2i * (beta - delta) / (nu0 + kappa2)) /
        #              sum(r2i / (nu0 + kappa2)), 0)
        updated[5] = np.maximum(l2, 0)

        # Perform numeric optimization to minimize variances of negative and positive
        # components (kappa).
        # TODO: Consider scenarios where optimization doesn't converge.
        updated[6] = minimize(func, params[6], args=(means[1], resp[1]), **args).x[0]
        updated[7] = minimize(func, params[7], args=(means[2], resp[2]), **args).x[0]

        # Avoid zero variances.
        # updated[6] = max(updated[6], eps)
        # updated[7] = max(updated[7], eps)

        # Loss (epsilon)
        # epsilon = sqrt((pi0_new - pi0)^2 + (pi1_new - pi1)^2 + (pi2_new - pi2)^2 +
        #                (delta_new - delta)^2 + (l1_new - l1)^2 + (l2_new - l2)^2 +
        #                (kappa1_new - kappa1)^2 + (kappa2_new - kappa2)^2)
        loss = np.linalg.norm(updated - params)

        params[:] = updated
        epoch += 1

    return _estimate_bias_var(beta, var_hat, params)


def _init_bias_params(beta):
    """Initialize parameters for iterative bias estimation.

    Parameters
    ----------
    beta : ndarray of shape (n_features,)
        Estimated coefficients before correction.

    Returns
    -------
    delta : float
        Initial global bias term.
    l1, l2 : float
        Initial means of negative and positive components, respectively.
    kappa1, kappa2 : float
        Initial variances of negative and positive components, respectively.

    """
    edges = np.quantile(beta, (0.125, 0.25, 0.75, 0.875))

    # Estimate delta (mean of values between q1 and q3)
    if np.any(mask := (beta >= edges[1]) & (beta <= edges[2])):
        delta = np.mean(beta[mask])
    else:
        delta = np.mean(beta)

    # Estimate l1
    if np.any(mask := beta < edges[0]):
        l1 = np.mean(beta_ := beta[mask])
        if beta_.size > 1:
            kappa1 = np.var(beta_, ddof=1, mean=l1) or 1.0
        else:
            kappa1 = 1.0
    else:
        l1 = np.min(beta)
        kappa1 = 1.0

    # Estimate l2
    if np.any(mask := beta > edges[3]):
        l2 = np.mean(beta_ := beta[mask])
        if beta_.size > 1:
            kappa2 = np.var(beta_, ddof=1, mean=l2) or 1.0
        else:
            kappa2 = 1.0
    else:
        l2 = np.max(beta)
        kappa2 = 1.0

    return delta, l1, l2, kappa1, kappa2


def _estimate_bias_var(beta, var_hat, params):
    """Estimate variance of sampling bias according to EM-optimized parameters.

    This function calculates the weighted least squares (WLS) estimator of bias.

    Parameters
    ----------
    beta : ndarray of shape (n_features,)
        Estimated coefficients (log-fold change before correction).
    var_hat : ndarray of shape (n_features,)
        Estimated variances of regression coefficients.
    params : ndarray of shape (8,)
        Model parameters optimized by EM.

    Returns
    -------
    delta_em : float
        EM estimator of bias.
    delta_wls : float
        WLS estimator of bias.
    var_delta : float
        Estimated variance of bias.

    """
    pi0, pi1, pi2, delta, l1, l2, kappa1, kappa2 = params

    # The EM estimator of bias
    delta_em = delta

    # Assign features to Gaussian components (C)
    q1, q2 = np.quantile(beta, (pi1, 1.0 - pi2))
    C1 = np.flatnonzero(beta < q1)
    C2 = np.flatnonzero(beta >= q2)

    # Numerator of the WLS estimator
    nu_inv = var_hat.copy()
    nu_inv[C1] += kappa1
    nu_inv[C2] += kappa2
    nu_inv[:] = 1.0 / nu_inv
    wls_denom = np.sum(nu_inv)

    # Denominator of the WLS estimator
    beta_ = beta.copy()
    beta_[C1] -= l1
    beta_[C2] -= l2
    nu_inv *= beta_
    wls_numer = np.sum(nu_inv)

    # Estimate the variance of bias
    wls_denom_inv = 1.0 / wls_denom
    delta_wls = wls_numer * wls_denom_inv
    var_delta = np.nan_to_num(wls_denom_inv)

    # TODO: var_delta will be used if conserve=True to account for the variance of
    # delta_hat
    return delta_em, delta_wls, var_delta


def _sample_fractions(data, dmat, beta_hat):
    """Estimate sampling fractions.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Data table. Zero-handled. Log-transformed.
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.
    beta_hat : ndarray of shape (n_features, n_covariates)
        Corrected coefficients.

    Returns
    -------
    theta_hat : ndarray of shape (n_samples,)
        Sampling fractions.

    """
    return np.mean(data - dmat @ beta_hat.T, axis=1)


def _calc_statistics(beta_hat, var_hat, method="holm"):
    """Calculate statistical significance while correcting for multiple testing.

    Parameters
    ----------
    beta_hat : ndarray of shape (n_features, n_covariates)
        Estimated coefficients post correction.
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances.

    Returns
    -------
    se_hat : ndarray of shape (n_features, n_covariates)
        Estimated standard errors.
    W : ndarray of shape (n_features, n_covariates)
        Test statistics.
    p : ndarray of shape (n_features, n_covariates)
        p-values.
    q : ndarray of shape (n_features, n_covariates)
        Adjusted p-values.

    """
    se_hat = var_hat**0.5
    W = beta_hat / se_hat
    pval = 2.0 * norm.sf(abs(W), loc=0, scale=1)
    func = _check_p_adjust(method)
    qval = np.apply_along_axis(func, 0, pval)
    return se_hat, W, pval, qval


def struc_zero(table, metadata, grouping, neg_lb=False):
    """Identify features with structural zeros.

    The function returns a boolean matrix of features with structural zeros, i.e.,
    observerd zeros due to systematical absence.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing count or proportional abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    metadata : pd.DataFrame or 2-D array_like
        The metadata for the model. Rows correspond to samples and columns correspond
        to covariates in the model. Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    grouping : str
        A metadata column name indicating the assignment of samples to groups.
    neg_lb : bool, optional
        Determine whether to use negative lower bound when calculating sample
        proportion. Default is False. Generally, it is recommended to set `neg_lb=True`
        when the sample size per group is relatively large.

    Returns
    -------
    pd.DataFrame
        A table of whether the features are structural zeros in groups (True: structural
        zero, False: not structural zero).

    References
    ----------
    .. [1] Lin, H. and Peddada, S.D., 2020. Analysis of compositions of microbiomes
       with bias correction. Nature communications, 11(1), p.3514.

    Examples
    --------
    >>> from skbio.stats.composition import ancombc
    >>> import pandas as pd

    Generate a DataFrame with 10 samples and 6 features with 0's in specific groups:

    >>> table = pd.DataFrame([[ 7,  1,  0, 11,  3,  1],
    ...                       [ 1,  1,  0, 13, 13,  0],
    ...                       [11,  5,  0,  1,  4,  1],
    ...                       [ 2,  2,  0, 16,  4,  0],
    ...                       [ 0,  1,  0,  0,  6,  0],
    ...                       [14,  8,  7,  9,  0,  5],
    ...                       [ 0,  7,  4,  1,  0, 26],
    ...                       [ 8,  1,  4, 28,  0, 10],
    ...                       [ 2,  2,  2,  4,  0,  5],
    ...                       [ 6,  4, 10,  1,  0,  9]],
    ...                      index=[f's{i}' for i in range(10)],
    ...                      columns=[f'f{i}' for i in range(6)])

    Then create a grouping vector. In this example, there is a treatment group
    and a placebo group.

    >>> metadata = pd.DataFrame(
    ...     {'grouping': ['treatment'] * 5 + ['placebo'] * 5},
    ...     index=[f's{i}' for i in range(10)])

    ``struc_zero`` function will detect features with structural zero. Features that
    are identified as structural zeros in given group are not used in furthur
    analysis such as ``ancombc`` and  ``dirmult_ttest``.
    Setting `neg_lb=True` declares that the true prevalence of a feature in a group is
    not significantly different from zero.

    >>> result = struc_zero(table, metadata, grouping="grouping", neg_lb=True)
    >>> result
                  f0     f1     f2     f3     f4     f5
    placebo    False  False  False  False   True  False
    treatment  False  False   True  False  False   True

    """
    # Validate feature table and metadata
    matrix, samples, features = _ingest_table(table)
    metadata = _check_metadata(metadata, matrix, samples)
    metadata = _type_cast_to_float(metadata)

    unique_groups, group_indices, group_counts = np.unique(
        metadata[grouping], return_inverse=True, return_counts=True
    )

    # Create a boolean matrix to indicate whether the value in table is 0 or not
    tmp = np.nan_to_num(matrix) != 0

    n_groups = len(unique_groups)
    n_features = tmp.shape[1]

    # Initialize group sum matrix
    group_sums = np.zeros((n_groups, n_features))
    np.add.at(group_sums, group_indices, tmp.astype(int))

    # Calculate sample sizes of the groups for each feature
    sample_size = group_counts[:, np.newaxis]

    # Calculate sample proportions of the groups for each feature
    p_hat = group_sums / sample_size

    # Calculate the lower bound of a 95% confidence interval for sample proportion
    if neg_lb:
        p_hat = p_hat - 1.96 * (p_hat * (1 - p_hat) / sample_size) ** 0.5

    zero_idx = p_hat <= 0

    # Output structural zero as a DataFrame
    res = pd.DataFrame(zero_idx, columns=features, index=unique_groups)
    return res


def _global_test(dmat, grouping, beta_hat, vcov_hat, alpha=0.05, p_adjust="holm"):
    """Output results from ANCOM-BC global test to determine feature that are
    differentially abundant between at least 2 groups across 3 or more different
    groups.

    Parameters
    ----------
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.
    grouping : str
        The group variable of interests in metadata.
    beta_hat : ndarray of shape (n_features, n_covariates)
        Corrected coefficients.
    vcov_hat : ndarray of shape (n_features, n_covariates, n_covariates)
        Estimated covariance matrices.
    alpha : float, optional
        Significance level for the statistical tests. Must be in the range of (0, 1).
        Default is 0.05.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are Holm-
        Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-
        Hochberg ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported
        by statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.

    Returns
    -------
    res_W : ndarray of shape (n_features,)
        W of global test.
    res_p : ndarray of shape (n_features,)
        p values of global test.
    qval : ndarray of shape (n_features,)
        adjusted p value of global test.
    reject : ndarray of shape (n_features,)
        if the variable is differentially abundant.

    """
    # Slices of columns in the dmat that the terms in the grouping is mapped to
    s = dmat.design_info.term_name_slices[grouping]

    # Get the index of the terms in the grouping
    group_ind = np.array(range(*s.indices(s.stop)))

    # Subset beta_hat and vcov_hat by grouping indices
    beta_hat_sub = beta_hat[:, group_ind]
    vcov_hat_sub = vcov_hat[:, group_ind][:, :, group_ind]

    # Inverse the subset of vcov_hat
    vcov_hat_sub_inv = np.linalg.pinv(vcov_hat_sub)

    dof = group_ind.size
    A = np.identity(dof)

    # for each feature, calcualte test statistics W by the following formula:
    # W = (A @ beta_hat_sub).T @ inv(A @ vcov_hat_sub @ A.T) @ (A @ beta_hat_sub)
    term = np.einsum("ik,jk->ji", A, beta_hat_sub)
    W_global = np.einsum("ni,nij,ni->n", term, vcov_hat_sub_inv, term)

    # Derive p-values from W statistics
    p_lower = chi2.cdf(W_global, dof)
    p_upper = chi2.sf(W_global, dof)
    pval = 2 * np.minimum(p_lower, p_upper)

    # Correct p-values
    func = _check_p_adjust(p_adjust)
    qval = np.apply_along_axis(func, 0, pval)
    reject = qval <= alpha

    return W_global, pval, qval, reject
