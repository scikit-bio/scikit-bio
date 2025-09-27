# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize
from patsy import dmatrix
from skbio.table._tabular import _ingest_table
from skbio.stats.composition import (
    _check_composition,
    _check_metadata,
    _check_p_adjust,
    _type_cast_to_float,
)


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

    # regression coefficients
    V = Vh.T
    dmat_inv = (V * S_inv) @ U.T
    beta = dmat_inv @ data

    # inverse Gram matrix
    gmat_inv = (V * S_inv**2) @ Vh

    # per-sample mean residuals (theta)
    diff = data - dmat @ beta
    theta = np.mean(diff, axis=1, keepdims=True)

    # centered residuals
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
    return var_hat, beta, theta.reshape(-1)


def _bias_params_init(beta):
    """Initialize parameters for iterative bias estimation.

    Parameters
    ----------
    beta : ndarray of shape (n_features, n_covariates)
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
    edges = np.quantile(beta, [0.125, 0.25, 0.75, 0.875])

    # estimate delta (mean of values between q1 and q3)
    if np.any(mask := (beta >= edges[1]) & (beta <= edges[2])):
        delta = np.mean(beta[mask])
    else:
        delta = np.mean(beta)

    # estimate l1
    if np.any(mask := beta < edges[0]):
        l1 = np.mean(beta_ := beta[mask])
        if beta_.size > 1:
            kappa1 = np.var(beta_, ddof=1, mean=l1) or 1.0
        else:
            kappa1 = 1.0
    else:
        l1 = np.min(beta)
        kappa1 = 1.0

    # estimate l2
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


def _estimate_bias_em(beta, var_hat, tol=1e-5, max_iter=100):
    """Estimate sampling bias through an expectation-maximization (EM) algorithm.

    This function models the observed coefficients (log-fold changes) for a given
    covariate as a Gaussian mixture distribution with three components: 0) null,
    1) negative, 2) positive. It aims to estimate a global bias term (delta) that
    affects all features.

    Parameters
    ----------
    beta : ndarray of shape (n_features, n_covariates)
        Estimated coefficients (log-fold change before correction).
    var_hat : ndarray of shape (n_features, n_covariates)
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
    # The original R code has `na.rm = TRUE` in many commands. This is not necessary
    # in the current implementation, because the pre-correction coefficients (beta)
    # is guaranteed to not contain NaN values.

    # Mask NaN values (deemed unnecessary; left here for future examination).
    # beta = beta[~np.isnan(beta)]

    # Initial model parameters
    pi0, pi1, pi2 = 0.75, 0.125, 0.125  # proportions of components (pi)
    delta, l1, l2, kappa1, kappa2 = _bias_params_init(beta)
    params = np.array([pi0, pi1, pi2, delta, l1, l2, kappa1, kappa2])
    updated = np.empty(8)

    # Pre-allocate memory for intermediates. Each array has three rows, representing
    # the three components (0, 1, 2), and columns representing individual features.
    n_feats = beta.shape[0]
    shape = (3, n_feats)
    nu_inv = np.empty(shape)  # inverse of variances
    stdevs = np.empty(shape)  # standard deviations
    ratios = np.empty(shape)  # coefficients / variances

    # mean coefficients
    means = np.empty(3)

    # posterior probabilities of feature-component assignments (EM's responsibilities)
    resp = np.empty(shape)

    # just a 2-row array to store random data
    intm = np.empty((2, n_feats))

    # Initialize intermediates. The 1st row is constant, representing pre-correction
    # estimates, whereas the 2nd and 3rd rows are to be modified during iteration.
    nu_inv[0] = 1.0 / var_hat
    stdevs[0] = var_hat**0.5
    ratios[0] = beta / var_hat

    # Nelder-Mead simplex algorithm for variance estimation
    def func(x, loc, resp):
        log_pdf = norm.logpdf(beta, loc=loc, scale=(var_hat + x) ** 0.5)
        return -np.sum(resp * log_pdf)

    # optimizer arguments
    args = dict(method="Nelder-Mead", bounds=[(0, None)])

    # expectation-maximization (E-M) iterations
    loss, epoch = np.inf, 0
    while loss > tol and epoch < max_iter:
        # update intermediates (2nd and 3rd rows only)
        intm[:] = var_hat + params[6:8][:, None]  # variances (kappa)
        nu_inv[1:] = 1.0 / intm
        stdevs[1:] = intm**0.5
        ratios[1:] = (beta - params[4:6][:, None]) * nu_inv[1:]  # means (l)

        ### E-step ###
        # mean coefficients
        delta = means[0] = params[3]  # global error (delta)
        means[1:] = delta + params[4:6]

        # posterior probabilities = mean probability density functions weighted by
        # component fractions
        resp[:] = norm.pdf(beta, means[:, None], stdevs)
        resp *= params[:3][:, None]
        resp /= np.sum(resp, axis=0, keepdims=True)

        ### M-step ###
        # proportions of components (pi)
        np.mean(resp, axis=1, out=updated[:3])

        # Gaussian mixture modeling of global error (delta)
        # The following code produces the same result as:
        #   updated[3] = np.sum(resp * ratios) / np.sum(resp * nu_inv)
        # But it avoids creating intermediate arrays.
        updated[3] = np.tensordot(resp, ratios) / np.tensordot(resp, nu_inv)

        # negative and positive components relative to delta (l)
        intm[:] = resp[1:] * nu_inv[1:]
        denom = np.sum(intm, axis=1)
        intm *= beta - delta
        numer = np.sum(intm, axis=1)
        l1_new, l2_new = numer / denom
        updated[4] = np.minimum(l1_new, 0)
        updated[5] = np.maximum(l2_new, 0)

        # Perform numeric optimization to minimize variances of negative and positive
        # components (kappa).
        # TODO: Consider scenarios where optimization won't converge.
        updated[6] = minimize(func, params[6], args=(means[1], resp[1]), **args).x[0]
        updated[7] = minimize(func, params[7], args=(means[2], resp[2]), **args).x[0]

        # loss (epsilon)
        loss = np.linalg.norm(updated - params)

        params[:] = updated
        epoch += 1

    return _post_estimate_bias(beta, var_hat, params)


def _post_estimate_bias(beta, var_hat, params):
    """Estimate sample bias according to EM-optimized parameters."""
    pi0, pi1, pi2, delta, l1, l2, kappa1, kappa2 = params
    nu = var_hat.copy()

    # The EM estimator of bias
    delta_em = delta

    # The weighted least squares (WLS) estimator of bias
    q1, q2 = np.quantile(beta, [pi1, 1.0 - pi2])
    C0 = np.where((beta >= q1) & (beta < q2))[0]
    C1 = np.where(beta < q1)[0]
    C2 = np.where(beta >= q2)[0]

    # Numerator of the WLS estimator
    nu[C1] += kappa1
    nu[C2] += kappa2
    nu_inv = 1.0 / nu
    wls_denom = np.sum(nu_inv)

    # Denominator of the WLS estimator
    nu_inv[C0] *= beta[C0]
    nu_inv[C1] *= (beta - l1)[C1]
    nu_inv[C2] *= (beta - l2)[C2]
    wls_nume = np.sum(nu_inv)

    # Estimate the variance of bias
    wls_denom_inv = 1.0 / wls_denom
    delta_wls = wls_nume * wls_denom_inv
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


def ancombc(
    table,
    metadata,
    formula,
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
    pd.DataFrame
        A table of features and covariates, their log-fold changes and other relevant
        statistics.

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

    See Also
    --------
    ancom

    Notes
    -----
    The input data table for ANCOM-BC must contain only positive numbers. One needs to
    remove zero values by, e.g., adding a pseudocount of 1.0.

    References
    ----------
    .. [1] Lin, H. and Peddada, S.D., 2020. Analysis of compositions of microbiomes
       with bias correction. Nature communications, 11(1), p.3514.

    Examples
    --------
    >>> from skbio.stats._ancombc import ancombc
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

    Now run ``ancom`` to determine if there are any features that are
    significantly different in abundance between the treatment and the placebo
    groups. The first DataFrame that is returned contains the ANCOM test
    results, and the second contains the percentile abundance data for each
    feature in each group.

    >>> result = ancombc(table + 1, metadata, 'grouping')

    """
    # Note: A pseudocount should have been added to the table by the user prior to
    # calling this function.
    matrix, samples, features = _ingest_table(table)
    _check_composition(np, matrix, nozero=True)  # TODO
    n_feats = matrix.shape[1]

    # log-transform count matrix
    matrix = np.log(matrix)

    # validate metadata
    metadata = _check_metadata(metadata, matrix, samples)

    # cast metadata to numbers where applicable
    metadata = _type_cast_to_float(metadata)

    # Create a design matrix based on metadata and formula.
    dmat = dmatrix(formula, metadata)

    # Obtain the list of covariates by selecting the relevant columns
    covars = dmat.design_info.column_names
    n_covars = len(covars)

    # validate parameters
    if not 0 < alpha < 1:
        raise ValueError("`alpha`=%f is not within 0 and 1." % alpha)

    # Estimate initial model parameters.
    var_hat, beta, _ = _estimate_params(matrix, dmat)

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
    # pandas' nullable boolean type
    res["Signif"] = pd.Series(reject.ravel(), dtype="boolean")
    return res
