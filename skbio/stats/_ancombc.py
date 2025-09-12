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
from statsmodels.stats.multitest import multipletests


def _est_params(data, dmat):
    """Estimate parameters.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Data table. Zero-handled. Log-transformed.
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.

    Returns
    -------
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances.
    beta_hat : ndarray of shape (n_features, n_covariates)
        Estimated coefficients.
    theta : ndarray of shape (n_samples,)
        Residuals of estimated data.

    """
    n_samples, n_features = data.shape
    n_covariates = dmat.shape[1]

    beta_hat = np.empty((n_covariates, n_features))
    y_crt_hat = np.empty((n_samples, n_features))
    var_hat = np.empty((n_features, n_covariates))

    for i in range(n_features):
        try:
            beta = np.linalg.lstsq(dmat, data[:, i], rcond=None)[0]
        except np.linalg.LinAlgError:  # try case to generate the error
            continue

        beta_hat[:, i] = beta[:n_covariates]
        y_crt_hat[:, i] = dmat.dot(beta)

    XtX = dmat.T.dot(dmat)
    try:
        XtX_inv = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        XtX_inv = np.linalg.pinv(XtX)

    # calculate residual
    diff = data - y_crt_hat
    theta = np.mean(diff, axis=1, keepdims=True)
    eps = diff - theta

    for i in range(n_features):  # need to optimize
        sigma2_xxT = np.zeros((n_covariates, n_covariates))
        for j in range(n_samples):
            sigma2_xxT_j = np.outer(eps[j, i] ** 2 * dmat[j], dmat[j])
            # sigma2_xxT_j[is.na(sigma2_xxT_j)] = 0.1
            sigma2_xxT += sigma2_xxT_j

        # estimated variance is the diagonal of the covariance matrix
        var_hat[i] = np.diagonal(XtX_inv @ sigma2_xxT @ XtX_inv)

    return var_hat, beta_hat, theta.reshape(-1)


def _bias_em(beta, var_hat, tol=1e-5, max_iter=100):
    """Estimate bias.

    Parameters
    ----------
    beta : ndarray of shape (n_features,)
        Data table. Zero-handled. Log-transformed.
    var_hat : ndarray of shape (n_features,)
        Design matrix.
    tol : float, optional
        Tolerance of EM iteration.
    max_iter : int
        Max number of iteration of EM iteration.

    Returns
    -------
    delta_em : float
        EM estimator of bias.
    delta_wls : float
        WLS estimator of bias.
    var_delta : float
        Estimated variances of bias.

    """
    beta = beta.copy()
    nu = var_hat.copy()

    # Initialization
    pi0_0 = 0.75
    pi1_0 = 0.125
    pi2_0 = 0.125

    # Extract valid numbers (not NaN) from beta.
    beta_ = beta[~np.isnan(beta)]

    # keep round temporarily and check results of em after more iterations
    quantiles = np.quantile(beta_, [0.125, 0.25, 0.75, 0.875]).round(6)

    if np.any(mask := (beta_ >= quantiles[1]) & (beta_ <= quantiles[2])):
        delta_0 = np.mean(beta_[mask])
    else:
        delta_0 = np.mean(beta_)

    if np.any(mask := beta_ < quantiles[0]):
        l1_0 = np.mean(beta_[mask])
        kappa1_0 = np.var(beta_[mask], ddof=1)
        if np.isnan(kappa1_0) or kappa1_0 == 0.0:
            kappa1_0 = 1.0
    else:
        l1_0 = np.min(beta_)
        kappa1_0 = 1.0

    if np.any(mask := beta_ > quantiles[3]):
        l2_0 = np.mean(beta_[mask])
        kappa2_0 = np.var(beta_[mask], ddof=1)
        if np.isnan(kappa2_0) or kappa2_0 == 0.0:
            kappa2_0 = 1.0
    else:
        l2_0 = np.max(beta_)
        kappa2_0 = 1.0

    # E-M algorithm
    pi0_vec = [pi0_0]
    pi1_vec = [pi1_0]
    pi2_vec = [pi2_0]
    delta_vec = [delta_0]
    l1_vec = [l1_0]
    l2_vec = [l2_0]
    kappa1_vec = [kappa1_0]
    kappa2_vec = [kappa2_0]

    # E-M iteration
    epoch = 0
    loss = 100
    beta_nu = beta / nu
    while loss > tol and epoch < max_iter:
        # Current value of paras
        pi0 = pi0_vec[-1]
        pi1 = pi1_vec[-1]
        pi2 = pi2_vec[-1]
        delta = delta_vec[-1]
        l1 = l1_vec[-1]
        l2 = l2_vec[-1]
        kappa1 = kappa1_vec[-1]
        kappa2 = kappa2_vec[-1]

        # E-step
        pdf0 = norm.pdf(beta, delta, nu**0.5)
        pdf1 = norm.pdf(beta, delta + l1, (nu + kappa1) ** 0.5)
        pdf2 = norm.pdf(beta, delta + l2, (nu + kappa2) ** 0.5)

        prod0, prod1, prod2 = pi0 * pdf0, pi1 * pdf1, pi2 * pdf2
        deno = prod0 + prod1 + prod2
        r0i, r1i, r2i = prod0 / deno, prod1 / deno, prod2 / deno

        # M-step
        pi0_new = np.nanmean(r0i)
        pi1_new = np.nanmean(r1i)
        pi2_new = np.nanmean(r2i)

        nu0_kappa1 = nu + kappa1
        nu0_kappa2 = nu + kappa2
        beta_delta = beta - delta
        delta_new = np.nansum(
            r0i * beta_nu
            + r1i * (beta - l1) / nu0_kappa1
            + r2i * (beta - l2) / nu0_kappa2
        ) / np.nansum(r0i / nu + r1i / nu0_kappa1 + r2i / nu0_kappa2)
        l1_new = np.nanmin(
            [np.nansum(r1i * beta_delta / nu0_kappa1) / np.nansum(r1i / nu0_kappa1), 0]
        )
        l2_new = np.nanmax(
            [np.nansum(r2i * beta_delta / nu0_kappa2) / np.nansum(r2i / nu0_kappa2), 0]
        )

        # Nelder-Mead simplex algorithm for kappa1 and kappa2
        def obj_kappa(x, ll, ri):
            log_pdf = norm.logpdf(beta, loc=delta + ll, scale=(nu + x) ** 0.5)
            return -np.sum(ri * log_pdf)

        kappa1_new = minimize(
            fun=obj_kappa,
            x0=kappa1,
            args=(l1, r1i),
            method="Nelder-Mead",
            bounds=[(0, None)],
        ).x[0]
        kappa2_new = minimize(
            fun=obj_kappa,
            x0=kappa2,
            args=(l2, r2i),
            method="Nelder-Mead",
            bounds=[(0, None)],
        ).x[0]

        # Merge to the paras vectors/matrices
        pi0_vec.append(pi0_new)
        pi1_vec.append(pi1_new)
        pi2_vec.append(pi2_new)
        delta_vec.append(delta_new)
        l1_vec.append(l1_new)
        l2_vec.append(l2_new)
        kappa1_vec.append(kappa1_new)
        kappa2_vec.append(kappa2_new)

        # Calculate the new epsilon (loss)
        loss = (
            (pi0_new - pi0) ** 2
            + (pi1_new - pi1) ** 2
            + (pi2_new - pi2) ** 2
            + (delta_new - delta) ** 2
            + (l1_new - l1) ** 2
            + (l2_new - l2) ** 2
            + (kappa1_new - kappa1) ** 2
            + (kappa2_new - kappa2) ** 2
        ) ** 0.5
        epoch += 1

    # The EM estimator of bias
    delta_em = delta_new

    # The weighted least squares (WLS) estimator of bias
    pi1 = pi1_new
    pi2 = pi2_new
    l1 = l1_new
    l2 = l2_new
    kappa1 = kappa1_new
    kappa2 = kappa2_new

    q1, q2 = np.quantile(beta, [pi1, 1 - pi2]).round(6)
    C0 = np.where((beta >= q1) & (beta < q2))[0]
    C1 = np.where(beta < q1)[0]  # add round to make the results exactly the same with r
    C2 = np.where(beta >= q2)[0]

    # Numerator of the WLS estimator
    nu[C1] += kappa1
    nu[C2] += kappa2
    nu_inv = 1 / nu
    wls_deno = np.sum(nu_inv)

    # Denominator of the WLS estimator
    nu_inv[C0] *= beta[C0]
    nu_inv[C1] *= (beta - l1)[C1]
    nu_inv[C2] *= (beta - l2)[C2]
    wls_nume = np.sum(nu_inv)

    # Estimate the variance of bias
    with np.errstate(divide="ignore", invalid="ignore"):
        wls_deno_inv = 1 / wls_deno

    delta_wls = wls_nume * wls_deno_inv
    var_delta = np.nan_to_num(wls_deno_inv)

    return delta_em, delta_wls, var_delta


def _correct_coefficients(beta, delta_em):
    beta_hat = beta.copy().T
    beta_hat = beta_hat - delta_em
    return beta_hat


def _sample_fractions(data, dmat, beta, delta_em):
    n_samples, n_features = data.shape
    beta_hat = beta.copy().T
    beta_hat = beta_hat - delta_em
    theta_hat = np.empty((n_samples, n_features))

    for i in range(n_features):
        theta_hat[:, i] = data[:, i] - dmat @ beta_hat[i]
    theta_hat = np.mean(theta_hat, axis=1)

    return theta_hat


def _multi_test(dmat, beta_hat, var_hat, alpha=0.5):
    n_covariates = dmat.shape[1]
    se_hat = var_hat**0.5

    W = beta_hat / se_hat
    p = 2.0 * norm.sf(abs(W), loc=0, scale=1)
    q = np.empty(p.shape)
    diff_abn = np.empty(p.shape)
    for i in range(n_covariates):
        res = multipletests(p[:, i], method="holm", alpha=alpha)
        q[:, i] = res[1]
        diff_abn[:, i] = res[0]

    return se_hat, W, p, q, diff_abn
