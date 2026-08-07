# ----------------------------------------------------------------------------
# ANCOM-BC2 multi-group test helper functions.
#
# Extracted from _ancombc.py for modularity. The implementation includes:
# - Global F-test
# - Pairwise directional test with mdFDR
# - Dunnett's test with mdFDR
# - Trend test with constrained estimation
# - Utility functions: _var_diff, _combn_fun, _combn_fun2,
#   _safe_inverse_spd, _constrain_est, _l_infty
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from itertools import combinations
from scipy.stats import norm, chi2, t
from scipy.optimize import minimize

from ._utils import _check_p_adjust


def _var_diff(vcov_sub):
    """Compute variance of difference for a pair from a covariance sub-matrix.

    For a 2x2 covariance matrix [[var_a, cov_ab], [cov_ab, var_b]],
    var(a - b) = var_a + var_b - 2*cov_ab = trace - 2*off_diag_sum.

    For a general kxk matrix, var(sum of signed differences) =
    sum of all elements with appropriate signs. For the simple pairwise
    case with contrast [1, -1], this is trace - 2*sum(off-diagonal).
    """
    vcov_sub = np.asarray(vcov_sub)
    # For the pairwise case: var(a-b) = var(a) + var(b) - 2*cov(a,b)
    # = sum(diag) - 2*sum(off-diagonal upper/lower)
    # = sum(all elements) - 2*sum(off-diagonal)
    # More precisely: c' V c where c = [1, -1] for 2x2
    if vcov_sub.shape == (2, 2):
        return vcov_sub[0, 0] + vcov_sub[1, 1] - 2 * vcov_sub[0, 1]
    else:
        # General case: for contrast vector c, var = c' V c
        # Default contrast is [1, -1, 0, ...] etc.
        off_diag = vcov_sub.sum() - np.trace(vcov_sub)
        return np.trace(vcov_sub) - off_diag


def _combn_fun(vec, func, sep="_"):
    """Apply func to all pairs of elements in vec.

    Returns a dict with individual elements and pairwise differences.
    """
    result = {}
    # Individual elements
    for i, name in enumerate(vec.index if hasattr(vec, "index") else range(len(vec))):
        result[str(name)] = vec[i]

    # Pairwise differences
    names = list(vec.index if hasattr(vec, "index") else range(len(vec)))
    for i, j in combinations(range(len(vec)), 2):
        pair_name = f"{names[j]}{sep}{names[i]}"
        result[pair_name] = vec[j] - vec[i]

    return result


def _combn_fun2(mat, sep="_"):
    """Compute variance of differences for all pairs from a covariance matrix.

    Returns a dict with diagonal variances and pairwise variance of differences.
    """
    result = {}
    names = list(mat.index if hasattr(mat, "index") else range(mat.shape[0]))

    # Diagonal variances
    for i, name in enumerate(names):
        result[str(name)] = mat.iloc[i, i] if hasattr(mat, "iloc") else mat[i, i]

    # Pairwise variance of differences
    for i, j in combinations(range(len(names)), 2):
        pair_name = f"{names[j]}{sep}{names[i]}"
        sub = (
            mat.iloc[[j, i], [j, i]] if hasattr(mat, "iloc") else mat[[j, i]][:, [j, i]]
        )
        result[pair_name] = _var_diff(sub)

    return result


def _ancombc_global_F(
    x, group, beta_hat, vcov_hat, dof=None, p_adj_method="holm", alpha=0.05
):
    """ANCOM-BC2 global test using F-test or chi-square test.

    Parameters
    ----------
    x : pd.DataFrame or ndarray
        Design matrix.
    group : str
        Group variable name.
    beta_hat : ndarray of shape (n_taxa, n_covariates)
        Bias-corrected coefficients.
    vcov_hat : list of ndarray of shape (n_covariates, n_covariates)
        Per-taxon covariance matrices.
    dof : ndarray of shape (n_taxa, n_covariates) or None
        Degrees of freedom. If None, use chi-square test.
    p_adj_method : str
        P-value adjustment method.
    alpha : float
        Significance level.

    Returns
    -------
    pd.DataFrame with columns: W, p_val, q_val, diff_abn
    """
    n_tax = beta_hat.shape[0]
    covariates = (
        x.columns
        if isinstance(x, pd.DataFrame)
        else [f"V{i}" for i in range(x.shape[1])]
    )

    # Identify group-related covariates (excluding interactions)
    group_ind = np.array([group in c and ":" not in c for c in covariates])
    beta_hat_sub = beta_hat[:, group_ind]

    W_global = np.full(n_tax, np.nan)
    p_global = np.ones(n_tax)

    for i in range(n_tax):
        beta_i = beta_hat_sub[i]
        vcov_i = vcov_hat[i][np.ix_(group_ind, group_ind)]
        k = int(np.sum(group_ind))

        try:
            A = np.eye(k)
            # W = beta' (A vcov A')^{-1} beta
            AvA = A @ vcov_i @ A.T
            W = float(beta_i @ np.linalg.pinv(AvA) @ beta_i)

            if dof is not None:
                # F-test
                dof_i = np.unique(dof[i, group_ind])
                if len(dof_i) == 1:
                    dof_i = dof_i[0]
                else:
                    dof_i = np.min(dof_i)
                from scipy.stats import f as f_dist

                p = 2 * min(
                    f_dist.cdf(W, dfn=k, dfd=dof_i), f_dist.sf(W, dfn=k, dfd=dof_i)
                )
            else:
                # Chi-square test
                p = 2 * min(chi2.cdf(W, df=k), chi2.sf(W, df=k))

            W_global[i] = W
            p_global[i] = p
        except Exception:
            pass  # Keep NaN/1.0 defaults

    # Multiple testing correction
    func = _check_p_adjust(p_adj_method)
    q_global = func(p_global)
    q_global = np.where(np.isnan(q_global), 1.0, q_global)
    diff_abn = q_global <= alpha

    result = pd.DataFrame(
        {
            "W": W_global,
            "p_val": p_global,
            "q_val": q_global,
            "diff_abn": diff_abn,
        }
    )
    return result


def _ancombc_pair(
    x, group, beta_hat, var_hat, vcov_hat, dof, fwer_ctrl_method="holm", alpha=0.05
):
    """ANCOM-BC2 pairwise directional test.

    For each pair of group levels, compute the difference in coefficients
    and its variance, then apply mdFDR correction.
    """
    covariates = (
        x.columns
        if isinstance(x, pd.DataFrame)
        else [f"V{i}" for i in range(x.shape[1])]
    )

    # Subset group-related covariates
    group_ind = np.array([group in c and ":" not in c for c in covariates])
    beta_hat_sub = beta_hat[:, group_ind]
    group_covars = [c for c, g in zip(covariates, group_ind) if g]

    # Compute pairwise differences and their variances
    n_tax = beta_hat.shape[0]
    n_group = int(np.sum(group_ind))

    # Generate all pairwise comparisons
    pair_names = []
    pair_indices = []
    for i, j in combinations(range(n_group), 2):
        pair_names.append(f"{group_covars[j]}_{group_covars[i]}")
        pair_indices.append((j, i))

    # Also include individual group effects
    all_names = list(group_covars) + pair_names
    n_comp = len(all_names)

    beta_pair = np.zeros((n_tax, n_comp))
    var_pair = np.zeros((n_tax, n_comp))

    # Individual effects
    for k in range(n_group):
        beta_pair[:, k] = beta_hat_sub[:, k]
        var_pair[:, k] = var_hat[:, np.where(group_ind)[0][k]]

    # Pairwise differences
    for k, (j, i) in enumerate(pair_indices):
        beta_pair[:, n_group + k] = beta_hat_sub[:, j] - beta_hat_sub[:, i]
        for tax in range(n_tax):
            vcov_sub = vcov_hat[tax][
                np.ix_(np.where(group_ind)[0][[j, i]], np.where(group_ind)[0][[j, i]])
            ]
            var_pair[tax, n_group + k] = _var_diff(vcov_sub)

    se_pair = np.sqrt(np.maximum(var_pair, 0))
    W_pair = beta_pair / se_pair

    # Get dof for group covariates
    if dof is not None:
        dof_group = dof[:, group_ind]
    else:
        dof_group = None

    # Apply mdFDR correction
    p_q = _mdfdr_pairwise(
        W=W_pair,
        dof=dof_group,
        fwer_ctrl_method=fwer_ctrl_method,
        x=x,
        group=group,
        beta_hat=beta_hat,
        vcov_hat=vcov_hat,
        alpha=alpha,
        dof_global=dof,
    )

    p_hat = p_q["p_val"]
    q_hat = p_q["q_val"]
    diff_pair = q_hat <= alpha

    return {
        "beta": beta_pair,
        "se": se_pair,
        "W": W_pair,
        "p_val": p_hat,
        "q_val": q_hat,
        "diff_abn": diff_pair,
        "comp_names": all_names,
    }


def _mdfdr_pairwise(
    W, dof, fwer_ctrl_method, x, group, beta_hat, vcov_hat, alpha, dof_global=None
):
    """Mixed directional FDR correction for pairwise tests.

    1. Run global test to screen significant taxa (count R).
    2. Only consider R significant taxa for pairwise p-values.
    3. Adjust at level R * alpha / d.
    """
    n_tax, n_comp = W.shape

    # Step 1: Global test screening
    res_screen = _ancombc_global_F(
        x=x,
        group=group,
        beta_hat=beta_hat,
        vcov_hat=vcov_hat,
        dof=dof_global,
        p_adj_method="BH",
        alpha=alpha,
    )

    R = max(int(res_screen["diff_abn"].sum()), 1)
    screen_ind = res_screen["diff_abn"].values

    # Step 2: Compute p-values from t-distribution
    if dof is not None:
        # Per-taxon, per-comparison dof
        # For pairwise comparisons, approximate dof
        dof_pair = np.full_like(W, 999.0)
        if dof.ndim == 2:
            # Use minimum dof across group covariates as approximation
            for i in range(n_tax):
                dof_pair[i, :] = np.nanmin(dof[i, :])
        p_val = 2 * t.sf(np.abs(W), df=dof_pair)
    else:
        p_val = 2 * norm.sf(np.abs(W))

    # Zero out p-values for taxa not significant in global test
    p_val = p_val * screen_ind[:, np.newaxis]
    p_val[p_val == 0] = 1.0
    p_val = np.where(np.isnan(p_val), 1.0, p_val)

    # Step 3: Adjust p-values
    func = _check_p_adjust(fwer_ctrl_method)
    q_val = np.zeros_like(p_val)
    for j in range(n_comp):
        q_val[:, j] = func(p_val[:, j])

    return {"p_val": p_val, "q_val": q_val}


def _dunn_global(x, group, W, B, dof, p_adj_method="holm", alpha=0.05):
    """Dunnett's global test for mdFDR.

    Bootstrap-based: generate null W from t-distribution, take max |W|.
    """
    n_tax = W.shape[0]
    covariates = (
        x.columns
        if isinstance(x, pd.DataFrame)
        else [f"V{i}" for i in range(x.shape[1])]
    )
    group_ind = np.array([group in c and ":" not in c for c in covariates])
    n_group = int(np.sum(group_ind))

    # Observed global statistic: max |W| per taxon
    W_global = np.max(np.abs(W), axis=1)

    # Bootstrap null distribution
    W_global_null = np.zeros((n_tax, B))
    rng = np.random.default_rng(42)

    for b in range(B):
        # Generate null W from t-distribution
        if dof is not None:
            # Use per-taxon, per-group dof
            W_null = np.zeros_like(W)
            for i in range(n_tax):
                for j in range(W.shape[1]):
                    df_val = dof[i, j] if not np.isnan(dof[i, j]) else 999
                    W_null[i, j] = rng.standard_t(df_val)
        else:
            W_null = rng.standard_normal(size=W.shape)

        W_global_null[:, b] = np.max(np.abs(W_null), axis=1)

    # P-values from bootstrap
    p_global = np.mean(W_global_null > W_global[:, np.newaxis], axis=1)

    # Multiple testing correction
    func = _check_p_adjust(p_adj_method)
    q_global = func(p_global)
    q_global = np.where(np.isnan(q_global), 1.0, q_global)
    diff_abn = q_global <= alpha

    return pd.DataFrame(
        {
            "W": W_global,
            "p_val": p_global,
            "q_val": q_global,
            "diff_abn": diff_abn,
        }
    )


def _ancombc_dunn(
    x, group, beta_hat, var_hat, dof, B=100, fwer_ctrl_method="holm", alpha=0.05
):
    """ANCOM-BC2 Dunnett's type of test.

    Compare each group to the reference group with mdFDR correction.
    """
    covariates = (
        x.columns
        if isinstance(x, pd.DataFrame)
        else [f"V{i}" for i in range(x.shape[1])]
    )
    group_ind = np.array([group in c and ":" not in c for c in covariates])

    beta_hat_dunn = beta_hat[:, group_ind]
    var_hat_dunn = var_hat[:, group_ind]
    se_hat_dunn = np.sqrt(np.maximum(var_hat_dunn, 0))
    W_dunn = beta_hat_dunn / se_hat_dunn

    if dof is not None:
        dof_dunn = dof[:, group_ind]
    else:
        dof_dunn = None

    # mdFDR correction
    p_q = _mdfdr_dunnett(
        W=W_dunn,
        dof=dof_dunn,
        fwer_ctrl_method=fwer_ctrl_method,
        x=x,
        group=group,
        B=B,
        alpha=alpha,
    )

    return {
        "beta": beta_hat_dunn,
        "se": se_hat_dunn,
        "W": W_dunn,
        "p_val": p_q["p_val"],
        "q_val": p_q["q_val"],
        "diff_abn": p_q["q_val"] <= alpha,
        "comp_names": [c for c, g in zip(covariates, group_ind) if g],
    }


def _mdfdr_dunnett(W, dof, fwer_ctrl_method, x, group, B, alpha):
    """mdFDR correction for Dunnett's test."""
    n_tax = W.shape[0]

    # Step 1: Global test screening via bootstrap
    res_screen = _dunn_global(
        x=x, group=group, W=W, B=B, dof=dof, p_adj_method="BH", alpha=alpha
    )
    R = max(int(res_screen["diff_abn"].sum()), 1)
    screen_ind = res_screen["diff_abn"].values

    # Step 2: Compute p-values
    if dof is not None:
        p_val = 2 * t.sf(np.abs(W), df=dof)
    else:
        p_val = 2 * norm.sf(np.abs(W))

    # Zero out for non-significant taxa
    p_val = p_val * screen_ind[:, np.newaxis]
    p_val[p_val == 0] = 1.0
    p_val = np.where(np.isnan(p_val), 1.0, p_val)

    # Step 3: Adjust
    func = _check_p_adjust(fwer_ctrl_method)
    q_val = np.zeros_like(p_val)
    for j in range(p_val.shape[1]):
        q_val[:, j] = func(p_val[:, j])

    return {"p_val": p_val, "q_val": q_val}


def _safe_inverse_spd(A, ridge=1e-8):
    """Safe inverse of a symmetric positive-definite matrix."""
    A = (A + A.T) / 2
    try:
        L = np.linalg.cholesky(A)
        return np.linalg.solve(A, np.eye(A.shape[0]))
    except np.linalg.LinAlgError:
        return np.linalg.inv(A + ridge * np.eye(A.shape[0]))


def _constrain_est(beta_hat, vcov_hat, contrast):
    """Constrained estimation via quadratic programming.

    Minimize (beta - beta_opt)' P (beta - beta_opt) subject to
    contrast @ beta_opt >= 0 (each row of contrast defines one inequality).

    Uses scipy.optimize.minimize with SLSQP method, replacing R's
    quadprog::solve.QP.
    """
    P = _safe_inverse_spd(vcov_hat)
    n = len(beta_hat)
    contrast = np.asarray(contrast)

    # Objective: 0.5 * x' D x - d' x  (quadratic form)
    Dmat = 2 * P
    Dmat = (Dmat + Dmat.T) / 2
    # Add small ridge for numerical stability (matches R code: diag(Dmat) += 1e-10)
    Dmat += 1e-10 * np.eye(n)
    dvec = 2 * P @ beta_hat

    def objective(x):
        return 0.5 * x @ Dmat @ x - dvec @ x

    def grad(x):
        return Dmat @ x - dvec

    # Each row of contrast defines one inequality: contrast[i] @ x >= 0
    # SLSQP expects a single constraint function returning a vector
    n_constraints = contrast.shape[0]
    constraints = {
        "type": "ineq",
        "fun": lambda x: contrast @ x,  # Returns vector; each element must be >= 0
        "jac": lambda x: contrast,
    }

    try:
        result = minimize(
            objective,
            x0=beta_hat,
            jac=grad,
            method="SLSQP",
            constraints=constraints,
            options={"ftol": 1e-12, "maxiter": 1000},
        )
        return result.x
    except Exception:
        return np.zeros(n)


def _l_infty(beta_opt, node):
    """Compute l_infinity norm for a pattern at specified node.

    Matches R's .l_infty function:
    l = max(abs(beta_opt[node]), abs(beta_opt[node] - beta_opt[length(beta_opt)]))

    Note: node is 0-based in Python (1-based in R). The last element
    beta_opt[-1] corresponds to beta_opt[length(beta_opt)] in R.
    """
    beta_opt = np.asarray(beta_opt)
    return max(
        abs(beta_opt[node]),
        abs(beta_opt[node] - beta_opt[-1]),
    )


def _ancombc_trend(
    x,
    group,
    beta_hat,
    var_hat,
    vcov_hat,
    p_adj_method="holm",
    alpha=0.05,
    trend_contrast=None,
    trend_node=None,
    trend_B=100,
):
    """ANCOM-BC2 trend test (pattern analysis).

    Uses constrained optimization to test ordered patterns in group effects.
    """
    n_tax = beta_hat.shape[0]
    covariates = (
        x.columns
        if isinstance(x, pd.DataFrame)
        else [f"V{i}" for i in range(x.shape[1])]
    )
    group_ind = np.array([group in c and ":" not in c for c in covariates])
    n_group = int(np.sum(group_ind))

    beta_hat_sub = beta_hat[:, group_ind]
    var_hat_sub = var_hat[:, group_ind]
    vcov_hat_sub = [v[np.ix_(group_ind, group_ind)] for v in vcov_hat]

    # Default contrast: all increasing pattern
    if trend_contrast is None:
        # Create default contrast matrix for monotone increasing trend
        C = np.zeros((n_group - 1, n_group))
        for i in range(n_group - 1):
            C[i, i] = -1
            C[i, i + 1] = 1
        trend_contrast = {"increasing": C}
        trend_node = {"increasing": n_group - 1}

    n_trend = len(trend_contrast)
    trend_names = list(trend_contrast.keys())

    # Constrained estimation for each taxon and each trend pattern
    beta_hat_opt_all = np.zeros((n_tax, n_group * n_trend))
    l_vals = np.zeros((n_tax, n_trend))

    for i in range(n_tax):
        for t_idx, (tname, contrast) in enumerate(trend_contrast.items()):
            C = np.asarray(contrast)
            beta_opt = _constrain_est(beta_hat_sub[i], vcov_hat_sub[i], C)
            start_col = t_idx * n_group
            beta_hat_opt_all[i, start_col : start_col + n_group] = beta_opt

            # Compute l_infinity norm
            node = trend_node[tname]
            l_vals[i, t_idx] = _l_infty(beta_opt, node)

    # W_trend = max l_infinity across patterns
    W_trend = np.max(l_vals, axis=1)
    opt_trend_idx = np.argmax(l_vals, axis=1)

    # Select the optimal trend's beta for each taxon
    beta_hat_trend = np.zeros((n_tax, n_group))
    for i in range(n_tax):
        t_idx = opt_trend_idx[i]
        start_col = t_idx * n_group
        beta_hat_trend[i] = beta_hat_opt_all[i, start_col : start_col + n_group]

    # Bootstrap null distribution
    rng = np.random.default_rng(42)
    W_trend_null = np.zeros((n_tax, trend_B))

    ident_mat = np.eye(n_group)
    var_hat_sub_dup = var_hat_sub  # (n_tax, n_group)

    for b in range(trend_B):
        # Generate null beta from N(0, I)
        beta_null = rng.standard_normal(size=(n_tax, n_group))

        # Constrained estimation under null
        l_null = np.zeros((n_tax, n_trend))
        for i in range(n_tax):
            for t_idx, (tname, contrast) in enumerate(trend_contrast.items()):
                C = np.asarray(contrast)
                beta_null_opt = _constrain_est(beta_null[i], ident_mat, C)
                # Scale by sqrt of variance
                beta_null_opt *= np.sqrt(np.maximum(var_hat_sub_dup[i], 0))
                node = trend_node[tname]
                l_null[i, t_idx] = _l_infty(beta_null_opt, node)

        W_trend_null[:, b] = np.max(l_null, axis=1)

    # P-values from bootstrap
    p_trend = np.mean(W_trend_null > W_trend[:, np.newaxis], axis=1)

    # Multiple testing correction
    func = _check_p_adjust(p_adj_method)
    q_trend = func(p_trend)
    q_trend = np.where(np.isnan(q_trend), 1.0, q_trend)
    diff_trend = q_trend <= alpha

    return {
        "beta": beta_hat_trend,
        "se": np.sqrt(np.maximum(var_hat_sub, 0)),
        "W": W_trend,
        "p_val": p_trend,
        "q_val": q_trend,
        "diff_abn": diff_trend,
    }
