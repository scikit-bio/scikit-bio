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

from typing import Optional
from collections.abc import Mapping
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import norm, chi2, t
from scipy.optimize import minimize
from patsy import dmatrix

from skbio.table._tabular import _ingest_table
from skbio.stats.composition import clr
from ._base import _check_composition
from ._utils import _check_metadata, _check_p_adjust, _type_cast_to_float

# Convergence tolerance for Nelder-Mead optimization in _estimate_bias_em selected
# according to the benchmark. SciPy's default is 1e-4, which would cause different
# results to R version.
NM_TOL = 1e-8


class ANCOMBCResult:
    """Results for ANCOM-BC and ANCOM-BC2 analyses.

    This class contains the primary differential abundance results. Post-hoc analyses
    (global test, multi-group comparisons, sensitivity analysis, etc.) are available as
    methods that compute on-demand using stored intermediate data.

    Attributes
    ----------
    res : pd.DataFrame
        Primary results with (FeatureID, Covariate) multi-index. Columns are: Log2(FC),
        SE, W, pvalue, qvalue, Signif.
    method : {"ANCOM-BC", "ANCOM-BC2"}
        Differential abundance method used for the analysis.

    Methods
    -------
    global_test
        Global test for differential abundance across >= 3 groups.
    dunnett_test
        Dunnett's test: each group vs. reference, with mdFDR correction.
    pairwise_test
        Pairwise directional test between all group pairs, with mdFDR.
    trend_test
        Trend test for ordered patterns in group effects.
    sensitivity_analysis
        Pseudo-count sensitivity analysis for robustness assessment.

    See Also
    --------
    ancombc : ANCOM-BC function.
    ancombc2 : ANCOM-BC2 function.
    """

    _private_defaults = {
        "_dmat": None,
        "_beta_hat": None,
        "_var_hat": None,
        "_vcov_hat": None,
        "_dof": None,
        "_tax_name": None,
        "_fix_eff": None,
        "_alpha": 0.05,
        "_p_adjust": "holm",
        "_s0_perc": 0.05,
        "_max_iter": 100,
        "_tol": 1e-5,
        "_pseudo": 0,
        "_O1": None,
        "_O2": None,
    }

    def __init__(
        self,
        res: pd.DataFrame,
        method: str,
        **kwargs,
    ):
        unexpected = set(kwargs).difference(self._private_defaults)
        if unexpected:
            names = ", ".join(sorted(unexpected))
            raise TypeError(f"Unexpected ANCOMBCResult argument(s): {names}")
        if method not in {"ANCOM-BC", "ANCOM-BC2"}:
            raise ValueError("`method` must be either 'ANCOM-BC' or 'ANCOM-BC2'.")

        self.res = res
        self._method = method
        for name, default in self._private_defaults.items():
            setattr(self, name, kwargs.get(name, default))

    @property
    def res(self) -> pd.DataFrame:
        return self._res

    @res.setter
    def res(self, value: pd.DataFrame):
        self._res = value

    @property
    def method(self) -> str:
        return self._method

    def __getitem__(self, key):
        return getattr(self, key)

    def keys(self):
        """Return a list of attribute names that are not private and are not None."""
        return [
            name
            for name in (
                "res",
                "method",
            )
            if getattr(self, name) is not None
        ]

    # formated output
    def __repr__(self):
        n_taxa = len(self.res.index.get_level_values("FeatureID").unique())
        n_cov = len(self.res.index.get_level_values("Covariate").unique())
        n_signif = int(self.res["Signif"].sum())
        return (
            f"ANCOMBCResult(method={self.method!r}, "
            f"n_taxa={n_taxa}, n_covariates={n_cov}, "
            f"n_signif={n_signif})"
        )

    def _resolve_post_hoc_params(self, alpha, p_adjust):
        if alpha == "inherit":
            alpha = self._alpha
        if p_adjust == "inherit":
            p_adjust = self._p_adjust
        return alpha, p_adjust

    def global_test(
        self, group: str, alpha: float | str = "inherit", p_adjust: str = "inherit"
    ) -> pd.DataFrame:
        """Perform global test for differential abundance across groups.

        The global test identifies features that are differentially abundant
        between at least two groups across three or more groups.

        Parameters
        ----------
        group : str
            Metadata column defining sample groups.
        alpha : float or "inherit", optional
            Significance level, or the value used by :func:`ancombc` or
            :func:`ancombc2`. Default is "inherit".
        p_adjust : str, optional
            Multiple-testing correction method, or "inherit" to use the value
            supplied upstream. Default is "inherit".

        Returns
        -------
        pd.DataFrame
            DataFrame indexed by FeatureID with columns: W, pvalue, qvalue,
            Signif.

        Raises
        ------
        ValueError
            If no group variable was specified.

        Notes
        -----
        """
        alpha, p_adjust = self._resolve_post_hoc_params(alpha, p_adjust)
        if self.method == "ANCOM-BC":
            W_g, pval, qval, reject = _global_test(
                self._dmat,
                group,
                self._beta_hat,
                self._vcov_hat,
                alpha,
                p_adjust,
            )
            result = pd.DataFrame(
                {
                    "W": W_g,
                    "pvalue": pval,
                    "qvalue": qval,
                    "Signif": reject,
                },
                index=self.res.index.get_level_values("FeatureID").unique(),
            )
            result.index.name = "FeatureID"
        else:
            raw = _ancombc_global_F(
                dmat=self._dmat,
                group=group,
                beta_hat=self._beta_hat,
                vcov_hat=self._vcov_hat,
                dof=self._dof,
                p_adj_method=p_adjust,
                alpha=alpha,
            )
            result = raw.copy()
            result.index = self._tax_name
            result = result.rename(
                columns={
                    "p_val": "pvalue",
                    "q_val": "qvalue",
                    "diff_abn": "Signif",
                }
            )
            result["Signif"] = result["Signif"].astype(bool)

        return result

    def dunnett_test(
        self,
        group: str,
        alpha: float | str = "inherit",
        p_adjust: str = "inherit",
        B: int = 100,
    ) -> pd.DataFrame:
        """Perform Dunnett's test (each group vs. reference) with mdFDR.

        Parameters
        ----------
        group : str
            Metadata column defining sample groups.
        alpha : float or "inherit", optional
            Significance level, or the value supplied upstream. Default is
            "inherit".
        p_adjust : str, optional
            P-value adjustment method, or "inherit" to use the value supplied
            upstream. Default is "inherit".
        B : int, optional
            Number of bootstrap iterations. Default is 100.

        Returns
        -------
        pd.DataFrame
            DataFrame with (FeatureID, Comparison) multi-index and columns:
            Log2(FC), SE, W, pvalue, qvalue, Signif.

        """
        alpha, p_adjust = self._resolve_post_hoc_params(alpha, p_adjust)

        raw = _ancombc_dunn(
            dmat=self._dmat,
            group=group,
            beta_hat=self._beta_hat,
            var_hat=self._var_hat,
            dof=self._dof,
            B=B,
            fwer_ctrl_method=p_adjust,
            alpha=alpha,
        )
        comp_names = raw["comp_names"]
        n_comp = len(comp_names)
        n_tax = len(self._tax_name)
        result = pd.DataFrame(
            {
                "FeatureID": [x for x in self._tax_name for _ in range(n_comp)],
                "Comparison": comp_names * n_tax,
                "Log2(FC)": raw["beta"].ravel(),
                "SE": raw["se"].ravel(),
                "W": raw["W"].ravel(),
                "pvalue": raw["p_val"].ravel(),
                "qvalue": raw["q_val"].ravel(),
                "Signif": raw["diff_abn"].ravel(),
            }
        )
        result.set_index(["FeatureID", "Comparison"], inplace=True)
        return result

    def pairwise_test(
        self,
        group: str,
        alpha: float | str = "inherit",
        p_adjust: str = "inherit",
    ) -> pd.DataFrame:
        """Perform pairwise directional test between all group pairs.

        Uses mixed directional FDR (mdFDR) correction via bootstrap.

        Parameters
        ----------
        group : str
            Metadata column defining sample groups.
        alpha : float or "inherit", optional
            Significance level, or the value supplied upstream. Default is
            "inherit".
        p_adjust : str, optional
            P-value adjustment method, or "inherit" to use the value supplied
            upstream. Default is "inherit".

        Returns
        -------
        pd.DataFrame
            DataFrame with (FeatureID, Comparison) multi-index and columns:
            Log2(FC), SE, W, pvalue, qvalue, Signif.

        """
        alpha, p_adjust = self._resolve_post_hoc_params(alpha, p_adjust)

        raw = _ancombc_pair(
            dmat=self._dmat,
            group=group,
            beta_hat=self._beta_hat,
            var_hat=self._var_hat,
            vcov_hat=self._vcov_hat,
            dof=self._dof,
            fwer_ctrl_method=p_adjust,
            alpha=alpha,
        )
        comp_names = raw["comp_names"]
        n_comp = len(comp_names)
        n_tax = len(self._tax_name)
        result = pd.DataFrame(
            {
                "FeatureID": [x for x in self._tax_name for _ in range(n_comp)],
                "Comparison": comp_names * n_tax,
                "Log2(FC)": raw["beta"].ravel(),
                "SE": raw["se"].ravel(),
                "W": raw["W"].ravel(),
                "pvalue": raw["p_val"].ravel(),
                "qvalue": raw["q_val"].ravel(),
                "Signif": raw["diff_abn"].ravel(),
            }
        )
        result.set_index(["FeatureID", "Comparison"], inplace=True)
        return result

    def trend_test(
        self,
        group: str,
        alpha: float | str = "inherit",
        p_adjust: str = "inherit",
        trend_contrast=None,
        trend_node=None,
        trend_B: int = 100,
    ) -> pd.DataFrame:
        """Perform trend test for ordered patterns in group effects.

        Uses constrained optimization to test monotone increasing/decreasing
        patterns in group-level effects.

        Parameters
        ----------
        group : str
            Metadata column defining sample groups.
        alpha : float or "inherit", optional
            Significance level, or the value supplied upstream. Default is
            "inherit".
        p_adjust : str, optional
            Multiple-testing correction method, or "inherit" to use the value
            supplied upstream. Default is "inherit".
        trend_contrast, trend_node : dict, optional
            Trend-test contrast matrices and their node indices.
        trend_B : int, optional
            Number of bootstrap iterations. Default is 100.

        Returns
        -------
        pd.DataFrame
            DataFrame indexed by FeatureID with columns: W, pvalue, qvalue,
            Signif.

        """
        alpha, p_adjust = self._resolve_post_hoc_params(alpha, p_adjust)

        raw = _ancombc_trend(
            dmat=self._dmat,
            group=group,
            beta_hat=self._beta_hat,
            var_hat=self._var_hat,
            vcov_hat=self._vcov_hat,
            p_adj_method=p_adjust,
            alpha=alpha,
            trend_contrast=trend_contrast,
            trend_node=trend_node,
            trend_B=trend_B,
        )
        result = pd.DataFrame(
            {
                "W": raw["W"],
                "pvalue": raw["p_val"],
                "qvalue": raw["q_val"],
                "Signif": raw["diff_abn"],
            },
            index=self._tax_name,
        )
        result.index.name = "FeatureID"
        return result

    def sensitivity_analysis(
        self,
        group: Optional[str] = None,
        global_test: bool = False,
        pairwise: bool = False,
        dunnett: bool = False,
        trend: bool = False,
        alpha: float | str = "inherit",
        p_adjust: str = "inherit",
        mdfdr_B: int = 100,
        trend_contrast=None,
        trend_node=None,
        trend_B: int = 100,
    ) -> dict:
        """Perform pseudo-count sensitivity analysis.

        Re-runs the core estimation with pseudo-count values [0.1, 0.5, 1]
        and compares q-values to assess robustness of significance calls to
        the choice of pseudo-count.

        Set a multi-group test option to True to assess its sensitivity.

        Parameters
        ----------
        group : str, optional
            Metadata column defining sample groups. Required for multi-group
            sensitivity analysis.
        global_test, pairwise, dunnett, trend : bool, optional
            Whether to assess sensitivity for the corresponding multi-group
            test. Defaults are False.
        alpha : float or "inherit", optional
            Significance level, or the value supplied upstream. Default is
            "inherit".
        p_adjust : str, optional
            P-value adjustment method, or "inherit" to use the value supplied
            upstream. Default is "inherit".
        mdfdr_B, trend_B : int, optional
            Bootstrap iterations for the Dunnett and trend tests. Defaults are 100.
        trend_contrast, trend_node : dict, optional
            Trend-test contrast matrices and their node indices.

        Returns
        -------
        dict
            Dictionary with keys:

            - ``passed_ss_prim``: bool ndarray, primary results robustness
            - ``passed_ss_global``: bool ndarray or None
            - ``passed_ss_pair``: bool ndarray or None
            - ``passed_ss_dunn``: bool ndarray or None
            - ``passed_ss_trend``: bool ndarray or None
            - ``ss_tab_prim``: float ndarray, proportion of pseudo-counts where
              q > alpha.

        Raises
        ------
        ValueError
            If this is an ANCOM-BC result (reestimate=False).

        Notes
        -----
        This method recomputes its result on every call.
        """
        if self.method == "ANCOM-BC":
            raise ValueError(
                "sensitivity_analysis() is only available for ANCOM-BC2 results."
            )

        n_tax = len(self._tax_name)
        n_cov = len(self._fix_eff)
        alpha, p_adjust = self._resolve_post_hoc_params(alpha, p_adjust)

        if (global_test or pairwise or dunnett or trend) and group is None:
            raise ValueError("Multi-group sensitivity analysis requires `group`.")

        # Original q_hat (n_tax, n_cov)
        q_hat_orig = self.res["qvalue"].values.reshape(n_tax, n_cov)

        # Re-run core estimation with each pseudo-count
        pseudo_list = [0.1, 0.5, 1]
        ss_list = []
        for pseudo_count in pseudo_list:
            res_pseudo = _ancombc2_estimate(
                data=self._O1,
                aggregated_data=self._O2,
                dmat=self._dmat,
                p_adj_method=p_adjust,
                pseudo=pseudo_count,
                s0_perc=self._s0_perc,
                alpha=alpha,
                max_iter=self._max_iter,
                tol=self._tol,
            )
            # If multi-group tests are requested, run them for this pseudo
            if global_test:
                res_pseudo["res_global"] = _ancombc_global_F(
                    dmat=res_pseudo["dmat"],
                    group=group,
                    beta_hat=res_pseudo["beta_hat"],
                    vcov_hat=res_pseudo["vcov_hat"],
                    dof=res_pseudo["dof"],
                    p_adj_method=p_adjust,
                    alpha=alpha,
                )
            if pairwise:
                res_pseudo["res_pair"] = _ancombc_pair(
                    dmat=res_pseudo["dmat"],
                    group=group,
                    beta_hat=res_pseudo["beta_hat"],
                    var_hat=res_pseudo["var_hat"],
                    vcov_hat=res_pseudo["vcov_hat"],
                    dof=res_pseudo["dof"],
                    fwer_ctrl_method=p_adjust,
                    alpha=alpha,
                )
            if dunnett:
                res_pseudo["res_dunn"] = _ancombc_dunn(
                    dmat=res_pseudo["dmat"],
                    group=group,
                    beta_hat=res_pseudo["beta_hat"],
                    var_hat=res_pseudo["var_hat"],
                    dof=res_pseudo["dof"],
                    B=mdfdr_B,
                    fwer_ctrl_method=p_adjust,
                    alpha=alpha,
                )
            if trend:
                res_pseudo["res_trend"] = _ancombc_trend(
                    dmat=res_pseudo["dmat"],
                    group=group,
                    beta_hat=res_pseudo["beta_hat"],
                    var_hat=res_pseudo["var_hat"],
                    vcov_hat=res_pseudo["vcov_hat"],
                    p_adj_method=p_adjust,
                    alpha=alpha,
                    trend_contrast=trend_contrast,
                    trend_node=trend_node,
                    trend_B=trend_B,
                )
            ss_list.append(res_pseudo)

        # Combine original and pseudo-count results
        all_q = [q_hat_orig] + [r["q_hat"] for r in ss_list]
        q_3d = np.stack(all_q, axis=-1)  # (n_tax, n_cov, n_pseudo)
        ss_tab_prim = np.mean(q_3d > alpha, axis=-1)
        passed_ss_prim = (ss_tab_prim == 0) | (ss_tab_prim == 1)

        # Global test sensitivity
        passed_ss_global = None
        if global_test:
            q_globals = [self.global_test(group, alpha, p_adjust)["qvalue"].values]
            for r in ss_list:
                if r.get("res_global") is not None:
                    q_globals.append(r["res_global"]["q_val"])
                else:
                    q_globals.append(np.ones(n_tax))
            q_global_3d = np.stack(q_globals, axis=-1)
            ss_global = np.mean(q_global_3d > alpha, axis=-1)
            passed_ss_global = (ss_global == 0) | (ss_global == 1)

        # Pairwise sensitivity
        passed_ss_pair = None
        if pairwise:
            q_pair_orig = self.pairwise_test(group, alpha, p_adjust)[
                "qvalue"
            ].values.reshape(n_tax, -1)
            q_pairs = [q_pair_orig]
            for r in ss_list:
                if r.get("res_pair") is not None:
                    q_pairs.append(r["res_pair"]["q_val"])
                else:
                    q_pairs.append(np.ones_like(q_pair_orig))
            q_pair_3d = np.stack(q_pairs, axis=-1)
            ss_pair = np.mean(q_pair_3d > alpha, axis=-1)
            passed_ss_pair = (ss_pair == 0) | (ss_pair == 1)

        # Dunnett sensitivity
        passed_ss_dunn = None
        if dunnett:
            q_dunn_orig = self.dunnett_test(group, alpha, p_adjust, mdfdr_B)[
                "qvalue"
            ].values.reshape(n_tax, -1)
            q_dunns = [q_dunn_orig]
            for r in ss_list:
                if r.get("res_dunn") is not None:
                    q_dunns.append(r["res_dunn"]["q_val"])
                else:
                    q_dunns.append(np.ones_like(q_dunn_orig))
            q_dunn_3d = np.stack(q_dunns, axis=-1)
            ss_dunn = np.mean(q_dunn_3d > alpha, axis=-1)
            passed_ss_dunn = (ss_dunn == 0) | (ss_dunn == 1)

        # Trend sensitivity uses global sensitivity.
        passed_ss_trend = None
        if trend:
            q_trend_orig = self.trend_test(
                group, alpha, p_adjust, trend_contrast, trend_node, trend_B
            )["qvalue"].values
            q_trends = [q_trend_orig] + [r["res_trend"]["q_val"] for r in ss_list]
            q_trend_3d = np.stack(q_trends, axis=-1)
            ss_trend = np.mean(q_trend_3d > alpha, axis=-1)
            passed_ss_trend = (ss_trend == 0) | (ss_trend == 1)

        result = {
            "passed_ss_prim": passed_ss_prim,
            "passed_ss_global": passed_ss_global,
            "passed_ss_pair": passed_ss_pair,
            "passed_ss_dunn": passed_ss_dunn,
            "passed_ss_trend": passed_ss_trend,
            "ss_tab_prim": ss_tab_prim,
        }
        return result


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
        Transposed from SVD output to match downstream convention.
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
    data = np.asarray(data)
    dmat = np.asarray(dmat)
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
    intm_mat = np.einsum("ji,jp,jq->ipq", eps**2, dmat, dmat, optimize=True)
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
    # Note: `norm.logpdf` doesn't have an `out` parameter. To further optimize this,
    # one needs to manually implement the under-the-hood algorithm.
    def func(x, loc, resp):
        log_pdf = norm.logpdf(beta, loc=loc, scale=(var_hat + x) ** 0.5)
        return -np.dot(resp, log_pdf)

    # Optimizer arguments. The Nelder-Mead simplex algorithm is used, which is
    # consistent with the original R code.
    # Note: The Nelder–Mead method doesn't actually enforce bounds during optimization.
    # It merely clips at the bounds. For further consideration.
    # Tight convergence criteria (xatol=1e-12, fatol=1e-12) produce more precise kappa
    # estimates every EM step, leading to smaller accumulated errors.
    args = dict(
        method="Nelder-Mead",
        bounds=[(0, None)],
        options={"xatol": NM_TOL, "fatol": NM_TOL},
    )

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
    r"""Identify features with structural zeros.

    .. versionadded:: 0.7.1

    Structural zeros refer to features that are systematically absent from certain
    sample groups. Consequently, the observed feature frequencies are all zeros, or
    mostly zeros, due to variability in technical factors. This function tests
    whether the proportion of observed zeros is close to zero, which suggests the
    absence of a feature in a given sample group.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing count or proportional abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    metadata : pd.DataFrame or 2-D array_like
        Metadata of the samples. Rows correspond to samples and columns correspond
        to covariates (attributes). Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    grouping : str
        A metadata column name indicating the assignment of samples to groups.
    neg_lb : bool, optional
        Determine whether to use negative lower bound when calculating sample
        proportions. Default is False. Generally, it is recommended to set it as True
        when the sample size per group is relatively large.

    Returns
    -------
    pd.DataFrame of bool of shape (n_features, n_groups)
        A table indicating whether each feature (row) is a structural zero in each
        group (column) (True: structural zero, False: not structural zero).

    Notes
    -----
    The structural zero test was initially proposed and implemented in the ANCOM-II
    method [1]_. It was adopted to the ANCOM-BC method [2]_ as a recommended method to
    complement test results. See :func:`ancombc` for how to use this function along
    with the ANCOM-BC test. Nevertheless, this function is generally useful with or
    without explicit statistical tests of feature abundances.

    A feature found to be a structural zero in a group should be automatically
    considered as differentially (less) abundant compared with other groups in which
    this feature is not a structural zero. Meanwhile, this feature should be excluded
    from subsequent analyses that involves this group. If a feature is identified as a
    structural zero in all groups, this feature should be removed entirely from
    downstream analyses.

    Note that the structural zero test should be applied to the original table before
    adding a pseudocount, which will otherwise mask all zeros and invalidate this test.

    References
    ----------
    .. [1] Kaul, A., Mandal, S., Davidov, O., & Peddada, S. D. (2017). Analysis of
       microbiome data in the presence of excess zeros. Frontiers in Microbiology, 8,
       2114.

    .. [2] Lin, H. and Peddada, S.D., 2020. Analysis of compositions of microbiomes
       with bias correction. Nature Communications, 11(1), p.3514.

    Examples
    --------
    >>> from skbio.stats.composition import struc_zero
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

    The ``struc_zero`` function will identify features with structural zeros. Features
    that are identified as structural zeros in given groups should not be used in
    further analyses such as ``ancombc`` and  ``dirmult_ttest``.

    Setting ``neg_lb=True`` declares that the true prevalence of a feature in a group
    is not significantly different from zero.

    >>> result = struc_zero(table, metadata, grouping='grouping', neg_lb=True)
    >>> result
        placebo  treatment
    f0    False      False
    f1    False      False
    f2    False       True
    f3    False      False
    f4     True      False
    f5    False       True

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
    return pd.DataFrame(zero_idx.T, index=features, columns=unique_groups)


def _global_test(dmat, grouping, beta_hat, vcov_hat, alpha=0.05, p_adjust="holm"):
    """Perform ANCOM-BC global test.

    The global test is to determine features that are differentially abundant between
    at least 2 sample groups across 3 or more groups.

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
    W_global : ndarray of float of shape (n_features,)
        W-statistics of global test.
    pval : ndarray of float of shape (n_features,)
        p-values of global test.
    qval : ndarray of float of shape (n_features,)
        Adjusted p-values of global test.
    reject : ndarray of bool of shape (n_features,)
        If the variable is differentially abundant.

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

    # for each feature, calculate test statistics W by the following formula:
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


def _ancombc2_estimate(
    data,
    aggregated_data,
    dmat,
    p_adj_method="holm",
    pseudo=0,
    s0_perc=0.05,
    alpha=0.05,
    max_iter=100,
    tol=1e-5,
):
    """ANCOM-BC2 core estimation"""
    # data preprocessing
    O1 = data + pseudo
    O2 = aggregated_data + pseudo

    n_tax = O2.shape[1]
    tax_name = O2.columns.tolist()

    fix_eff = dmat.design_info.column_names
    n_fix_eff = len(fix_eff)

    # CLR transformation
    # Per-taxon centering: y = log(x) - mean(log(x)) per taxon (row).
    O1_vals = O1.values.astype(float)
    y1 = clr(O1_vals, axis=0)

    # Initial estimates via SVD-based OLS (_estimate_params)
    var_hat1, beta1, _, _ = _estimate_params(data=y1, dmat=dmat)
    beta1 = beta1.T  # (n_feats, n_covariates)

    # EM bias correction
    bias1 = np.empty((n_fix_eff, 3))
    for i in range(n_fix_eff):
        beta_col = beta1[:, i]
        var_col = var_hat1[:, i]
        # Remove NaN pairs
        mask = ~(np.isnan(beta_col) | np.isnan(var_col))
        if mask.sum() > 0:
            bias1[i] = _estimate_bias_em(
                beta_col[mask], var_col[mask], tol=tol, max_iter=max_iter
            )
        else:
            bias1[i] = [0.0, 0.0, 0.0]

    delta_em = bias1[:, 0]
    delta_wls = bias1[:, 1]
    var_delta = bias1[:, 2]

    # Correct coefficients and compute sampling fractions
    beta1_corrected = beta1 - delta_em[np.newaxis, :]

    # Compute theta_hat (sampling fractions)
    resid_all = y1 - (dmat @ beta1_corrected.T)  # (n_samp, n_feats)
    theta_hat_arr = np.nanmean(resid_all, axis=1)  # mean over features/taxa (columns)

    # Handle NaN in theta (samples with all-NaN residuals)
    nan_theta = np.isnan(theta_hat_arr)
    if np.any(nan_theta):
        theta_hat_arr[nan_theta] = 0.0

    # Final estimates with bias-corrected theta
    O2_vals = O2.values.astype(float)
    y2 = clr(O2_vals, axis=0)

    # Subtract theta (sampling fractions) from data (y) before estimation.
    # _estimate_params returns beta as (n_covariates, n_features),
    # transpose to (n_features, n_covariates) for downstream consistency.
    y2_crt = y2 - theta_hat_arr[:, np.newaxis]  # (n_samples, n_feats)
    var_hat, beta_hat, _, beta_covmat2 = _estimate_params(data=y2_crt, dmat=dmat)
    beta_hat = beta_hat.T  # (n_feats, n_covariates)
    # Convert beta_covmat to list format for compatibility
    vcov_hat = [beta_covmat2[i] for i in range(n_tax)]
    # Compute degrees of freedom (per-taxon, based on valid samples)
    dof = np.full((n_tax, n_fix_eff), np.nan)
    for i in range(n_tax):
        n_valid = np.sum(~np.isnan(y2_crt[:, i]))
        if n_valid > n_fix_eff:
            dof[i, :] = n_valid - n_fix_eff

    # Variance adjustment
    # Account for variance of delta
    var_hat = (
        var_hat
        + var_delta[np.newaxis, :]
        + 2 * np.sqrt(np.abs(var_hat * var_delta[np.newaxis, :]))
    )

    # SAM-like fudge factor
    if s0_perc is None:
        s02 = np.zeros(n_fix_eff)
    else:
        s02 = np.nanquantile(var_hat, s0_perc, axis=0)

    var_hat = var_hat + s02[np.newaxis, :]
    var_hat[np.isnan(beta_hat)] = np.nan

    se_hat = np.sqrt(np.maximum(var_hat, 0))
    se_hat[np.isnan(var_hat)] = np.nan

    # Update vcov_hat with adjusted variances
    for i in range(n_tax):
        np.fill_diagonal(vcov_hat[i], var_hat[i])

    # Step 7: Primary results
    W = np.where(np.isnan(se_hat), np.nan, beta_hat / se_hat)

    # P-values from t-test or Z-test
    if dof is not None:
        p_hat = 2 * t.sf(np.abs(W), df=dof)
    else:
        p_hat = 2 * norm.sf(np.abs(W))
    p_hat = np.where(np.isnan(p_hat), 1.0, p_hat)

    # Multiple testing correction per covariate
    func = _check_p_adjust(p_adj_method)
    q_hat = np.zeros_like(p_hat)
    for j in range(n_fix_eff):
        q_hat[:, j] = func(p_hat[:, j])

    diff_abn = q_hat <= alpha

    # Bias-corrected log table
    y_bias_crt = pd.DataFrame(
        y2 - theta_hat_arr[:, np.newaxis], index=O2.index, columns=O2.columns
    )

    return {
        "bias_correct_log_table": y_bias_crt,
        "beta_hat": beta_hat,
        "se_hat": se_hat,
        "W": W,
        "p_hat": p_hat,
        "q_hat": q_hat,
        "diff_abn": diff_abn,
        "vcov_hat": vcov_hat,
        "var_hat": var_hat,
        "dof": dof,
        "fix_eff": fix_eff,
        "tax_name": tax_name,
        "dmat": dmat,
        "y2": y2,
        "res_global": None,
        "res_pair": None,
        "res_dunn": None,
        "res_trend": None,
    }


def _aggregate_features(table, aggregator, has_feature_ids):
    """Aggregate feature columns according to user-supplied aggregate IDs."""
    if aggregator is None:
        return table

    features = table.columns
    if callable(aggregator):
        if not has_feature_ids:
            raise ValueError("A callable aggregator requires named features.")
        aggregate_ids = [aggregator(feature) for feature in features]
    elif isinstance(aggregator, Mapping) or isinstance(aggregator, pd.Series):
        if not has_feature_ids:
            raise ValueError("A mapping aggregator requires named features.")
        try:
            aggregate_ids = [aggregator[feature] for feature in features]
        except KeyError as error:
            raise ValueError(
                f"Aggregator does not define feature {error.args[0]!r}."
            ) from error
    else:
        aggregate_ids = np.asarray(aggregator, dtype=object)
        if aggregate_ids.ndim != 1 or len(aggregate_ids) != len(features):
            raise ValueError(
                "A sequence aggregator must be one-dimensional and have one entry "
                "per feature."
            )

    aggregate_ids = np.asarray(aggregate_ids, dtype=object)
    if pd.isna(aggregate_ids).any():
        raise ValueError("Aggregator must assign every feature an aggregate ID.")

    return table.T.groupby(aggregate_ids, sort=False).sum().T


def _ancombc(
    table,
    metadata,
    reestimate=False,
    formula=None,
    aggregator=None,
    p_adjust="holm",
    pseudo=0,
    s0_perc=0.05,
    alpha=0.05,
    max_iter=100,
    tol=1e-5,
):
    """Private ANCOM-BC/BC2 core function."""
    # Note: A pseudocount should have been added to the table by the user prior to
    # calling this function.
    matrix, samples, features = _ingest_table(table)
    # NOTE: ANCOM-BC does not handle zeros in the input table. The user should have
    # added a pseudocount. ANCOM-BC2 should be able to handle zeros.
    # TODO: Add zero-handling in ANCOM-BC2.
    _check_composition(np, matrix, nozero=not reestimate)
    n_feats = matrix.shape[1]
    has_feature_ids = features is not None
    if features is None:
        features = np.arange(n_feats)

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
        raise ValueError(f"`alpha`={alpha} is not within 0 and 1.")

    if not reestimate:
        # Log-transform count matrix.
        matrix = np.log(matrix)

        # Estimate initial model parameters.
        var_hat, beta, _, vcov_hat = _estimate_params(
            matrix, dmat
        )  # TODO: use a helper function for initial estimates

        # Estimate and correct for sampling bias via expectation-maximization (EM).
        # beta: (n_covariates, n_features); iterate over covariates (rows).
        bias = np.empty((n_covars, 3))
        for i in range(n_covars):
            bias[i] = _estimate_bias_em(
                beta[i], var_hat[:, i], tol=tol, max_iter=max_iter
            )
        delta_em = bias[:, 0]
        delta_wls = bias[:, 1]

        # Correct coefficients (logFC) according to estimated bias.
        # beta: (n_covariates, n_features); transpose to (n_features, n_covariates)
        # then subtract delta_em broadcast over columns.
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
        res.set_index(["FeatureID", "Covariate"], inplace=True)

        return ANCOMBCResult(
            res=res,
            method="ANCOM-BC",
            _dmat=dmat,
            _beta_hat=beta_hat,
            _var_hat=var_hat,
            _vcov_hat=vcov_hat,
            _tax_name=features,
            _alpha=alpha,
            _p_adjust=p_adjust,
        )

    else:
        # Data preprocessing (no filtering. user pre-filters)
        O1 = pd.DataFrame(matrix, index=samples, columns=features)
        O2 = _aggregate_features(O1, aggregator, has_feature_ids)

        # Step 3: Run ANCOMBC2 core estimation
        res_main = _ancombc2_estimate(
            data=O1,
            aggregated_data=O2,
            dmat=dmat,
            p_adj_method=p_adjust,
            pseudo=pseudo,
            s0_perc=s0_perc,
            alpha=alpha,
            max_iter=max_iter,
            tol=tol,
        )

        # Format primary results into a DataFrame with MultiIndex (FeatureID, Covariate)
        tax_name = res_main["tax_name"]
        fix_eff = res_main["fix_eff"]
        n_tax = len(tax_name)
        n_cov = len(fix_eff)

        beta_hat = res_main["beta_hat"]
        se_hat = res_main["se_hat"]
        W = res_main["W"]
        p_hat = res_main["p_hat"]
        q_hat = res_main["q_hat"]
        diff_abn = res_main["diff_abn"]

        res_dict = {
            "FeatureID": [x for x in tax_name for _ in range(n_cov)],
            "Covariate": list(fix_eff) * n_tax,
            "Log2(FC)": beta_hat.ravel(),
            "SE": se_hat.ravel(),
            "W": W.ravel(),
            "pvalue": p_hat.ravel(),
            "qvalue": q_hat.ravel(),
        }
        res_dict["Signif"] = diff_abn.ravel()
        res = pd.DataFrame(res_dict)
        res.set_index(["FeatureID", "Covariate"], inplace=True)

        return ANCOMBCResult(
            res=res,
            method="ANCOM-BC2",
            _dmat=res_main["dmat"],
            _beta_hat=res_main["beta_hat"],
            _var_hat=res_main["var_hat"],
            _vcov_hat=res_main["vcov_hat"],
            _dof=res_main["dof"],
            _tax_name=tax_name,
            _fix_eff=fix_eff,
            _alpha=alpha,
            _p_adjust=p_adjust,
            _s0_perc=s0_perc,
            _max_iter=max_iter,
            _tol=tol,
            _pseudo=pseudo,
            _O1=O1,
            _O2=O2,
        )


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
        A matrix containing strictly positive count or proportional abundance data of
        the samples. See :ref:`supported formats <table_like>`.

    metadata : pd.DataFrame or 2-D array_like
        Metadata of the samples. Rows correspond to samples and columns correspond
        to covariates (attributes). Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    formula : str or generic Formula object
        A formula defining the model using factors included in the metadata columns.
        Refer to `Patsy's documentation
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
        Method to correct *p*-values for multiple comparisons. Options are
        Holm-Bonferroni ("holm" or "holm-bonferroni") (default), Benjamini-Hochberg
        ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported by
        statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.

    Returns
    -------
    :class:`ANCOMBCResult`
        Result object with primary results and post-hoc analysis methods.

    See Also
    --------
    ancombc2 : ANCOM-BC2 with multi-group tests and sensitivity analysis.
    struc_zero : Standalone structural zero detection.

    References
    ----------
    .. [1] Lin, H. and Peddada, S.D., 2020. Analysis of compositions of microbiomes
       with bias correction. Nature Communications, 11(1), p.3514.

    Examples
    --------
    >>> from skbio.stats.composition import ancombc
    >>> import pandas as pd

    **A basic example**

    Let's create a data table with six samples and seven features (e.g., these may be
    microbial taxa):

    >>> samples = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    >>> features = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7']
    >>> table = pd.DataFrame(
    ...     [[ 2,  0,  4,  7,  0,  0,  1],
    ...      [ 1,  1,  0,  6,  5,  1, 10],
    ...      [ 3,  0,  2,  9,  6,  0,  1],
    ...      [ 0, 12,  1,  2,  0,  3,  2],
    ...      [ 2,  2, 37,  0,  0,  7,  3],
    ...      [10,  1,  0,  0,  4,  4,  3]],
    ...     index=samples, columns=features)

    Then create a sampling grouping vector. In this example, there is a "healthy" group
    and a "sick" group.

    >>> grouping = ['healthy', 'healthy', 'healthy', 'sick', 'sick', 'sick']
    >>> metadata = pd.Series(grouping, index=samples, name='status').to_frame()

    Now run ``ancombc`` to determine if there are any features that are significantly
    different in abundance between the healthy and the sick groups. Note that a
    pseudocount of 1 is manually added to the table to remove zero values.

    >>> result = ancombc(table + 1, metadata, formula='status').res
    >>> result.round(3)
                              Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F1        Intercept         -0.045  0.218 -0.207   0.836   1.000   False
              status[T.sick]     0.126  0.589  0.214   0.831   1.000   False
    F2        Intercept         -0.873  0.137 -6.381   0.000   0.000    True
              status[T.sick]     1.241  0.538  2.307   0.021   0.105   False
    F3        Intercept         -0.202  0.475 -0.425   0.671   1.000   False
              status[T.sick]     0.561  0.964  0.582   0.561   1.000   False
    F4        Intercept          1.005  0.136  7.399   0.000   0.000    True
              status[T.sick]    -1.723  0.392 -4.399   0.000   0.000    True
    F5        Intercept          0.141  0.422  0.335   0.738   1.000   False
              status[T.sick]    -0.690  0.625 -1.104   0.270   1.000   False
    F6        Intercept         -0.873  0.137 -6.381   0.000   0.000    True
              status[T.sick]     1.480  0.160  9.255   0.000   0.000    True
    F7        Intercept          0.157  0.401  0.391   0.695   1.000   False
              status[T.sick]     0.049  0.405  0.121   0.904   1.000   False

    The covariate "status[T.sick]" stands for the effect of the "sick" group relative
    to the reference group, "healthy" (the first group in alphabetical order is
    automatically selected as the reference group). "Log2(FC)" represents the
    :func:`clr`-transformed fold change of abundance (positive/negative: more/less
    abundant in "sick" than in "healthy", respectively). A "True" in the "Signif"
    column indicates a significantly differentially abundant feature-covariate pair.
    This example shows that two features differ by healthy/sick status.

    >>> result.query('Covariate != "Intercept" & Signif == True').round(3)
                              Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F4        status[T.sick]    -1.723  0.392 -4.399     0.0     0.0    True
    F6        status[T.sick]     1.480  0.160  9.255     0.0     0.0    True

    **An advanced example**

    Now we will create a complex dataset with 15 samples grouped into three disease
    status: "mild", "moderate" and "severe", plus age as a confounder. Our goal is
    to identify features that are differentially abundant between disease status.

    >>> samples = [f'S{i}' for i in range(1, 16)]
    >>> features = [f'F{i}' for i in range(1, 9)]
    >>> data = [[ 2,  0,  7,  0,  0,  2,  3,  2],
    ...         [ 0,  2,  0,  1,  1,  2,  0,  0],
    ...         [ 3,  1,  0,  9,  0,  1,  1,  0],
    ...         [ 2,  0,  1,  8,  1,  0, 11, 46],
    ...         [ 1,  0,  1,  1,  0,  0,  2,  2],
    ...         [ 0,  3, 22,  1,  1,  1,  3,  0],
    ...         [ 1,  7, 16,  1,  0,  0,  2,  2],
    ...         [ 0,  5,  6,  1,  2,  1,  0,  1],
    ...         [ 1,  7,  0,  2,  1,  0,  3,  2],
    ...         [ 0,  6,  4,  2,  0,  2,  2,  1],
    ...         [ 3, 13,  7,  0,  0,  0,  3,  9],
    ...         [ 1,  8,  5,  1,  0,  0,  0,  0],
    ...         [ 0,  5, 14,  1,  0,  1,  0,  1],
    ...         [ 5, 26,  3,  2,  0,  3,  1,  3],
    ...         [ 0, 18,  7,  0,  0,  3,  1,  0]]
    >>> table = pd.DataFrame(data, index=samples, columns=features)
    >>> status = ['mild'] * 5 + ['moderate'] * 5 + ['severe'] * 5
    >>> age = [39, 19, 20, 31, 15, 37, 27, 47, 26, 23, 39, 48, 46, 33, 36]
    >>> metadata = pd.DataFrame({'status': status, 'age': age}, index=samples)

    Run ``ancombc``. This time, we specify two factors: "status" and "age" in the
    formula, such that the function will test the individual effects of each factor
    while controlling for the other. Additionally, we instruct the function to perform
    a *global test* on "status", which identifies features that are differentially
    abundant between at least two status.

    >>> res = ancombc(table + 1, metadata, formula='status + age')
    >>> res_main = res.res
    >>> res_main.round(3)
                                  Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F1        Intercept             -0.006  0.346 -0.016   0.987   1.000   False
              status[T.moderate]    -0.337  0.245 -1.376   0.169   1.000   False
              status[T.severe]       0.278  0.350  0.792   0.428   1.000   False
              age                    0.001  0.011  0.117   0.907   1.000   False
    F2        Intercept              0.072  0.446  0.160   0.873   1.000   False
              status[T.moderate]     1.906  0.259  7.371   0.000   0.000    True
              status[T.severe]       2.935  0.304  9.649   0.000   0.000    True
              age                   -0.022  0.013 -1.697   0.090   0.538   False
    F3        Intercept             -1.690  0.570 -2.967   0.003   0.021    True
              status[T.moderate]     1.010  0.578  1.747   0.081   0.565   False
              status[T.severe]       0.717  0.517  1.388   0.165   1.000   False
              age                    0.063  0.022  2.871   0.004   0.033    True
    F4        Intercept              0.439  0.483  0.910   0.363   1.000   False
              status[T.moderate]    -0.046  0.405 -0.113   0.910   1.000   False
              status[T.severe]      -0.244  0.536 -0.455   0.649   1.000   False
              age                   -0.004  0.019 -0.199   0.842   1.000   False
    F5        Intercept             -1.227  0.393 -3.125   0.002   0.014    True
              status[T.moderate]     0.274  0.293  0.934   0.350   1.000   False
              status[T.severe]      -0.323  0.320 -1.011   0.312   1.000   False
              age                    0.027  0.014  2.001   0.045   0.318   False
    F6        Intercept             -0.439  0.440 -0.999   0.318   1.000   False
              status[T.moderate]     0.114  0.432  0.264   0.792   1.000   False
              status[T.severe]       0.375  0.544  0.691   0.490   1.000   False
              age                    0.008  0.016  0.464   0.643   1.000   False
    F7        Intercept              0.201  0.472  0.426   0.670   1.000   False
              status[T.moderate]     0.082  0.302  0.270   0.787   1.000   False
              status[T.severe]      -0.264  0.401 -0.658   0.510   1.000   False
              age                    0.004  0.016  0.256   0.798   1.000   False
    F8        Intercept             -0.168  0.526 -0.319   0.750   1.000   False
              status[T.moderate]    -0.402  0.584 -0.688   0.491   1.000   False
              status[T.severe]      -0.299  0.728 -0.411   0.681   1.000   False
              age                    0.022  0.018  1.222   0.222   1.000   False

    We found that feature "F2" is significantly differentially (more) abundant in
    "moderate" and "severe" groups compared with "mild", which serves as the reference
    group. Additionally, feature "F3" is separately correlated with age.

    >>> res_main.query('Covariate != "Intercept" & Signif == True').round(3)
                                  Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F2        status[T.moderate]     1.906  0.259  7.371   0.000   0.000    True
              status[T.severe]       2.935  0.304  9.649   0.000   0.000    True
    F3        age                    0.063  0.022  2.871   0.004   0.033    True

    The global test result suggests that "F2" and "F4" are differentially abundant
    between two of the three groups (though it doesn't tell which groups).

    >>> res_global = res.global_test(group='status')
    >>> res_global.round(3)
                    W  pvalue  qvalue  Signif
    FeatureID
    F1          1.855   0.791   1.000   False
    F2         80.771   0.000   0.000    True
    F3          2.925   0.463   1.000   False
    F4         -0.093   0.000   0.000    True
    F5          1.068   0.827   1.000   False
    F6          0.121   0.117   0.704   False
    F7          0.220   0.208   1.000   False
    F8          0.485   0.430   1.000   False

    **Structural zero test**

    The structural zero test identifies features that are systematically absent from
    certain sample groups. This test is an option of the R command ``ancombc``. In
    scikit-bio, :func:`struc_zero` is a standalone function, as it is generally useful
    with or without ANCOM-BC.

    >>> from skbio.stats.composition import struc_zero
    >>> res_zero = struc_zero(table, metadata, 'status')
    >>> res_zero
         mild  moderate  severe
    F1  False     False   False
    F2  False     False   False
    F3  False     False   False
    F4  False     False   False
    F5  False     False    True
    F6  False     False   False
    F7  False     False   False
    F8  False     False   False

    The result reveals that feature "F5" is a structural zero in the "severe" groups,
    as all or most of its values are zero. Although the ANCOM-BC test itself didn't
    identify "F5", we should consider it as differentially (less) abundant than in the
    other two groups.

    We can use this additional information to update the global test result. The rule
    is that a feature should be considered as globally differentially abundant if it is
    a structural zero in at least one but not all groups.

    >>> signif_global = res_global['Signif'] | (
    ...     ~res_zero.all(axis=1) & res_zero.any(axis=1))
    >>> signif_global
    FeatureID
    F1    False
    F2     True
    F3    False
    F4     True
    F5     True
    F6    False
    F7    False
    F8    False
    dtype: bool

    The main ANCOM-BC result can also be updated with structural zeros.

    >>> signif_main = res_main.query(
    ...     'Covariate.str.startswith("status[T.")')['Signif'].unstack()
    >>> signif_main.columns = signif_main.columns.str.removeprefix(
    ...     'status[T.').str.removesuffix(']')
    >>> signif_zero = res_zero.loc[signif_main.index, signif_main.columns]
    >>> signif_main |= signif_zero
    >>> signif_main.columns.name = None
    >>> signif_main
               moderate  severe
    FeatureID
    F1            False   False
    F2             True    True
    F3            False   False
    F4            False   False
    F5            False    True
    F6            False   False
    F7            False   False
    F8            False   False

    """
    return _ancombc(
        table=table,
        metadata=metadata,
        reestimate=False,
        formula=formula,
        max_iter=max_iter,
        tol=tol,
        alpha=alpha,
        p_adjust=p_adjust,
    )


def ancombc2(
    table,
    metadata,
    formula,
    aggregator=None,
    p_adjust="holm",
    pseudo=0,
    s0_perc=0.05,
    alpha=0.05,
    max_iter=100,
    tol=1e-5,
):
    r"""Perform differential abundance test using ANCOM-BC2.

    Analysis of Compositions of Microbiomes with Bias Correction 2
    (ANCOM-BC2) [1]_ is a differential abundance testing method that extends
    ANCOM-BC with several improvements:

    - **SAM-like fudge factor** (``s0_perc``): Adds a small constant to
      the denominator of the test statistic to stabilize inference for
      rare taxa with extremely small standard errors.
    - **Pseudo-count sensitivity analysis**: Assesses robustness of results
      to the choice of pseudo-count added to zero counts.
    - **Multi-group comparisons**: Supports global test, pairwise
      directional test, Dunnett's type of test, and trend test.
    - **Mixed directional FDR** (mdFDR): Controls false discovery rate
      in multi-group settings using a bootstrap-based approach.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Count or proportional abundance data. See :ref:`supported formats
        <table_like>`.
    metadata : pd.DataFrame or 2-D array_like
        Sample metadata.
    formula : str
        Formula defining the model using factors in metadata columns.
    aggregator : callable, mapping, or 1-D array_like, optional
        Rule for aggregating features before final regression. A callable maps
        each feature ID to an aggregate ID. A mapping, including a
        :class:`pandas.Series`, maps feature IDs to aggregate IDs. These forms
        require named features. A one-dimensional sequence provides one
        aggregate ID per feature in table order. Features assigned the same
        aggregate ID are summed. By default, no aggregation is performed.
    p_adjust : str, optional
        Multiple testing correction method. Default is "holm".
    pseudo : float, optional
        Pseudo-count to add to zeros. Default is 0.
    s0_perc : float, optional
        SAM-like fudge factor percentile. Default is 0.05.
    alpha : float, optional
        Significance level. Default is 0.05.
    max_iter : int, optional
        Maximum EM iterations. Default is 100.
    tol : float, optional
        EM convergence tolerance. Default is 1e-5.

    Returns
    -------
    :class:`ANCOMBCResult`
        Result object with primary results and post-hoc analysis methods.

    See Also
    --------
    ancombc : ANCOM-BC without multi-group tests or sensitivity analysis.
    struc_zero : Standalone structural zero detection.

    References
    ----------
    .. [1] Lin, H. and Peddada, S.D., 2024. Multigroup analysis of
       compositions of microbiomes with covariate adjustments and repeated
       measures. Nature Methods, 21(1), 83–91.

    Examples
    --------
    >>> from skbio.stats.composition import ancombc2
    >>> import pandas as pd

    >>> samples = [f'S{i}' for i in range(1, 16)]
    >>> features = [f'F{i}' for i in range(1, 9)]
    >>> data = [[ 2,  0,  7,  0,  0,  2,  3,  2],
    ...         [ 0,  2,  0,  1,  1,  2,  0,  0],
    ...         [ 3,  1,  0,  9,  0,  1,  1,  0],
    ...         [ 2,  0,  1,  8,  1,  0, 11, 46],
    ...         [ 1,  0,  1,  1,  0,  0,  2,  2],
    ...         [ 0,  3, 22,  1,  1,  1,  3,  0],
    ...         [ 1,  7, 16,  1,  0,  0,  2,  2],
    ...         [ 0,  5,  6,  1,  2,  1,  0,  1],
    ...         [ 1,  7,  0,  2,  1,  0,  3,  2],
    ...         [ 0,  6,  4,  2,  0,  2,  2,  1],
    ...         [ 3, 13,  7,  0,  0,  0,  3,  9],
    ...         [ 1,  8,  5,  1,  0,  0,  0,  0],
    ...         [ 0,  5, 14,  1,  0,  1,  0,  1],
    ...         [ 5, 26,  3,  2,  0,  3,  1,  3],
    ...         [ 0, 18,  7,  0,  0,  3,  1,  0]]
    >>> table = pd.DataFrame(data, index=samples, columns=features)
    >>> status = ['mild'] * 5 + ['moderate'] * 5 + ['severe'] * 5
    >>> age = [39, 19, 20, 31, 15, 37, 27, 47, 26, 23, 39, 48, 46, 33, 36]
    >>> metadata = pd.DataFrame({'status': status, 'age': age}, index=samples)

    >>> res = ancombc2(table + 1, metadata, formula='status + age')
    >>> res_main = res.res
    >>> res_main.round(3)
                                  Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F1        Intercept              0.342  0.723  0.473   0.645   1.000   False
              status[T.moderate]    -0.337  0.517 -0.652   0.528   1.000   False
              status[T.severe]       0.278  0.673  0.413   0.688   1.000   False
              age                    0.001  0.024  0.054   0.958   1.000   False
    F2        Intercept             -0.541  0.796 -0.680   0.511   1.000   False
              status[T.moderate]     1.906  0.527  3.618   0.004   0.032    True
              status[T.severe]       2.935  0.639  4.590   0.001   0.006    True
              age                   -0.022  0.025 -0.877   0.399   1.000   False
    F3        Intercept             -2.242  0.893 -2.509   0.029   0.232   False
              status[T.moderate]     1.010  0.788  1.282   0.226   1.000   False
              status[T.severe]       0.717  0.803  0.893   0.391   1.000   False
              age                    0.063  0.032  1.955   0.077   0.612   False
    F4        Intercept              0.580  0.824  0.703   0.497   1.000   False
              status[T.moderate]    -0.046  0.639 -0.071   0.945   1.000   False
              status[T.severe]      -0.244  0.819 -0.298   0.771   1.000   False
              age                   -0.004  0.029 -0.126   0.902   1.000   False
    F5        Intercept             -0.502  0.757 -0.663   0.521   1.000   False
              status[T.moderate]     0.274  0.552  0.496   0.630   1.000   False
              status[T.severe]      -0.323  0.651 -0.497   0.629   1.000   False
              age                    0.027  0.025  1.069   0.308   1.000   False
    F6        Intercept             -0.045  0.791 -0.057   0.956   1.000   False
              status[T.moderate]     0.114  0.662  0.172   0.866   1.000   False
              status[T.severe]       0.375  0.825  0.455   0.658   1.000   False
              age                    0.008  0.028  0.275   0.788   1.000   False
    F7        Intercept              0.291  0.816  0.357   0.728   1.000   False
              status[T.moderate]     0.082  0.559  0.146   0.887   1.000   False
              status[T.severe]      -0.264  0.711 -0.371   0.718   1.000   False
              age                    0.004  0.027  0.150   0.883   1.000   False
    F8        Intercept             -0.118  0.858 -0.138   0.893   1.000   False
              status[T.moderate]    -0.402  0.793 -0.507   0.622   1.000   False
              status[T.severe]      -0.299  0.984 -0.304   0.767   1.000   False
              age                    0.022  0.029  0.763   0.462   1.000   False

    """
    return _ancombc(
        table=table,
        metadata=metadata,
        reestimate=True,
        formula=formula,
        aggregator=aggregator,
        p_adjust=p_adjust,
        pseudo=pseudo,
        s0_perc=s0_perc,
        alpha=alpha,
        max_iter=max_iter,
        tol=tol,
    )


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
    dmat, group, beta_hat, vcov_hat, dof=None, p_adj_method="holm", alpha=0.05
):
    """ANCOM-BC2 global test using F-test or chi-square test.

    Parameters
    ----------
    dmat : patsy.DesignMatrix
        Design matrix generated from the fitted formula.
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
    covariates = dmat.design_info.column_names

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
    dmat, group, beta_hat, var_hat, vcov_hat, dof, fwer_ctrl_method="holm", alpha=0.05
):
    """ANCOM-BC2 pairwise directional test.

    For each pair of group levels, compute the difference in coefficients
    and its variance, then apply mdFDR correction.
    """
    covariates = dmat.design_info.column_names

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
        dmat=dmat,
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
    W, dof, fwer_ctrl_method, dmat, group, beta_hat, vcov_hat, alpha, dof_global=None
):
    """Mixed directional FDR correction for pairwise tests.

    1. Run global test to screen significant taxa (count R).
    2. Only consider R significant taxa for pairwise p-values.
    3. Adjust at level R * alpha / d.
    """
    n_tax, n_comp = W.shape

    # Step 1: Global test screening
    res_screen = _ancombc_global_F(
        dmat=dmat,
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


def _dunn_global(dmat, group, W, B, dof, p_adj_method="holm", alpha=0.05):
    """Dunnett's global test for mdFDR.

    Bootstrap-based: generate null W from t-distribution, take max |W|.
    """
    n_tax = W.shape[0]
    covariates = dmat.design_info.column_names
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
    dmat, group, beta_hat, var_hat, dof, B=100, fwer_ctrl_method="holm", alpha=0.05
):
    """ANCOM-BC2 Dunnett's type of test.

    Compare each group to the reference group with mdFDR correction.
    """
    covariates = dmat.design_info.column_names
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
        dmat=dmat,
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


def _mdfdr_dunnett(W, dof, fwer_ctrl_method, dmat, group, B, alpha):
    """mdFDR correction for Dunnett's test."""
    n_tax = W.shape[0]

    # Step 1: Global test screening via bootstrap
    res_screen = _dunn_global(
        dmat=dmat, group=group, W=W, B=B, dof=dof, p_adj_method="BH", alpha=alpha
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
    dmat,
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
    covariates = dmat.design_info.column_names
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
