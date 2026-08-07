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
from dataclasses import dataclass, field, fields
from typing import Optional
from scipy.stats import norm, chi2, t
from scipy.optimize import minimize
from patsy import dmatrix

from skbio.table._tabular import _ingest_table
from skbio.stats.composition import clr
from ._base import _check_composition
from ._utils import _check_metadata, _check_p_adjust, _type_cast_to_float

# Multi-group test helpers (global, pairwise, Dunnett's, trend) are in
# _ancombc2_tests.py for modularity.
from ._ancombc2_tests import (
    _var_diff,
    _combn_fun,
    _combn_fun2,
    _ancombc_global_F,
    _ancombc_pair,
    _mdfdr_pairwise,
    _dunn_global,
    _ancombc_dunn,
    _mdfdr_dunnett,
    _safe_inverse_spd,
    _constrain_est,
    _l_infty,
    _ancombc_trend,
)

# Convergence tolerance for Nelder-Mead optimization in _estimate_bias_em selected
# according to the benchmark. SciPy's default is 1e-4, which would cause different
# results to R version.
NM_TOL = 1e-8


@dataclass
class ANCOMBCResult:
    """Results for ANCOM-BC and ANCOM-BC2 analyses.

    This class contains the primary differential abundance results. Post-hoc
    analyses (structural zeros, global test, multi-group comparisons,
    sensitivity analysis) are available as methods that compute on-demand
    using stored intermediate data
    Returned by both :func:`ancombc` and :func:`ancombc2`.

    Attributes
    ----------
    res : pd.DataFrame
        Primary results with (FeatureID, Covariate) multi-index.
        For ANCOM-BC (reestimate=False): columns Log2(FC), SE, W, pvalue,
        qvalue, Signif.
        For ANCOM-BC2 (reestimate=True): columns lfc, SE, W, pvalue, qvalue,
        Signif.
    delta_em : ndarray
        EM-estimated biases per covariate.
    delta_wls : ndarray
        WLS-estimated biases per covariate.
    samp_frac : ndarray
        Estimated sampling fractions per sample.
    feature_table : pd.DataFrame or None
        Preprocessed feature table. None for ANCOM-BC.
    reestimate : bool
        Whether ANCOM-BC2 (True) or ANCOM-BC (False) was used.

    Methods
    -------
    structural_zeros()
        Detect structural zeros (features systematically absent per group).
    global_test()
        Global test for differential abundance across >= 3 groups.
    dunnett_test()
        Dunnett's test: each group vs. reference, with mdFDR correction.
    pairwise_test()
        Pairwise directional test between all group pairs, with mdFDR.
    trend_test()
        Trend test for ordered patterns in group effects.
    sensitivity_analysis()
        Pseudo-count sensitivity analysis for robustness assessment.

    See Also
    --------
    ancombc : ANCOM-BC function (reestimate=False).
    ancombc2 : ANCOM-BC2 function (reestimate=True).
    """

    # Primary results
    res: pd.DataFrame

    # Diagnostic values
    delta_em: np.ndarray = field(default=None, repr=False)
    delta_wls: np.ndarray = field(default=None, repr=False)
    samp_frac: np.ndarray = field(default=None, repr=False)

    # Preprocessing outputs
    feature_table: Optional[pd.DataFrame] = field(default=None, repr=False)

    # Determine whether ANCOM-BC2 (True) or ANCOM-BC (False) is used
    reestimate: bool = field(default=False, repr=False)

    # Private intermediate data
    _table: object = field(default=None, repr=False)
    _metadata: object = field(default=None, repr=False)
    _group: Optional[str] = field(default=None, repr=False)
    _grouping: Optional[str] = field(default=None, repr=False)
    _dmat: object = field(default=None, repr=False)
    _x_df: object = field(default=None, repr=False)
    _beta_hat: object = field(default=None, repr=False)
    _var_hat: object = field(default=None, repr=False)
    _vcov_hat: object = field(default=None, repr=False)
    _dof: object = field(default=None, repr=False)
    _tax_name: object = field(default=None, repr=False)
    _fix_eff: object = field(default=None, repr=False)
    _alpha: float = field(default=0.05, repr=False)
    _p_adjust: str = field(default="holm", repr=False)
    _neg_lb: bool = field(default=False, repr=False)
    _mdfdr_fwer_ctrl_method: str = field(default="holm", repr=False)
    _mdfdr_B: int = field(default=100, repr=False)
    _trend_contrast: object = field(default=None, repr=False)
    _trend_node: object = field(default=None, repr=False)
    _trend_B: int = field(default=100, repr=False)
    _s0_perc: float = field(default=0.05, repr=False)
    _max_iter: int = field(default=100, repr=False)
    _tol: float = field(default=1e-5, repr=False)
    _pseudo: float = field(default=0, repr=False)
    _O1: object = field(default=None, repr=False)
    _O2: object = field(default=None, repr=False)
    _formula: Optional[str] = field(default=None, repr=False)

    # private caches
    _structural_zeros: object = field(default=None, repr=False)
    _global_test: object = field(default=None, repr=False)
    _dunnett_test: object = field(default=None, repr=False)
    _pairwise_test: object = field(default=None, repr=False)
    _trend_test: object = field(default=None, repr=False)
    _sensitivity: object = field(default=None, repr=False)

    def __getitem__(self, key):
        return getattr(self, key)

    def keys(self):
        """Return a list of attribute names that are not private and are not None."""
        return [
            f.name
            for f in fields(self)
            if not f.name.startswith("_") and getattr(self, f.name) is not None
        ]

    # formated output
    def __repr__(self):
        n_taxa = len(self.res.index.get_level_values("FeatureID").unique())
        n_cov = len(self.res.index.get_level_values("Covariate").unique())
        n_signif = int(self.res["Signif"].sum())
        label = "ANCOM-BC2" if self.reestimate else "ANCOM-BC"
        return (
            f"ANCOMBCResult(method={label!r}, "
            f"n_taxa={n_taxa}, n_covariates={n_cov}, "
            f"n_signif={n_signif})"
        )

    def structural_zeros(self) -> pd.DataFrame:
        """Detect structural zeros in the original table.

        Structural zeros are features systematically absent from certain
        sample groups. This method tests whether the proportion of observed
        zeros is close to zero, suggesting absence of a feature in a group.

        Returns
        -------
        pd.DataFrame of bool
            DataFrame of shape (n_features, n_groups) indicating whether
            each feature is a structural zero in each group.

        Raises
        ------
        ValueError
            If no group variable was specified (``group`` was None when
            :func:`ancombc2` was called).

        Notes
        -----
        Results are cached: repeated calls return the same object without
        recomputation.
        """
        if self._structural_zeros is not None:
            return self._structural_zeros
        if self._group is None:
            raise ValueError(
                "structural_zeros requires a group variable; "
                "pass group=... to ancombc2()."
            )
        _struc_zero_fn = globals()["struc_zero"]
        result = _struc_zero_fn(self._table, self._metadata, self._group, self._neg_lb)
        self._structural_zeros = result
        return result

    def global_test(self) -> pd.DataFrame:
        """Perform global test for differential abundance across groups.

        The global test identifies features that are differentially abundant
        between at least two groups across three or more groups.

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
        Results are cached: repeated calls return the same object without
        recomputation.
        """
        if self._global_test is not None:
            return self._global_test

        if not self.reestimate:
            if self._grouping is None:
                raise ValueError(
                    "global_test requires a grouping variable; "
                    "pass grouping=... to ancombc()."
                )
            W_g, pval, qval, reject = _global_test(
                self._dmat,
                self._grouping,
                self._beta_hat,
                self._vcov_hat,
                self._alpha,
                self._p_adjust,
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
            if self._group is None:
                raise ValueError(
                    "global_test requires a group variable; "
                    "pass group=... to ancombc2()."
                )
            raw = _ancombc_global_F(
                x=self._x_df,
                group=self._group,
                beta_hat=self._beta_hat,
                vcov_hat=self._vcov_hat,
                dof=self._dof,
                p_adj_method=self._p_adjust,
                alpha=self._alpha,
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

        self._global_test = result
        return result

    def dunnett_test(self) -> pd.DataFrame:
        """Perform Dunnett's test (each group vs. reference) with mdFDR.

        Returns
        -------
        pd.DataFrame
            DataFrame with (FeatureID, Comparison) multi-index and columns:
            lfc, SE, W, pvalue, qvalue, Signif.

        Raises
        ------
        ValueError
            If this is an ANCOM-BC result (reestimate=False), or if no
            group variable was specified.
        """
        if self._dunnett_test is not None:
            return self._dunnett_test
        if not self.reestimate:
            raise ValueError("dunnett_test() is only available for ANCOM-BC2 results.")
        if self._group is None:
            raise ValueError(
                "dunnett_test requires a group variable; pass group=... to ancombc2()."
            )

        raw = _ancombc_dunn(
            x=self._x_df,
            group=self._group,
            beta_hat=self._beta_hat,
            var_hat=self._var_hat,
            dof=self._dof,
            B=self._mdfdr_B,
            fwer_ctrl_method=self._mdfdr_fwer_ctrl_method,
            alpha=self._alpha,
        )
        comp_names = raw["comp_names"]
        n_comp = len(comp_names)
        n_tax = len(self._tax_name)
        result = pd.DataFrame(
            {
                "FeatureID": [x for x in self._tax_name for _ in range(n_comp)],
                "Comparison": comp_names * n_tax,
                "lfc": raw["beta"].ravel(),
                "SE": raw["se"].ravel(),
                "W": raw["W"].ravel(),
                "pvalue": raw["p_val"].ravel(),
                "qvalue": raw["q_val"].ravel(),
                "Signif": raw["diff_abn"].ravel(),
            }
        )
        result.set_index(["FeatureID", "Comparison"], inplace=True)
        self._dunnett_test = result
        return result

    def pairwise_test(self) -> pd.DataFrame:
        """Perform pairwise directional test between all group pairs.

        Uses mixed directional FDR (mdFDR) correction via bootstrap.

        Returns
        -------
        pd.DataFrame
            DataFrame with (FeatureID, Comparison) multi-index and columns:
            lfc, SE, W, pvalue, qvalue, Signif.

        Raises
        ------
        ValueError
            If this is an ANCOM-BC result (reestimate=False), or if no
            group variable was specified.
        """
        if self._pairwise_test is not None:
            return self._pairwise_test
        if not self.reestimate:
            raise ValueError("pairwise_test() is only available for ANCOM-BC2 results.")
        if self._group is None:
            raise ValueError(
                "pairwise_test requires a group variable; pass group=... to ancombc2()."
            )

        raw = _ancombc_pair(
            x=self._x_df,
            group=self._group,
            beta_hat=self._beta_hat,
            var_hat=self._var_hat,
            vcov_hat=self._vcov_hat,
            dof=self._dof,
            fwer_ctrl_method=self._mdfdr_fwer_ctrl_method,
            alpha=self._alpha,
        )
        comp_names = raw["comp_names"]
        n_comp = len(comp_names)
        n_tax = len(self._tax_name)
        result = pd.DataFrame(
            {
                "FeatureID": [x for x in self._tax_name for _ in range(n_comp)],
                "Comparison": comp_names * n_tax,
                "lfc": raw["beta"].ravel(),
                "SE": raw["se"].ravel(),
                "W": raw["W"].ravel(),
                "pvalue": raw["p_val"].ravel(),
                "qvalue": raw["q_val"].ravel(),
                "Signif": raw["diff_abn"].ravel(),
            }
        )
        result.set_index(["FeatureID", "Comparison"], inplace=True)
        self._pairwise_test = result
        return result

    def trend_test(self) -> pd.DataFrame:
        """Perform trend test for ordered patterns in group effects.

        Uses constrained optimization to test monotone increasing/decreasing
        patterns in group-level effects.

        Returns
        -------
        pd.DataFrame
            DataFrame indexed by FeatureID with columns: W, pvalue, qvalue,
            Signif.

        Raises
        ------
        ValueError
            If this is an ANCOM-BC result (reestimate=False), or if no
            group variable was specified.
        """
        if self._trend_test is not None:
            return self._trend_test
        if not self.reestimate:
            raise ValueError("trend_test() is only available for ANCOM-BC2 results.")
        if self._group is None:
            raise ValueError(
                "trend_test requires a group variable; pass group=... to ancombc2()."
            )

        raw = _ancombc_trend(
            x=self._x_df,
            group=self._group,
            beta_hat=self._beta_hat,
            var_hat=self._var_hat,
            vcov_hat=self._vcov_hat,
            p_adj_method=self._p_adjust,
            alpha=self._alpha,
            trend_contrast=self._trend_contrast,
            trend_node=self._trend_node,
            trend_B=self._trend_B,
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
        self._trend_test = result
        return result

    def sensitivity_analysis(self) -> dict:
        """Perform pseudo-count sensitivity analysis.

        Re-runs the core estimation with pseudo-count values [0.1, 0.5, 1]
        and compares q-values to assess robustness of significance calls to
        the choice of pseudo-count.

        If multi-group tests (global, pairwise, Dunnett, trend) have already
        been computed via their respective methods, sensitivity is also
        assessed for those tests by re-running them with each pseudo-count.

        Returns
        -------
        dict
            Dictionary with keys:
            - ``passed_ss_prim``: bool ndarray, primary results robustness
            - ``passed_ss_global``: bool ndarray or None
            - ``passed_ss_pair``: bool ndarray or None
            - ``passed_ss_dunn``: bool ndarray or None
            - ``passed_ss_trend``: bool ndarray or None
            - ``ss_tab_prim``: float ndarray, proportion of pseudo-counts
              where q > alpha

        Raises
        ------
        ValueError
            If this is an ANCOM-BC result (reestimate=False).
        """
        if self._sensitivity is not None:
            return self._sensitivity
        if not self.reestimate:
            raise ValueError(
                "sensitivity_analysis() is only available for ANCOM-BC2 results."
            )

        n_tax = len(self._tax_name)
        n_cov = len(self._fix_eff)
        alpha = self._alpha

        # Determine which multi-group tests to include
        has_global = self._global_test is not None
        has_pair = self._pairwise_test is not None
        has_dunn = self._dunnett_test is not None
        has_trend = self._trend_test is not None
        if has_trend:
            has_global = True  # Trend requires global test

        # Original q_hat (n_tax, n_cov)
        q_hat_orig = self.res["qvalue"].values.reshape(n_tax, n_cov)

        # Re-run core estimation with each pseudo-count
        pseudo_list = [0.1, 0.5, 1]
        ss_list = []
        for pseudo_count in pseudo_list:
            res_pseudo = _ancombc2_estimate(
                data=self._O1,
                aggregate_data=self._O2,
                metadata=self._metadata,
                fix_formula=self._formula,
                p_adj_method=self._p_adjust,
                pseudo=pseudo_count,
                s0_perc=self._s0_perc,
                group=self._group,
                alpha=alpha,
                verbose=False,
                max_iter=self._max_iter,
                tol=self._tol,
            )
            # If multi-group tests are requested, run them for this pseudo
            if has_global and self._group is not None:
                res_pseudo["res_global"] = _ancombc_global_F(
                    x=res_pseudo["x_df"],
                    group=self._group,
                    beta_hat=res_pseudo["beta_hat"],
                    vcov_hat=res_pseudo["vcov_hat"],
                    dof=res_pseudo["dof"],
                    p_adj_method=self._p_adjust,
                    alpha=alpha,
                )
            if has_pair and self._group is not None:
                res_pseudo["res_pair"] = _ancombc_pair(
                    x=res_pseudo["x_df"],
                    group=self._group,
                    beta_hat=res_pseudo["beta_hat"],
                    var_hat=res_pseudo["var_hat"],
                    vcov_hat=res_pseudo["vcov_hat"],
                    dof=res_pseudo["dof"],
                    fwer_ctrl_method=self._mdfdr_fwer_ctrl_method,
                    alpha=alpha,
                )
            if has_dunn and self._group is not None:
                res_pseudo["res_dunn"] = _ancombc_dunn(
                    x=res_pseudo["x_df"],
                    group=self._group,
                    beta_hat=res_pseudo["beta_hat"],
                    var_hat=res_pseudo["var_hat"],
                    dof=res_pseudo["dof"],
                    B=self._mdfdr_B,
                    fwer_ctrl_method=self._mdfdr_fwer_ctrl_method,
                    alpha=alpha,
                )
            if has_trend and self._group is not None:
                res_pseudo["res_trend"] = _ancombc_trend(
                    x=res_pseudo["x_df"],
                    group=self._group,
                    beta_hat=res_pseudo["beta_hat"],
                    var_hat=res_pseudo["var_hat"],
                    vcov_hat=res_pseudo["vcov_hat"],
                    p_adj_method=self._p_adjust,
                    alpha=alpha,
                    trend_contrast=self._trend_contrast,
                    trend_node=self._trend_node,
                    trend_B=self._trend_B,
                )
            ss_list.append(res_pseudo)

        # Combine original and pseudo-count results
        all_q = [q_hat_orig] + [r["q_hat"] for r in ss_list]
        q_3d = np.stack(all_q, axis=-1)  # (n_tax, n_cov, n_pseudo)
        ss_tab_prim = np.mean(q_3d > alpha, axis=-1)
        passed_ss_prim = (ss_tab_prim == 0) | (ss_tab_prim == 1)

        # Global test sensitivity
        passed_ss_global = None
        if has_global and self._global_test is not None:
            q_globals = [self._global_test["qvalue"].values]
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
        if has_pair and self._pairwise_test is not None:
            q_pair_orig = self._pairwise_test["qvalue"].values.reshape(n_tax, -1)
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
        if has_dunn and self._dunnett_test is not None:
            q_dunn_orig = self._dunnett_test["qvalue"].values.reshape(n_tax, -1)
            q_dunns = [q_dunn_orig]
            for r in ss_list:
                if r.get("res_dunn") is not None:
                    q_dunns.append(r["res_dunn"]["q_val"])
                else:
                    q_dunns.append(np.ones_like(q_dunn_orig))
            q_dunn_3d = np.stack(q_dunns, axis=-1)
            ss_dunn = np.mean(q_dunn_3d > alpha, axis=-1)
            passed_ss_dunn = (ss_dunn == 0) | (ss_dunn == 1)

        # Trend sensitivity (uses global sensitivity)
        passed_ss_trend = None
        if has_trend and self._trend_test is not None:
            if passed_ss_global is not None:
                passed_ss_trend = passed_ss_global

        result = {
            "passed_ss_prim": passed_ss_prim,
            "passed_ss_global": passed_ss_global,
            "passed_ss_pair": passed_ss_pair,
            "passed_ss_dunn": passed_ss_dunn,
            "passed_ss_trend": passed_ss_trend,
            "ss_tab_prim": ss_tab_prim,
        }
        self._sensitivity = result
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
    aggregate_data,
    metadata,
    fix_formula,
    p_adj_method="holm",
    pseudo=0,
    s0_perc=0.05,
    group=None,
    alpha=0.05,
    max_iter=100,
    tol=1e-5,
):
    """ANCOM-BC2 core estimation"""
    # data preprocessing
    O1 = data + pseudo
    O2 = aggregate_data + pseudo

    n_tax = O2.shape[1]
    tax_name = O2.columns.tolist()

    # Build design matrix using patsy (consistent with existing ancombc)
    # TODO: use ancombc dmat
    x_df = dmatrix(f"~{fix_formula}", data=metadata, return_type="dataframe")
    fix_eff = x_df.columns.tolist()
    n_fix_eff = len(fix_eff)
    x_arr = x_df.values

    # CLR transformation
    # Per-taxon centering: y = log(x) - mean(log(x)) per taxon (row).
    O1_vals = O1.values.astype(float)
    y1 = clr(O1_vals, axis=0)

    # Initial estimates via SVD-based OLS (_estimate_params)
    var_hat1, beta1, _, _ = _estimate_params(data=y1, dmat=x_arr)
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
    resid_all = y1 - (x_arr @ beta1_corrected.T)  # (n_samp, n_feats)
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
    var_hat, beta_hat, _, beta_covmat2 = _estimate_params(data=y2_crt, dmat=x_arr)
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
        "feature_table": O2,
        "bias_correct_log_table": y_bias_crt,
        "samp_frac": theta_hat_arr,
        "delta_em": delta_em,
        "delta_wls": delta_wls,
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
        "x_df": x_df,
        "y2": y2,
        "res_global": None,
        "res_pair": None,
        "res_dunn": None,
        "res_trend": None,
    }


def _ancombc(
    table,
    metadata,
    reestimate=False,
    formula=None,
    grouping=None,
    aggregate_data=None,
    p_adjust="holm",
    pseudo=0,
    pseudo_sens=True,
    s0_perc=0.05,
    group=None,
    struc_zero=False,
    neg_lb=False,
    alpha=0.05,
    max_iter=100,
    tol=1e-5,
    mdfdr_fwer_ctrl_method="holm",
    mdfdr_B=100,
    trend_contrast=None,
    trend_node=None,
    trend_B=100,
    global_test=False,
    pairwise=False,
    dunnett=False,
    trend=False,
):
    """Private ANCOM-BC/BC2 core function.

    Dispatches on ``reestimate``:
    - ``reestimate=False``: ANCOM-BC path (log(matrix+1), no filtering,
      no sensitivity, Log2(FC) columns).
    - ``reestimate=True``: ANCOM-BC2 path (CLR centering, sensitivity,
      multi-group tests, lfc columns).
    """
    # Note: A pseudocount should have been added to the table by the user prior to
    # calling this function.
    matrix, samples, features = _ingest_table(table)
    _check_composition(np, matrix, nozero=True)
    n_feats = matrix.shape[1]
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

        # Estimate sampling fractions.
        samp_frac = _sample_fractions(matrix, dmat, beta_hat)

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
            delta_em=delta_em,
            delta_wls=delta_wls,
            samp_frac=samp_frac,
            reestimate=False,
            _table=pd.DataFrame(matrix, index=samples, columns=features),
            _metadata=metadata,
            _grouping=grouping,
            _dmat=dmat,
            _beta_hat=beta_hat,
            _vcov_hat=vcov_hat,
            _alpha=alpha,
            _p_adjust=p_adjust,
        )

    else:
        # Data preprocessing (no filtering. user pre-filters)
        O1 = pd.DataFrame(matrix, index=samples, columns=features)
        if aggregate_data is not None:
            agg_matrix, agg_samples, agg_features = _ingest_table(aggregate_data)
            _check_composition(np, agg_matrix, nozero=True)
            n_agg_feats = agg_matrix.shape[1]
            if agg_features is None:
                agg_features = np.arange(n_agg_feats)

        if aggregate_data is not None:
            O2 = pd.DataFrame(agg_matrix, index=agg_samples, columns=agg_features)
        else:
            O2 = O1

        # Step 3: Run ANCOMBC2 core estimation
        res_main = _ancombc2_estimate(
            data=O1,
            aggregate_data=O2,
            metadata=metadata,
            fix_formula=formula,
            p_adj_method=p_adjust,
            pseudo=pseudo,
            s0_perc=s0_perc,
            group=group,
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

        lfc = beta_hat

        res_dict = {
            "FeatureID": [x for x in tax_name for _ in range(n_cov)],
            "Covariate": list(fix_eff) * n_tax,
            "lfc": lfc.ravel(),  # TODO: should have the same results as ancombc
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
            delta_em=res_main["delta_em"],
            delta_wls=res_main["delta_wls"],
            samp_frac=res_main["samp_frac"],
            feature_table=res_main["feature_table"],
            reestimate=True,
            _table=pd.DataFrame(matrix, index=samples, columns=features),
            _metadata=metadata,
            _group=group,
            _x_df=res_main["x_df"],
            _beta_hat=res_main["beta_hat"],
            _var_hat=res_main["var_hat"],
            _vcov_hat=res_main["vcov_hat"],
            _dof=res_main["dof"],
            _tax_name=tax_name,
            _fix_eff=fix_eff,
            _alpha=alpha,
            _p_adjust=p_adjust,
            _neg_lb=neg_lb,
            _mdfdr_fwer_ctrl_method=mdfdr_fwer_ctrl_method,
            _mdfdr_B=mdfdr_B,
            _trend_contrast=trend_contrast,
            _trend_node=trend_node,
            _trend_B=trend_B,
            _s0_perc=s0_perc,
            _max_iter=max_iter,
            _tol=tol,
            _pseudo=pseudo,
            _O1=O1,
            _O2=O2,
            _formula=formula,
        )


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
    grouping : str, optional
        A metadata column name of interest for *global test*, which identifies features
        that are differentially abundant between at least two groups across three or
        more groups in that column. Must be one of the factors in ``formula``. Default
        is None, which skips global test.
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
    ANCOMBCResult
        Result object with primary results, global test (if grouping specified),
        and diagnostic values.

    See Also
    --------
    ancombc2 : ANCOM-BC2 with multi-group tests and sensitivity analysis.
    struc_zero : Standalone structural zero detection.

    References
    ----------
    .. [1] Lin, H. and Peddada, S.D., 2020. Analysis of compositions of microbiomes
       with bias correction. Nature Communications, 11(1), p.3514.

    """
    return _ancombc(
        table=table,
        metadata=metadata,
        reestimate=False,
        formula=formula,
        grouping=grouping,
        max_iter=max_iter,
        tol=tol,
        alpha=alpha,
        p_adjust=p_adjust,
    )


def ancombc2(
    table,
    metadata,
    formula,
    aggregate_data=None,
    p_adjust="holm",
    pseudo=0,
    pseudo_sens=True,
    s0_perc=0.05,
    group=None,
    struc_zero=False,
    neg_lb=False,
    alpha=0.05,
    max_iter=100,
    tol=1e-5,
    mdfdr_fwer_ctrl_method="holm",
    mdfdr_B=100,
    trend_contrast=None,
    trend_node=None,
    trend_B=100,
    global_test=False,
    pairwise=False,
    dunnett=False,
    trend=False,  # TODO: remove parameters, only use for post hoc anlaysis
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
    aggregate_data : table_like, optional
        Pre-aggregated data for regression. When provided, ``table`` is used
        for EM bias estimation and ``aggregate_data`` for final regression.
    p_adjust : str, optional
        Multiple testing correction method. Default is "holm".
    pseudo : float, optional
        Pseudo-count to add to zeros. Default is 0.
    pseudo_sens : bool, optional
        Whether to perform sensitivity analysis. Default is True.
    s0_perc : float, optional
        SAM-like fudge factor percentile. Default is 0.05.
    group : str, optional
        Group variable for multi-group tests. Default is None.
    struc_zero : bool, optional
        Whether to detect structural zeros. Default is False.
    neg_lb : bool, optional
        Whether to use negative lower bound for structural zeros. Default
        is False.
    alpha : float, optional
        Significance level. Default is 0.05.
    max_iter : int, optional
        Maximum EM iterations. Default is 100.
    tol : float, optional
        EM convergence tolerance. Default is 1e-5.
    mdfdr_fwer_ctrl_method : str, optional
        FWER control method for mdFDR. Default is "holm".
    mdfdr_B : int, optional
        Number of bootstrap iterations for mdFDR. Default is 100.
    trend_contrast : dict, optional
        Contrast matrix for trend test. Default is None.
    trend_node : list, optional
        Node list for trend test. Default is None.
    trend_B : int, optional
        Number of bootstrap iterations for trend test. Default is 100.
    global_test : bool, optional
        Whether to perform global test. Default is False.
    pairwise : bool, optional
        Whether to perform pairwise directional test. Default is False.
    dunnett : bool, optional
        Whether to perform Dunnett's test. Default is False.
    trend : bool, optional
        Whether to perform trend test. Default is False.

    Returns
    -------
    ANCOMBCResult
        Result object with primary results, post-hoc analyses (global,
        pairwise, Dunnett, trend, structural zeros), sensitivity analysis,
        and diagnostic values.

    See Also
    --------
    ancombc : ANCOM-BC without multi-group tests or sensitivity analysis.
    struc_zero : Standalone structural zero detection.

    References
    ----------
    .. [1] Lin, H. and Peddada, S.D., 2024. Multigroup analysis of
       compositions of microbiomes with covariate adjustments and repeated
       measures. Nature Methods, 21(1), 83–91.

    """
    return _ancombc(
        table=table,
        metadata=metadata,
        reestimate=True,
        formula=formula,
        aggregate_data=aggregate_data,
        p_adjust=p_adjust,
        pseudo=pseudo,
        pseudo_sens=pseudo_sens,
        s0_perc=s0_perc,
        group=group,
        struc_zero=struc_zero,
        neg_lb=neg_lb,
        alpha=alpha,
        max_iter=max_iter,
        tol=tol,
        mdfdr_fwer_ctrl_method=mdfdr_fwer_ctrl_method,
        mdfdr_B=mdfdr_B,
        trend_contrast=trend_contrast,
        trend_node=trend_node,
        trend_B=trend_B,
        global_test=global_test,
        pairwise=pairwise,
        dunnett=dunnett,
        trend=trend,
    )
