# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""Batched REML fitting of random-intercept linear mixed-effects models.

Fits :math:`y = X\beta + u_g + \epsilon` with one random intercept per group,
for many response vectors that share a single design, in one set of array
operations. This is the model that
:class:`~statsmodels.regression.mixed_linear_model.MixedLM` fits by default, and
the estimates, standard errors and log-likelihood here agree with it to
floating-point precision at the same variance ratio.

Notes
-----
Writing :math:`\theta` for the ratio of the group variance to the residual
variance, the marginal covariance of group :math:`g` is
:math:`\sigma^2 (I + \theta 1 1^T)`, whose inverse and log-determinant are
available in closed form. Profiling out :math:`\beta` and :math:`\sigma^2`
leaves a one-dimensional maximization over :math:`\theta`, in which every term
is a contraction of per-group sums of the design with per-group sums of the
response. The design half is therefore computed once and reused for every
response vector, and the maximization runs for all of them at once.

"""

from typing import NamedTuple

import numpy as np

# Search bounds and resolution for the maximization, in log(theta).
_LOG_LO = -20.0
_LOG_HI = 20.0
_N_GRID = 65
_N_REFINE = 45
_STEP = (_LOG_HI - _LOG_LO) / (_N_GRID - 1)
_LOG_GRID = np.linspace(_LOG_LO, _LOG_HI, _N_GRID)

# Golden-section constant of the batched line search.
_INV_PHI = 0.5 * (np.sqrt(5.0) - 1.0)

# Largest number of cells held in one temporary array. Features are fitted in
# blocks sized against this, bounding both the (n_groups, n_feats) reductions
# and the (n_feats, n_covars, n_covars) stacks.
_BLOCK_CELLS = 1 << 21


class _Summary(NamedTuple):
    """Per-feature sufficient statistics of a response matrix."""

    xty: np.ndarray  # X'y, shape (n_covars, n_feats)
    gsum: np.ndarray  # per-group sums of y, shape (n_groups, n_feats)
    gsum2: np.ndarray  # their squares, shape (n_groups, n_feats)
    yty: np.ndarray  # y'y, shape (n_feats,)


class _GridPoint(NamedTuple):
    """Design-only quantities at one fixed point of the search grid."""

    scaled_sums: np.ndarray  # (X_g'1 scaled by the group weight), (k, g)
    weight: np.ndarray  # theta / (1 + n_g theta), (g,)
    xvx_inv: np.ndarray  # (X'V^-1 X)^-1, (k, k)
    const: float  # the part of the log-likelihood free of the response
    usable: bool  # False if X'V^-1 X is singular at this theta


def _randint_applicable(
    re_formula, vc_formula, model_kwargs, fit_kwargs, fit_method, fit_converge
):
    """Report whether the batched fit can serve this model configuration.

    The batched fit covers only the default model, a single random intercept
    per group fitted by REML. Anything that selects a different random effects
    structure, a different optimizer, or a convergence filter that only an
    iterative optimizer can supply must be fitted one feature at a time.

    """
    return (
        re_formula is None
        and vc_formula is None
        and not model_kwargs
        and not fit_kwargs
        and fit_method is None
        and not fit_converge
    )


class _RandIntDesign:
    """Design-only quantities of a random-intercept mixed model.

    Holds every part of the REML profile likelihood that depends on the fixed
    effects design and the grouping but not on the response, so that it is
    computed once and reused for all features in all replicates.

    Parameters
    ----------
    exog : ndarray of shape (n_samples, n_covars)
        Fixed effects design matrix.
    groups : ndarray of shape (n_samples,)
        Group index of each sample, taking values in ``range(n_groups)``.
    n_groups : int
        Number of groups.

    """

    def __init__(self, exog, groups, n_groups):
        self.exog = np.ascontiguousarray(exog, dtype=np.float64)
        n_obs, self.n_covars = self.exog.shape
        self.n_groups = n_groups
        self.dof = float(n_obs - self.n_covars)
        self.sizes = np.bincount(groups, minlength=n_groups).astype(np.float64)

        # Rows ordered by group, so per-group sums are a `reduceat` reduction.
        self.order = np.argsort(groups, kind="stable")
        self.starts = np.searchsorted(groups[self.order], np.arange(n_groups))

        self.sums = np.add.reduceat(self.exog[self.order], self.starts, axis=0)
        k = self.n_covars
        self.outer = np.ascontiguousarray(
            (self.sums[:, :, None] * self.sums[:, None, :]).reshape(n_groups, k * k)
        )
        self.xtx_flat = (self.exog.T @ self.exog).ravel()

        # Repeated group sizes collapse the log-determinant of V to one
        # logarithm per distinct size rather than one per group.
        uniq, counts = np.unique(self.sizes, return_counts=True)
        self.uniq_sizes = uniq
        self.uniq_counts = counts.astype(np.float64)

        # Additive constant of the REML log-likelihood, free of theta. A design
        # with no residual degrees of freedom cannot be fitted at all; leave the
        # constant at zero and let the linear algebra raise.
        if self.dof > 0:
            self.llf_const = (
                0.5 * self.dof * (np.log(self.dof) - np.log(2 * np.pi) - 1.0)
            )
        else:
            self.llf_const = 0.0

        self.grid = [self._grid_point(np.exp(x)) for x in _LOG_GRID]

    def _grid_point(self, theta):
        """Precompute the design half of the objective at one fixed theta."""
        k = self.n_covars
        weight = theta / (1.0 + self.sizes * theta)
        xvx = (self.xtx_flat - weight @ self.outer).reshape(k, k)
        sign, logdet_x = np.linalg.slogdet(xvx)
        if sign <= 0:
            return _GridPoint(None, weight, None, -np.inf, False)
        logdet_v = self.uniq_counts @ np.log1p(self.uniq_sizes * theta)
        return _GridPoint(
            scaled_sums=np.ascontiguousarray((self.sums * weight[:, None]).T),
            weight=weight,
            xvx_inv=np.linalg.inv(xvx),
            const=-0.5 * logdet_v - 0.5 * logdet_x + self.llf_const,
            usable=True,
        )

    def summarize(self, resp):
        """Reduce a response matrix to per-feature sufficient statistics."""
        gsum = np.add.reduceat(resp[self.order], self.starts, axis=0)
        return _Summary(
            xty=self.exog.T @ resp,
            gsum=gsum,
            gsum2=gsum * gsum,
            yty=np.einsum("ij,ij->j", resp, resp),
        )


def _randint_llf(des, logdet_v, logdet_x, qform):
    """Assemble the REML profile log-likelihood from its three terms."""
    # A degenerate feature can drive `qform` to zero; the resulting -inf makes
    # the search reject that theta, and the caller screens the final fit.
    with np.errstate(divide="ignore", invalid="ignore"):
        return (
            -0.5 * logdet_v
            - 0.5 * des.dof * np.log(qform)
            - 0.5 * logdet_x
            + des.llf_const
        )


def _randint_profile(des, theta, stats):
    """Evaluate the objective and its by-products at per-feature ``theta``."""
    k, n_feats = des.n_covars, theta.shape[0]
    dinv = 1.0 / (1.0 + des.sizes[:, None] * theta[None, :])
    weight = theta[None, :] * dinv
    xvx = (des.xtx_flat[None, :] - np.ascontiguousarray(weight.T) @ des.outer).reshape(
        n_feats, k, k
    )

    # One feature whose theta makes X'V^-1 X singular must not abort the whole
    # batch: it is solved against the identity and then scored as impossible.
    sign, logdet_x = np.linalg.slogdet(xvx)
    bad = sign <= 0
    if bad.any():
        xvx = np.where(bad[:, None, None], np.eye(k), xvx)

    xvy = des.sums.T @ (weight * stats.gsum)
    np.subtract(stats.xty, xvy, out=xvy)
    yvy = stats.yty - np.einsum("gp,gp->p", weight, stats.gsum2)
    beta = np.linalg.solve(xvx, xvy.T[:, :, None])[:, :, 0].T
    qform = yvy - np.einsum("kp,kp->p", xvy, beta)
    logdet_v = des.uniq_counts @ np.log1p(np.outer(des.uniq_sizes, theta))
    llf = _randint_llf(des, logdet_v, logdet_x, qform)
    if bad.any():
        beta[:, bad] = np.nan
        llf = np.where(bad, -np.inf, llf)
    return llf, dinv, xvx, beta, qform


def _randint_scan(des, stats, n_feats):
    """Bracket the maximizer of each feature on the shared log-theta grid."""
    best = np.zeros(n_feats, dtype=int)
    best_llf = np.full(n_feats, -np.inf)
    for pos, point in enumerate(des.grid):
        if not point.usable:
            continue
        xvy = stats.xty - point.scaled_sums @ stats.gsum
        yvy = stats.yty - point.weight @ stats.gsum2
        qform = yvy - np.einsum("kp,kp->p", xvy, point.xvx_inv @ xvy)
        with np.errstate(divide="ignore", invalid="ignore"):
            cur = point.const - 0.5 * des.dof * np.log(qform)
        better = cur > best_llf
        best_llf[better] = cur[better]
        best[better] = pos
    return best


def _randint_fit(des, resp):
    """Fit one random-intercept model per column of ``resp``.

    Parameters
    ----------
    des : _RandIntDesign
        Precomputed design quantities.
    resp : ndarray of shape (n_samples, n_feats)
        Response vectors, one per feature.

    Returns
    -------
    beta : ndarray of shape (n_covars, n_feats)
        Fixed effects estimates.
    bse : ndarray of shape (n_covars, n_feats)
        Standard errors of the fixed effects.
    ok : ndarray of shape (n_feats,)
        Features whose fit produced usable statistics.
    theta : ndarray of shape (n_feats,)
        Fitted ratio of the group variance to the residual variance.
    llf : ndarray of shape (n_feats,)
        REML log-likelihood at the fit.

    """
    n_feats = resp.shape[1]
    block = max(1, _BLOCK_CELLS // max(des.n_groups, des.n_covars**2, 1))
    if n_feats <= block:
        return _randint_fit_block(des, resp)
    parts = [
        _randint_fit_block(des, resp[:, i : i + block])
        for i in range(0, n_feats, block)
    ]
    return tuple(np.concatenate(x, axis=-1) for x in zip(*parts))


def _randint_fit_block(des, resp):
    """Fit a block of features small enough to hold in one set of temporaries."""
    stats = des.summarize(resp)

    # Bracket the maximizer on a shared log-theta grid before refining it. A
    # flat or boundary-seeking likelihood therefore cannot strand the search at
    # an arbitrary starting point, which is how the iterative optimizers in
    # `MixedLM.fit` fail on these data.
    best = _randint_scan(des, stats, resp.shape[1])
    lo = _LOG_LO + (best - 1) * _STEP
    hi = _LOG_LO + (best + 1) * _STEP

    # Batched golden-section search on log theta.
    left = hi - _INV_PHI * (hi - lo)
    right = lo + _INV_PHI * (hi - lo)
    f_left = -_randint_profile(des, np.exp(left), stats)[0]
    f_right = -_randint_profile(des, np.exp(right), stats)[0]
    for _ in range(_N_REFINE):
        take_left = f_left < f_right
        hi = np.where(take_left, right, hi)
        lo = np.where(take_left, lo, left)
        span = _INV_PHI * (hi - lo)
        new_left = np.where(take_left, hi - span, right)
        new_right = np.where(take_left, left, lo + span)
        probe = np.where(take_left, new_left, new_right)
        f_probe = -_randint_profile(des, np.exp(probe), stats)[0]
        f_left, f_right = (
            np.where(take_left, f_probe, f_right),
            np.where(take_left, f_left, f_probe),
        )
        left, right = new_left, new_right

    return _randint_stats(des, np.exp(0.5 * (lo + hi)), stats)


def _randint_stats(des, theta, stats):
    """Fixed effects and their standard errors at the fitted ``theta``."""
    llf, dinv, xvx, beta, qform = _randint_profile(des, theta, stats)
    dof, sums, sizes = des.dof, des.sums, des.sizes
    k, n_feats = des.n_covars, theta.shape[0]

    with np.errstate(divide="ignore", invalid="ignore"):
        rsum = dinv * (stats.gsum - sums @ beta)
        ssum = sizes[:, None] * dinv
        dinv2 = dinv * dinv

        cross = sums.T @ (rsum * dinv)
        b_term = np.einsum("gp,gp->p", rsum, rsum)
        d_term = 2.0 * np.einsum("gp,gp,gp->p", ssum, rsum, rsum)
        hess_re = 0.5 * np.einsum("gp,gp->p", ssum, ssum) - 0.5 * dof * (
            d_term / qform - b_term**2 / qform**2
        )
        pmat = (np.ascontiguousarray(dinv2.T) @ des.outer).reshape(n_feats, k, k)
        fmat = 2.0 * (np.ascontiguousarray((ssum * dinv2).T) @ des.outer).reshape(
            n_feats, k, k
        )
        xvx_inv = np.linalg.inv(xvx)
        prod = xvx_inv @ pmat
        hess_re = hess_re + 0.5 * (
            np.einsum("pij,pji->p", prod, prod) - np.einsum("pij,pji->p", xvx_inv, fmat)
        )

        # Fixed effects block of the inverse observed information, obtained as
        # the Schur complement against the single variance parameter. This is
        # what `MixedLM` reports; it differs from scale * (X'V^-1 X)^-1 by
        # about 1%.
        prec = dof / qform
        schur = prec[:, None, None] * xvx + np.einsum(
            "kp,lp,p->pkl", cross, cross, prec**2 / hess_re
        )
        bse = np.sqrt(np.einsum("pkk->pk", np.linalg.inv(schur))).T

    # An indefinite Schur complement means the observed information is not
    # positive definite, which happens when theta is estimated at the boundary.
    # `MixedLM` returns NaN standard errors in the same situation. A zero
    # standard error would turn into an infinite test statistic, so screen it
    # out here rather than report a p-value of exactly zero.
    ok = (
        np.isfinite(beta).all(axis=0)
        & np.isfinite(bse).all(axis=0)
        & (bse > 0).all(axis=0)
    )
    return beta, bse, ok, theta, llf
