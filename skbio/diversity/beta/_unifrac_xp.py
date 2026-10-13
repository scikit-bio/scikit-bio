# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Array-API-generic UniFrac implementation.

This mirrors the per-node math of the fused Numba CUDA/HIP kernel in
``skbio.diversity.beta._unifrac_gpu`` (see ``_make_unifrac_kernel``), but is
written entirely in terms of ``xp`` operations (via
:func:`skbio.util._array.ingest_array`) rather than hardware-specific device
code. It therefore runs correctly on any array-API-compliant array --
plain NumPy on CPU, or CuPy/PyTorch on CUDA or ROCm -- with no GPU backend
(``numba-cuda``/``numba.hip``) required.

It serves as the fallback taken by ``engine='gpu'`` when no usable GPU
backend is detected, or when the fused kernel fails to build/run on the
running stack; see ``_unifrac_gpu.py``'s ``*_gpu_or_xp`` wrappers.

To avoid materializing an ``O(n_samples**2, n_nodes)`` tensor, the pairwise
loop runs over sample ``i`` in plain Python (cheap, ``n_samples``
iterations) and vectorizes over all ``j > i`` in one shot via broadcasting,
keeping peak memory at ``O(n_samples, n_nodes)`` per iteration. The
non-variance-adjusted unweighted and weighted paths instead reduce to a
pair of matrix products computed once up front (see
``_unweighted_unifrac_xp_core`` and ``_weighted_unifrac_xp_core``), since
their per-pair terms are separable across samples; variance adjustment
breaks that separability, so those paths keep the per-``i`` loop.
"""

import numpy as np

from skbio.util._array import ingest_array, _to_numpy


def _condensed_offset(i, n_samples):
    """Flat offset such that out[offset + k] is pair (i, i + 1 + k)."""
    return i * n_samples - (i * (i + 1)) // 2


def _variance_adjust_terms(xp, counts_i, counts_block, sample_totals, i):
    """Compute the reciprocal VAW term and a validity mask for row ``i``.

    Parameters
    ----------
    counts_i : array, shape (n_nodes,)
        Per-node counts for sample ``i``.
    counts_block : array, shape (n_remaining, n_nodes)
        Per-node counts for samples ``j > i``.
    sample_totals : array, shape (n_samples,)
        Per-sample sum of tip counts.
    i : int
        Row index.

    Returns
    -------
    inv_vaw : array, shape (n_remaining, n_nodes)
        ``1 / sqrt(mi * (m - mi))`` where valid, else ``0``.
    valid : array, shape (n_remaining, n_nodes), bool
        Where ``mi * (m - mi) > 0``.

    """
    m = sample_totals[i] + sample_totals[i + 1 :]  # (n_remaining,)
    mi = counts_i[None, :] + counts_block  # (n_remaining, n_nodes)
    # mi * (m - mi) is provably >= 0 for valid count data: a node's count is
    # a sum of a subset of its sample's tip counts, so mi (itself a sum of
    # two such non-negative node counts) can never exceed m (the matching
    # sum of the two samples' tip totals) -- no inner `where` is needed to
    # suppress a sqrt-of-negative warning here; `valid` below still gates
    # the (mi == m or mi == 0) boundary where the term is undefined.
    arg = mi * (m[:, None] - mi)
    valid = arg > 0
    # Guard the reciprocal itself (not `arg`'s sign, which needs no
    # guarding -- see above) against 0-division at invalid entries.
    safe_vaw = xp.where(valid, xp.sqrt(arg), xp.asarray(1.0, dtype=arg.dtype))
    inv_vaw = xp.where(valid, 1 / safe_vaw, xp.asarray(0.0, dtype=arg.dtype))
    return inv_vaw, valid


def _safe_divide(xp, numer, denom, valid):
    """``numer / denom`` where ``valid``, else ``0``; avoids 0-division warnings."""
    safe_denom = xp.where(valid, denom, xp.asarray(1.0, dtype=denom.dtype))
    return xp.where(valid, numer / safe_denom, xp.asarray(0.0, dtype=numer.dtype))


def _normalize(xp, numerator, denominator):
    """``numerator / denominator``, or ``0`` where the denominator vanishes.

    Matches the fused kernel's final ``0.0 if denominator == 0.0 else
    numerator / denominator``.
    """
    return _safe_divide(xp, numerator, denominator, denominator != 0)


def _weighted_unifrac_xp_core(
    xp,
    proportions,
    counts_by_node,
    sample_totals,
    branch_lengths,
    normalized,
    variance_adjust,
):
    n_samples = proportions.shape[0]
    n_pairs = n_samples * (n_samples - 1) // 2
    out = xp.zeros(n_pairs, dtype=proportions.dtype)

    # Non-VAW denominator is separable across samples: sum(bl * (p_i + p_j))
    # == B_i + B_j, where B = proportions @ branch_lengths. Computed once
    # here instead of being rebuilt every loop iteration. VAW rescales the
    # per-node terms by a (i, j, node)-dependent factor first, so this does
    # not hold once variance_adjust is True.
    if not variance_adjust:
        totals_by_sample = proportions @ branch_lengths  # (n_samples,)

    for i in range(n_samples - 1):
        p_u = proportions[i, :]
        p_v = proportions[i + 1 :, :]
        d = xp.abs(p_u[None, :] - p_v)

        if variance_adjust:
            s = p_u[None, :] + p_v
            inv_vaw, _ = _variance_adjust_terms(
                xp,
                counts_by_node[i, :],
                counts_by_node[i + 1 :, :],
                sample_totals,
                i,
            )
            # inv_vaw is already 0 where invalid, so multiplying by it (in
            # place of the old divide-through-_safe_divide) zeroes those
            # terms out directly, with no extra `where` needed.
            s = s * inv_vaw
            d = d * inv_vaw
            denominator = xp.sum(branch_lengths[None, :] * s, axis=-1)
        else:
            denominator = totals_by_sample[i] + totals_by_sample[i + 1 :]

        # Rejected by measurement: a matrix product against `d` here is
        # slower than this fused sum-of-product at these shapes, so the
        # numerator stays as-is regardless of variance_adjust.
        numerator = xp.sum(branch_lengths[None, :] * d, axis=-1)

        if normalized:
            result = _normalize(xp, numerator, denominator)
        else:
            result = numerator

        offset = _condensed_offset(i, n_samples)
        out[offset : offset + result.shape[0]] = result

    return out


def _unweighted_unifrac_xp_core(
    xp,
    counts_by_node,
    sample_totals,
    branch_lengths,
    normalized,
    variance_adjust,
):
    n_samples = counts_by_node.shape[0]
    n_pairs = n_samples * (n_samples - 1) // 2
    out = xp.zeros(n_pairs, dtype=branch_lengths.dtype)

    if not variance_adjust:
        # Presence/absence per (sample, node). The per-pair observed/differs
        # branch-length sums are then separable across samples:
        #   A_i = sum_node(U_i * bl)                  (per-sample total)
        #   G_ij = sum_node(U_i * U_j * bl)            (shared-node total)
        #   observed_ij = A_i + A_j - G_ij             (union)
        #   differs_ij  = A_i + A_j - 2 * G_ij         (symmetric difference)
        # so both reduce to two matrix products computed once, instead of
        # the O(n_nodes) work above being redone inside the pairwise loop.
        U = xp.astype(counts_by_node > 0, branch_lengths.dtype)
        A = U @ branch_lengths  # (n_samples,)
        # (n_samples, n_samples)
        G = (U * branch_lengths[None, :]) @ xp.matrix_transpose(U)

        for i in range(n_samples - 1):
            G_i = G[i, i + 1 :]
            denominator = A[i] + A[i + 1 :] - G_i
            numerator = A[i] + A[i + 1 :] - 2 * G_i

            if normalized:
                result = _normalize(xp, numerator, denominator)
            else:
                result = numerator

            offset = _condensed_offset(i, n_samples)
            out[offset : offset + result.shape[0]] = result

        return out

    for i in range(n_samples - 1):
        p_u = counts_by_node[i, :]
        p_v = counts_by_node[i + 1 :, :]
        present_u = p_u > 0
        present_v = p_v > 0
        observed = present_u[None, :] | present_v
        differs = present_u[None, :] != present_v

        inv_vaw, _ = _variance_adjust_terms(
            xp,
            p_u,
            p_v,
            sample_totals,
            i,
        )
        # inv_vaw is already 0 where invalid, so gating on `observed` alone
        # (instead of the old `observed & vaw_valid`) is enough.
        length_term = xp.where(
            observed,
            branch_lengths[None, :] * inv_vaw,
            xp.asarray(0.0, dtype=branch_lengths.dtype),
        )

        numerator = xp.sum(
            xp.where(differs, length_term, xp.asarray(0.0, dtype=length_term.dtype)),
            axis=-1,
        )
        denominator = xp.sum(length_term, axis=-1)

        if normalized:
            result = _normalize(xp, numerator, denominator)
        else:
            result = numerator

        offset = _condensed_offset(i, n_samples)
        out[offset : offset + result.shape[0]] = result

    return out


def _generalized_unifrac_xp_core(
    xp,
    proportions,
    counts_by_node,
    sample_totals,
    branch_lengths,
    alpha,
    variance_adjust,
):
    n_samples = proportions.shape[0]
    n_pairs = n_samples * (n_samples - 1) // 2
    out = xp.zeros(n_pairs, dtype=proportions.dtype)

    for i in range(n_samples - 1):
        p_u = proportions[i, :]
        p_v = proportions[i + 1 :, :]
        s = p_u[None, :] + p_v
        d = xp.abs(p_u[None, :] - p_v)

        if variance_adjust:
            inv_vaw, _ = _variance_adjust_terms(
                xp,
                counts_by_node[i, :],
                counts_by_node[i + 1 :, :],
                sample_totals,
                i,
            )
            # inv_vaw is already 0 where invalid, so multiplying by it (in
            # place of the old divide-through-_safe_divide) zeroes those
            # terms out directly; `s != 0` below then already captures
            # both "invalid" and "valid but s was 0" the same way the old
            # `vaw_valid & (s != 0)` gate did.
            s = s * inv_vaw
            d = d * inv_vaw

        gate = s != 0

        safe_s = xp.where(gate, s, xp.asarray(1.0, dtype=s.dtype))
        sum_pow = xp.where(
            gate,
            branch_lengths[None, :] * safe_s**alpha,
            xp.asarray(0.0, dtype=s.dtype),
        )
        # sum_pow is already exactly 0 outside `gate`, and d / safe_s is finite
        # everywhere, so the gated nodes drop out of this sum on their own.
        numerator = xp.sum(sum_pow * (d / safe_s), axis=-1)
        denominator = xp.sum(sum_pow, axis=-1)

        result = _normalize(xp, numerator, denominator)

        offset = _condensed_offset(i, n_samples)
        out[offset : offset + result.shape[0]] = result

    return out


def _prepare(counts, taxa, tree, validate, need_proportions=True):
    """Shared setup: node counts, tip totals, (optionally) proportions.

    ``need_proportions=False`` skips building the ``(n_samples, n_nodes)``
    proportions matrix for callers that only need ``counts_by_node`` (e.g.
    unweighted UniFrac, which works on presence/absence), returning ``None``
    in its place.
    """
    from skbio.diversity.beta._unifrac import _setup_multiple_unifrac, _get_tip_indices

    counts_by_node, tree_index, branch_lengths = _setup_multiple_unifrac(
        counts, taxa, tree, validate
    )
    counts_by_node = np.ascontiguousarray(counts_by_node, dtype=np.float64)
    branch_lengths = np.ascontiguousarray(branch_lengths, dtype=np.float64)
    tip_indices = _get_tip_indices(tree_index)
    sample_totals = counts_by_node[:, tip_indices].sum(axis=1)

    xp, counts_by_node, sample_totals, branch_lengths = ingest_array(
        counts_by_node, sample_totals, branch_lengths
    )
    if need_proportions:
        totals = sample_totals[:, None]
        proportions = _safe_divide(xp, counts_by_node, totals, totals > 0)
    else:
        proportions = None
    return xp, counts_by_node, sample_totals, branch_lengths, proportions


def _weighted_unifrac_xp(
    counts, taxa, tree, normalized, variance_adjust=False, validate=True
):
    """Compute the condensed weighted UniFrac distance vector, array-API-generic.

    Array-API-generic counterpart of
    :func:`skbio.diversity.beta._unifrac_gpu._weighted_unifrac_gpu`: computes
    the same quantity without a GPU backend (``numba-cuda``/``numba.hip``),
    working on any array-API-compliant array (including plain NumPy).

    """
    xp, counts_by_node, sample_totals, branch_lengths, proportions = _prepare(
        counts, taxa, tree, validate
    )
    return _to_numpy(
        _weighted_unifrac_xp_core(
            xp,
            proportions,
            counts_by_node,
            sample_totals,
            branch_lengths,
            normalized,
            variance_adjust,
        )
    )


def _unweighted_unifrac_xp(
    counts, taxa, tree, normalized, variance_adjust=False, validate=True
):
    """Compute the condensed unweighted UniFrac distance vector, array-API-generic.

    Array-API-generic counterpart of
    :func:`skbio.diversity.beta._unifrac_gpu._unweighted_unifrac_gpu`.

    """
    xp, counts_by_node, sample_totals, branch_lengths, _ = _prepare(
        counts, taxa, tree, validate, need_proportions=False
    )
    return _to_numpy(
        _unweighted_unifrac_xp_core(
            xp,
            counts_by_node,
            sample_totals,
            branch_lengths,
            normalized,
            variance_adjust,
        )
    )


def _generalized_unifrac_xp(
    counts, taxa, tree, alpha, variance_adjust=False, validate=True
):
    """Compute the condensed generalized UniFrac distance vector, array-API-generic.

    Array-API-generic counterpart of
    :func:`skbio.diversity.beta._unifrac_gpu._generalized_unifrac_gpu`.

    """
    xp, counts_by_node, sample_totals, branch_lengths, proportions = _prepare(
        counts, taxa, tree, validate
    )
    return _to_numpy(
        _generalized_unifrac_xp_core(
            xp,
            proportions,
            counts_by_node,
            sample_totals,
            branch_lengths,
            alpha,
            variance_adjust,
        )
    )
