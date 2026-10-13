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
keeping peak memory at ``O(n_samples, n_nodes)`` per iteration.
"""

import numpy as np

from skbio.util._array import ingest_array, _to_numpy


def _condensed_offset(i, n_samples):
    """Flat offset such that out[offset + k] is pair (i, i + 1 + k)."""
    return i * n_samples - (i * (i + 1)) // 2


def _variance_adjust_terms(xp, counts_i, counts_block, sample_totals, i):
    """Compute ``vaw`` and a validity mask for row ``i`` against ``j > i``.

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
    vaw : array, shape (n_remaining, n_nodes)
        ``sqrt(mi * (m - mi))``, undefined (but finite) where invalid.
    valid : array, shape (n_remaining, n_nodes), bool
        Where ``vaw > 0``.

    """
    m = sample_totals[i] + sample_totals[i + 1 :]  # (n_remaining,)
    mi = counts_i[None, :] + counts_block  # (n_remaining, n_nodes)
    arg = mi * (m[:, None] - mi)
    valid = arg > 0
    safe_arg = xp.where(valid, arg, xp.asarray(0.0, dtype=arg.dtype))
    vaw = xp.sqrt(safe_arg)
    return vaw, valid


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
    xp, proportions, counts_by_node, sample_totals, branch_lengths, normalized,
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
            vaw, valid = _variance_adjust_terms(
                xp, counts_by_node[i, :], counts_by_node[i + 1 :, :],
                sample_totals, i,
            )
            s = _safe_divide(xp, s, vaw, valid)
            d = _safe_divide(xp, d, vaw, valid)

        numerator = xp.sum(branch_lengths[None, :] * d, axis=-1)
        denominator = xp.sum(branch_lengths[None, :] * s, axis=-1)

        if normalized:
            result = _normalize(xp, numerator, denominator)
        else:
            result = numerator

        offset = _condensed_offset(i, n_samples)
        out[offset : offset + result.shape[0]] = result

    return out


def _unweighted_unifrac_xp_core(
    xp, counts_by_node, sample_totals, branch_lengths, normalized, variance_adjust,
):
    n_samples = counts_by_node.shape[0]
    n_pairs = n_samples * (n_samples - 1) // 2
    out = xp.zeros(n_pairs, dtype=branch_lengths.dtype)

    for i in range(n_samples - 1):
        p_u = counts_by_node[i, :]
        p_v = counts_by_node[i + 1 :, :]
        present_u = p_u > 0
        present_v = p_v > 0
        observed = present_u[None, :] | present_v
        differs = present_u[None, :] != present_v

        if variance_adjust:
            vaw, vaw_valid = _variance_adjust_terms(
                xp, p_u, p_v, sample_totals, i,
            )
            length_term = _safe_divide(
                xp, branch_lengths[None, :], vaw, observed & vaw_valid
            )
        else:
            length_term = xp.where(
                observed,
                branch_lengths[None, :],
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
    xp, proportions, counts_by_node, sample_totals, branch_lengths, alpha,
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
            vaw, vaw_valid = _variance_adjust_terms(
                xp, counts_by_node[i, :], counts_by_node[i + 1 :, :],
                sample_totals, i,
            )
            s = _safe_divide(xp, s, vaw, vaw_valid)
            d = _safe_divide(xp, d, vaw, vaw_valid)
            gate = vaw_valid & (s != 0)
        else:
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


def _prepare(counts, taxa, tree, validate):
    """Shared setup: node counts, tip totals, proportions, in ``xp`` terms."""
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
    totals = sample_totals[:, None]
    proportions = _safe_divide(xp, counts_by_node, totals, totals > 0)
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
    return _to_numpy(_weighted_unifrac_xp_core(
        xp, proportions, counts_by_node, sample_totals, branch_lengths,
        normalized, variance_adjust,
    ))


def _unweighted_unifrac_xp(
    counts, taxa, tree, normalized, variance_adjust=False, validate=True
):
    """Compute the condensed unweighted UniFrac distance vector, array-API-generic.

    Array-API-generic counterpart of
    :func:`skbio.diversity.beta._unifrac_gpu._unweighted_unifrac_gpu`.

    """
    xp, counts_by_node, sample_totals, branch_lengths, _ = _prepare(
        counts, taxa, tree, validate
    )
    return _to_numpy(_unweighted_unifrac_xp_core(
        xp, counts_by_node, sample_totals, branch_lengths, normalized,
        variance_adjust,
    ))


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
    return _to_numpy(_generalized_unifrac_xp_core(
        xp, proportions, counts_by_node, sample_totals, branch_lengths, alpha,
        variance_adjust,
    ))
