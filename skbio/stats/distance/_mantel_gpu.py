# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Mantel GPU backend: a fused single-source Numba CUDA/HIP kernel.

Analogous to the PERMANOVA GPU backend: it is taken when both distance matrices
are resident on a matching GPU device and ``engine="numba"`` is requested (CUDA
CuPy / PyTorch compile through ``numba.cuda``, ROCm CuPy / PyTorch through
``numba.hip``).
One thread block per permutation walks the upper triangle and accumulates the
permuted Pearson correlation, reducing across the block in shared memory. The
statistic and p-value are assembled on the host in the same RNG order as the
NumPy and array-API paths, so the result is identical to the cython, numba and
xp engines. Spearman uses the same kernel on rank-transformed matrices.

If the fused kernel cannot build or run on the current stack (see :mod:`._gpu`,
notably the AMD/ROCm init-order case), the caller catches it and falls back to
the array-API GPU path. NVIDIA is unaffected.
"""

import numpy as np

from warnings import warn

from scipy.stats import ConstantInputWarning, NearConstantInputWarning

from skbio.util._array import _get_backend_name
from ._gpu import _TPB, _get_kernel

_kernels = {}  # backend name -> compiled kernel (built on first use)


def _build_kernel(gpu):
    """Compile the permuted-Pearson kernel for ``gpu`` (a Numba GPU module)."""
    from numba import float64 as nb_f64

    @gpu.jit
    def _mantel_r(n_dims, mat, perm_order, ym_norm, mul, add, out):  # pragma: no cover
        # runs on the device; coverage.py cannot instrument compiled PTX
        p = gpu.blockIdx.x              # one block per permutation
        t = gpu.threadIdx.x
        tile = gpu.blockDim.x           # == _TPB
        acc_arr = gpu.shared.array(_TPB, nb_f64)

        # Each thread sums the contribution of a strided subset of rows; the row's
        # upper-triangle entries pair the (permuted) x value with the fixed y value
        # at the matching condensed index, normalized on the fly via mul/add.
        # Only the accumulator is float64; the per-element products stay in the
        # matrix dtype (float32 on consumer GPUs, whose fp64 throughput is poor).
        acc = nb_f64(0.0)
        row = t
        while row < n_dims - 1:
            vrow = perm_order[p, row]
            row_start = row * (n_dims - 1) - ((row - 1) * row) // 2
            icol = 0
            col = row + 1
            while col < n_dims:
                vcol = perm_order[p, col]
                y = ym_norm[row_start + icol]
                x = mat[vrow, vcol] * mul + add
                acc += y * x
                col += 1
                icol += 1
            row += tile

        # shared-memory tree reduction of acc across the block (128 -> 32 -> 8 -> 1)
        acc_arr[t] = acc
        gpu.syncthreads()
        if t < 32:
            s = 0.0
            k = t
            while k < tile:
                s += acc_arr[k]
                k += 32
            acc_arr[t] = s
        gpu.syncthreads()
        if t < 8:
            s = 0.0
            k = t
            while k < 32:
                s += acc_arr[k]
                k += 8
            acc_arr[t] = s
        gpu.syncthreads()
        if t == 0:
            s = 0.0
            for k in range(8):
                s += acc_arr[k]
            if s > 1.0:
                s = 1.0
            elif s < -1.0:
                s = -1.0
            out[p] = s

    return _mantel_r


def _perm_order(n, permutations, rng):
    """Identity permutation followed by ``permutations`` permutations of range(n).

    Consumes ``rng`` in the same order as ``_mantel_stats_pearson_flat`` and
    ``_mantel_stats_pearson_xp`` (identity first, then ``permutations`` calls to
    ``rng.permutation(n)``), so the p-value matches the other engines exactly.
    """
    out = np.empty((permutations + 1, n), dtype=np.int64)
    out[0] = np.arange(n)
    for r in range(1, permutations + 1):
        out[r] = rng.permutation(n)
    return out


def _run_mantel_gpu(gpu, x, y, permutations, rng, alternative, spearman=False):
    """Run the Mantel pearson/spearman test with the fused GPU kernel.

    Mirrors :func:`._mantel._mantel_stats_pearson_xp`'s preprocessing (upper
    triangle, normalization, optional rank transform) but replaces the
    per-permutation array-API loop with the fused kernel on the device-resident
    matrices. ``gpu`` is the Numba module from :func:`._gpu._numba_gpu_module_for`.
    Returns ``(orig_stat, p_value, n)`` to match :func:`._mantel.mantel`.
    """
    # Lazy import avoids a circular import (``_mantel`` imports this module).
    from ._mantel import _upper_tri_xp, _xp_rank_average, _device
    from skbio.util._array import ingest_array, _get_array

    xp, X, Y = ingest_array(x, y)
    if X.shape != Y.shape:
        raise ValueError("Distance matrices must have the same shape.")
    if (X.ndim != 2) or (X.shape[0] != X.shape[1]):
        raise ValueError("Distance matrix must be a square 2-D array.")
    n = X.shape[0]
    if n < 3:
        raise ValueError(
            "Distance matrices must have at least 3 matching IDs "
            "between them (i.e., minimum 3x3 in size)."
        )

    iu0, iu1 = _upper_tri_xp(xp, n, device=_device(X))
    x_flat = X[iu0, iu1]
    y_flat = Y[iu0, iu1]

    if spearman:
        x_flat = _xp_rank_average(xp, x_flat)
        y_flat = _xp_rank_average(xp, y_flat)
        # rebuild the ranked symmetric matrix so the kernel can gather permuted
        # pairs from it (one-time host assembly, as in the array-API path).
        xr = _get_array(x_flat, to_numpy=True)
        i0, j0 = _get_array(iu0, to_numpy=True), _get_array(iu1, to_numpy=True)
        full = np.zeros((n, n), dtype=np.float64)
        full[i0, j0] = xr
        full[j0, i0] = xr
        Xsrc = xp.asarray(full, device=_device(X))
    else:
        Xsrc = X

    # constant input -> correlation undefined (matches scipy/skbio behavior)
    if bool(xp.all(x_flat == x_flat[0])) or bool(xp.all(y_flat == y_flat[0])):
        warn(ConstantInputWarning())
        return np.nan, np.nan, n

    xmean = float(xp.mean(x_flat))
    normxm = float(xp.sqrt(xp.sum((x_flat - xmean) ** 2)))
    ymean = float(xp.mean(y_flat))
    normym = float(xp.sqrt(xp.sum((y_flat - ymean) ** 2)))
    ym_norm = (y_flat - ymean) / normym

    # near-constant input -> loss of precision in r (matches the NumPy path)
    threshold = 1e-13
    if (normxm < threshold * abs(xmean)) or (normym < threshold * abs(ymean)):
        warn(NearConstantInputWarning())

    orig_stat = float(xp.sum(((x_flat - xmean) / normxm) * ym_norm))
    orig_stat = max(min(orig_stat, 1.0), -1.0)

    if permutations == 0 or np.isnan(orig_stat):
        return orig_stat, np.nan, n

    mul = 1.0 / normxm
    add = -xmean / normxm
    # Keep the kernel's per-element math in the matrix dtype (float32 on consumer
    # GPUs); only its accumulator is float64. Match the scalars and y to Xsrc.
    if Xsrc.dtype == xp.float32:
        mul, add = np.float32(mul), np.float32(add)
        ym_norm = xp.astype(ym_norm, xp.float32)
    perm_order = _perm_order(n, permutations, rng)

    kernel = _get_kernel(_kernels, _build_kernel, gpu, _get_backend_name(xp))
    d_perm = gpu.to_device(perm_order)
    d_out = gpu.device_array(permutations + 1, dtype=np.float64)
    # Xsrc and ym_norm are already on the device; the kernel reads them in place
    # via their __cuda_array_interface__ (no host round-trip).
    kernel[permutations + 1, _TPB](n, Xsrc, d_perm, ym_norm, mul, add, d_out)
    gpu.synchronize()

    stats = d_out.copy_to_host()
    comp_stat = stats[0]        # identity permutation, matches perm_order[0]
    permuted_stats = stats[1:]

    if alternative == "two-sided":
        count_better = (np.absolute(permuted_stats) >= np.absolute(comp_stat)).sum()
    elif alternative == "greater":
        count_better = (permuted_stats >= comp_stat).sum()
    else:
        count_better = (permuted_stats <= comp_stat).sum()
    p_value = (count_better + 1) / (permutations + 1)

    return orig_stat, p_value, n
