# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""PERMANOVA GPU backend: a fused single-source Numba CUDA/HIP kernel.

This is the fast path taken when a DistanceMatrix is already on a GPU device and
``engine="numba"`` is requested: the same source compiles through ``numba.cuda``
on NVIDIA and ``numba.hip`` on AMD, and the kernel reads the device-resident
matrix in place (no host round-trip). It follows the coalesced design of
scikit-bio-binaries (``pmn_f_stat_sW_cuda_one``): one block per permutation,
threads walk the columns so global reads are coalesced, group labels cached in
shared memory, s_W reduced in shared memory.

Correctness-first stays the rule: this backend runs only when the device maps to
an available Numba GPU module (CUDA CuPy / PyTorch -> ``numba.cuda``, ROCm CuPy /
PyTorch -> ``numba.hip``). If the fused kernel cannot build or run on the current stack,
the caller catches it, marks the backend, and keeps the vendor-neutral xp path.
Permutations are drawn on the host in the same RNG order as the CPU Monte-Carlo
path, so the p-value is identical to the cython, numba and xp engines.

On AMD/ROCm the numba-hip kernel and torch-ROCm can conflict at kernel-compile
time on some stacks (an AMD init-order issue). Where they do, the kernel launch
raises and the matrix falls back to the array-API path; where the stack supports
both, the kernel runs in place. NVIDIA is unaffected.
"""

import numpy as np

from skbio.util import get_rng
from skbio.util._array import _get_backend_name
from ._gpu import _TPB, _numba_gpu_module_for, _get_kernel

_kernels = {}  # backend name -> compiled kernel (built on first use)


def _build_kernel(gpu):
    """Compile the within-group sum-of-squares kernel for ``gpu`` (cuda or hip)."""
    from numba import float64 as nb_f64, int64 as nb_i64

    @gpu.jit
    def _pmn_sW(n_dims, mat, groupings, inv_gs, out_sW):  # pragma: no cover
        # runs on the device; coverage.py cannot instrument compiled PTX
        gel = gpu.blockIdx.x            # one block per permutation
        icol = gpu.threadIdx.x
        tile = gpu.blockDim.x           # == _TPB
        s_W_arr = gpu.shared.array(_TPB, nb_f64)
        row_grouping = gpu.shared.array(_TPB, nb_i64)

        # Only the accumulators are float64; the per-element products stay in the
        # matrix dtype (float32 on consumer GPUs, whose fp64 throughput is poor).
        s_W = nb_f64(0.0)
        trow = 0
        while trow < n_dims - 1:
            tcol = trow + 1
            while tcol < n_dims:
                max_row = min(trow + tile, n_dims - 1)
                max_col = min(tcol + tile, n_dims)
                gpu.syncthreads()
                rc = trow + icol
                if rc < max_row:
                    row_grouping[icol] = groupings[gel, rc]
                gpu.syncthreads()

                local = nb_f64(0.0)
                row = trow
                while row < max_row:
                    min_col = max(tcol, row + 1)
                    gidx = row_grouping[row - trow]
                    col = min_col + icol      # thread -> column: coalesced read
                    if col < max_col:
                        if groupings[gel, col] == gidx:
                            v = mat[row, col]
                            local += v * v * inv_gs[gidx]
                    row += 1
                s_W += local
                tcol += tile
            trow += tile

        # shared-memory tree reduction of s_W across the block (128 -> 32 -> 8 -> 1)
        s_W_arr[icol] = s_W
        gpu.syncthreads()
        if icol < 32:
            acc = 0.0
            t = icol
            while t < tile:
                acc += s_W_arr[t]
                t += 32
            s_W_arr[icol] = acc
        gpu.syncthreads()
        if icol < 8:
            acc = 0.0
            t = icol
            while t < 32:
                acc += s_W_arr[t]
                t += 8
            s_W_arr[icol] = acc
        gpu.syncthreads()
        if icol == 0:
            acc = 0.0
            for t in range(8):
                acc += s_W_arr[t]
            out_sW[gel] = acc

    return _pmn_sW


def _permutation_batch(grouping, permutations, seed):
    """The observed grouping followed by ``permutations`` permutations of it.

    Uses ``get_rng(seed)`` and ``rng.permutation`` in the same order as the CPU
    Monte-Carlo driver (``_run_monte_carlo_stats``), so p-values match exactly.
    """
    rng = get_rng(seed)
    out = np.empty((permutations + 1, grouping.shape[0]), dtype=np.int64)
    out[0] = grouping
    for i in range(permutations):
        out[i + 1] = rng.permutation(grouping)
    return out


def _assemble_fp(s_W, s_T, sample_size, num_groups, permutations):
    """Pseudo-F and Monte Carlo p-value from the per-permutation ``s_W``.

    ``s_W[0]`` is the observed grouping; ``s_W[1:]`` the permutations. Pure NumPy
    (no GPU), so the assembly and p-value formula are covered by ordinary CI and
    match the cython/numba/xp result.
    """
    n, ng = sample_size, num_groups
    f = ((s_T - s_W) / (ng - 1)) / (s_W / (n - ng))
    obs = float(f[0])
    if permutations > 0:
        p_value = ((f[1:] >= f[0]).sum() + 1) / (permutations + 1)
    else:
        p_value = np.nan
    return obs, float(p_value)


def _run_permanova_gpu(gpu, distmat, grouping, column, permutations, seed, ids=None):
    """Run PERMANOVA with the fused GPU kernel on a device-resident matrix.

    Mirrors :func:`_permanova_array_api`'s preprocessing (same ``ids`` alignment
    and ``s_T``), then runs the fused kernel on the on-device matrix instead of
    the per-permutation array-API loop. ``gpu`` is the Numba module returned by
    :func:`_numba_gpu_module_for`.
    """
    from ._base import _preprocess_input_sng, _build_results
    from skbio.util._array import ingest_array

    xp, dm = ingest_array(distmat)
    sample_size = dm.shape[0]
    num_groups, grouping = _preprocess_input_sng(ids, sample_size, grouping, column)

    group_sizes = np.bincount(grouping)
    # Match the matrix dtype so the kernel's per-element math stays in that
    # precision (float32 on consumer GPUs); the kernel accumulates in float64.
    elem_dtype = np.float32 if dm.dtype == xp.float32 else np.float64
    inv_gs = (1.0 / group_sizes).astype(elem_dtype)
    # full 2-D matrix (array-API input is never condensed); halve to count each
    # unordered pair once, matching the DistanceMatrix full-matrix path.
    s_T = float(xp.sum(dm * dm, dtype=xp.float64) / sample_size / 2.0)

    groupings = _permutation_batch(grouping, permutations, seed)

    kernel = _get_kernel(_kernels, _build_kernel, gpu, _get_backend_name(xp))
    d_grp = gpu.to_device(groupings)
    d_inv = gpu.to_device(inv_gs)
    d_out = gpu.device_array(permutations + 1, dtype=np.float64)
    # dm is already on the device; the kernel reads it in place via the array's
    # __cuda_array_interface__ (no host round-trip).
    kernel[permutations + 1, _TPB](sample_size, dm, d_grp, d_inv, d_out)
    gpu.synchronize()

    stat, p_value = _assemble_fp(
        d_out.copy_to_host(), s_T, sample_size, num_groups, permutations
    )
    return _build_results(
        "PERMANOVA", "pseudo-F", sample_size, num_groups, stat, p_value, permutations
    )
