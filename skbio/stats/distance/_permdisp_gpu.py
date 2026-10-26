# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""PERMDISP GPU backend: a fused single-source Numba CUDA/HIP kernel.

This is the fast path taken when a DistanceMatrix is already on a GPU device and
``engine="numba"`` is requested: the same source compiles through
``numba_cuda_mlir.cuda`` or ``numba.cuda`` on NVIDIA and ``numba.hip`` on AMD. pcoa
returns its coordinates to the host either way, so the kernel reads that much
smaller array rather than the distance matrix itself, and only the per-permutation
statistics come back.

One block owns one (permutation, group) pair and reduces over that group's
samples in shared memory, so no thread holds a per-group or per-dimension array.

The median test runs the same modified Weiszfeld iteration as
``_cutils.geomedian_axis_one`` and the CPU Numba port, including breaking before
adopting the new estimate, so all three stop on the same iterate. Permutations
are drawn on the host in the same RNG order as the CPU paths.

This backend runs only when the device maps to an available Numba GPU module. If
the kernel cannot build or run, the caller marks the backend unavailable and keeps
the host Numba engine.
"""

import math

import numpy as np

from skbio.util import get_rng
from ._gpu import _get_kernel

_kernels = {}  # backend name -> compiled kernel (built on first use)

# _MAX_DIMS bounds the shared arrays the kernels declare, which need a
# compile-time constant shape on both numba.cuda and numba.hip. _MAX_GROUPS is
# a limit on the grid's group axis, not on a kernel array. Both are checked on
# the host before launch. 64 dimensions cover the 10 axes permdisp retains by
# default, and up to about 2000 would still fit the median kernel in 48 KB of
# shared memory. The group limit was not tuned. Each (permutation, group) block
# scans all samples, so GPU time grows with the number of groups much faster
# than the CPU engine's.
_MAX_GROUPS = 64
_MAX_DIMS = 64


def _kernel_supports_shape(dims, num_groups):
    """Whether the fused kernel can take an ordination of this shape.

    True if ``dims`` and ``num_groups`` are within ``_MAX_DIMS`` and
    ``_MAX_GROUPS``.
    """
    return dims <= _MAX_DIMS and num_groups <= _MAX_GROUPS


# Upper bound on the grouping buffer held at once, host side and device side,
# which would otherwise be (permutations + 1) x n. The result does not depend on
# it, because the permutations are drawn in the same order whatever the chunk
# size. It is a memory cap, not a tuned value.
_PERM_CHUNK_BYTES = 64 * 1024 * 1024


# Threads per block, one block per (permutation, group). It sizes the shared
# reduction array and the reductions halve it, so it must be a compile-time power
# of two. 64 sits in the middle of the range that ran within 2% on an RTX 4060
# (32 to 128), and 256 was 9% slower on the median test.
_TPB_BG = 64


def _build_centroid_kernel_bg(gpu):
    """Compile the block-per-(permutation, group) centroid kernel.

    Dimensions are looped on the outside, one shared-memory reduction each, so
    a thread only holds a scalar. ``samples`` is ``(dims, n)`` so the threads of
    a reduction read consecutive addresses. Each block writes its own slot of
    the outputs, so no atomics are needed; the host turns the partials into F.
    """
    # numpy's float64 rather than numba's, which numba-cuda-mlir does not accept.
    nb_f64 = np.float64

    @gpu.jit
    def _pd_centroid_bg(
        samples, n, dims, groupings, out_count, out_sum, out_ss
    ):  # pragma: no cover
        # runs on the device; coverage.py cannot instrument compiled PTX
        p = gpu.blockIdx.x
        g = gpu.blockIdx.y
        t = gpu.threadIdx.x
        nt = gpu.blockDim.x

        red = gpu.shared.array(_TPB_BG, nb_f64)
        centroid = gpu.shared.array(_MAX_DIMS, nb_f64)

        # count of samples in this (permutation, group).
        c = nb_f64(0.0)
        for i in range(t, n, nt):
            if groupings[p, i] == g:
                c += 1.0
        red[t] = c
        gpu.syncthreads()
        stride = nt // 2
        while stride > 0:
            if t < stride:
                red[t] += red[t + stride]
            gpu.syncthreads()
            stride //= 2
        count = red[0]
        gpu.syncthreads()
        if count == 0.0:
            if t == 0:
                out_count[p, g] = 0.0
                out_sum[p, g] = 0.0
                out_ss[p, g] = 0.0
            return

        # centroid, one reduction per dimension.
        for j in range(dims):
            acc = nb_f64(0.0)
            for i in range(t, n, nt):
                if groupings[p, i] == g:
                    acc += samples[j, i]
            red[t] = acc
            gpu.syncthreads()
            stride = nt // 2
            while stride > 0:
                if t < stride:
                    red[t] += red[t + stride]
                gpu.syncthreads()
                stride //= 2
            if t == 0:
                centroid[j] = red[0] / count
            gpu.syncthreads()

        # sum of distances to the centroid.
        acc = nb_f64(0.0)
        for i in range(t, n, nt):
            if groupings[p, i] == g:
                d2 = nb_f64(0.0)
                for j in range(dims):
                    diff = samples[j, i] - centroid[j]
                    d2 += diff * diff
                acc += math.sqrt(d2)
        red[t] = acc
        gpu.syncthreads()
        stride = nt // 2
        while stride > 0:
            if t < stride:
                red[t] += red[t + stride]
            gpu.syncthreads()
            stride //= 2
        sum_dist = red[0]
        gpu.syncthreads()
        mean_g = sum_dist / count

        # within-group sum of squared residuals.
        acc = nb_f64(0.0)
        for i in range(t, n, nt):
            if groupings[p, i] == g:
                d2 = nb_f64(0.0)
                for j in range(dims):
                    diff = samples[j, i] - centroid[j]
                    d2 += diff * diff
                resid = math.sqrt(d2) - mean_g
                acc += resid * resid
        red[t] = acc
        gpu.syncthreads()
        stride = nt // 2
        while stride > 0:
            if t < stride:
                red[t] += red[t + stride]
            gpu.syncthreads()
            stride //= 2

        if t == 0:
            out_count[p, g] = count
            out_sum[p, g] = sum_dist
            out_ss[p, g] = red[0]

    return _pd_centroid_bg


def _assemble_f_from_partials(counts, sums, ss, num_groups):
    """Turn per-(permutation, group) partials into one F per permutation.

    Mirrors the CPU driver: the grand mean is over samples rather than over
    groups, so it is the total distance divided by the total count.
    """
    total = counts.sum(axis=1)
    grand = sums.sum(axis=1) / total
    with np.errstate(invalid="ignore", divide="ignore"):
        means = np.where(counts > 0, sums / np.maximum(counts, 1.0), 0.0)
        ss_between = (counts * (means - grand[:, None]) ** 2).sum(axis=1)
        ss_within = ss.sum(axis=1)
        ms_between = ss_between / (num_groups - 1)
        ms_within = ss_within / (total - num_groups)
        out = np.where(
            ms_within == 0.0,
            np.where(ss_between == 0.0, np.nan, np.inf),
            ms_between / np.where(ms_within == 0.0, 1.0, ms_within),
        )
    return out


# Tolerance and iteration cap of the Weiszfeld iteration. They are the
# defaults of _cutils.geomedian_axis_one and the CPU Numba port, and must stay
# equal to them for the three engines to stop on the same iterate.
_GEO_EPS = 1e-7
_GEO_MAXITERS = 500

# Dimension-tile width for the Weiszfeld weighted sum. Each thread holds this
# many float64 accumulators, so it trades registers against repeated distance
# work. The best width follows the number of dimensions (8 at 5, 16 at 10, 32 at
# 20 on an RTX 4060), and 16 suits the 10 axes permdisp retains by default.
_JCH = 16


def _build_median_kernel_bg(gpu):
    """Compile the block-per-(permutation, group) spatial-median kernel.

    Same shape as the centroid kernel. The geometric median is the modified
    Weiszfeld iteration of ``_geomedian_nb``, with the same ``eps``,
    ``maxiters`` and stopping rule, so the GPU and CPU paths agree. The block
    exits on convergence rather than running a fixed number of iterations.
    """
    # numpy's float64 rather than numba's, which numba-cuda-mlir does not accept.
    nb_f64 = np.float64

    @gpu.jit
    def _pd_median_bg(
        samples, n, dims, groupings, out_count, out_sum, out_ss
    ):  # pragma: no cover
        # runs on the device; coverage.py cannot instrument compiled PTX
        p = gpu.blockIdx.x
        g = gpu.blockIdx.y
        t = gpu.threadIdx.x
        nt = gpu.blockDim.x

        red = gpu.shared.array(_TPB_BG, nb_f64)
        y = gpu.shared.array(_MAX_DIMS, nb_f64)
        y1 = gpu.shared.array(_MAX_DIMS, nb_f64)
        T = gpu.shared.array(_MAX_DIMS, nb_f64)
        flag = gpu.shared.array(2, nb_f64)  # 0: Dinvs, 1: nzeros

        # count.
        c = nb_f64(0.0)
        for i in range(t, n, nt):
            if groupings[p, i] == g:
                c += 1.0
        red[t] = c
        gpu.syncthreads()
        st = nt // 2
        while st > 0:
            if t < st:
                red[t] += red[t + st]
            gpu.syncthreads()
            st //= 2
        count = red[0]
        gpu.syncthreads()
        if count == 0.0:
            if t == 0:
                out_count[p, g] = 0.0
                out_sum[p, g] = 0.0
                out_ss[p, g] = 0.0
            return

        # initial estimate: the mean.
        for j in range(dims):
            acc = nb_f64(0.0)
            for i in range(t, n, nt):
                if groupings[p, i] == g:
                    acc += samples[j, i]
            red[t] = acc
            gpu.syncthreads()
            st = nt // 2
            while st > 0:
                if t < st:
                    red[t] += red[t + st]
                gpu.syncthreads()
                st //= 2
            if t == 0:
                y[j] = red[0] / count
            gpu.syncthreads()

        # modified Weiszfeld.
        if count > 1.0:
            for _ in range(_GEO_MAXITERS):
                # inverse distances and the count of samples on the estimate
                sdinv = nb_f64(0.0)
                nz = nb_f64(0.0)
                for i in range(t, n, nt):
                    if groupings[p, i] == g:
                        d2 = nb_f64(0.0)
                        for j in range(dims):
                            diff = samples[j, i] - y[j]
                            d2 += diff * diff
                        di = math.sqrt(d2)
                        if di > _GEO_EPS:
                            sdinv += 1.0 / di
                        else:
                            nz += 1.0
                red[t] = sdinv
                gpu.syncthreads()
                st = nt // 2
                while st > 0:
                    if t < st:
                        red[t] += red[t + st]
                    gpu.syncthreads()
                    st //= 2
                if t == 0:
                    flag[0] = red[0]
                gpu.syncthreads()
                red[t] = nz
                gpu.syncthreads()
                st = nt // 2
                while st > 0:
                    if t < st:
                        red[t] += red[t + st]
                    gpu.syncthreads()
                    st //= 2
                if t == 0:
                    flag[1] = red[0]
                gpu.syncthreads()
                dinvs = flag[0]
                nzeros = flag[1]

                # every sample sits on the estimate: it is the median
                if nzeros == count:
                    break

                # Weighted sum over the samples that are off the estimate.
                # The weight 1/di depends on the sample, not on the dimension,
                # so the dimension loop is tiled _JCH wide and the distance is
                # computed once per tile rather than once per dimension.
                for j0 in range(0, dims, _JCH):
                    accs = gpu.local.array(_JCH, nb_f64)
                    for jj in range(_JCH):
                        accs[jj] = nb_f64(0.0)
                    for i in range(t, n, nt):
                        if groupings[p, i] == g:
                            d2 = nb_f64(0.0)
                            for k in range(dims):
                                diff = samples[k, i] - y[k]
                                d2 += diff * diff
                            di = math.sqrt(d2)
                            if di > _GEO_EPS:
                                w = 1.0 / di
                                for jj in range(_JCH):
                                    if j0 + jj < dims:
                                        accs[jj] += w * samples[j0 + jj, i]
                    # One reduction per dimension of the tile. The tail tile
                    # still reduces every slot so that no thread skips a
                    # syncthreads, and only the in-range ones are written.
                    for jj in range(_JCH):
                        red[t] = accs[jj]
                        gpu.syncthreads()
                        st = nt // 2
                        while st > 0:
                            if t < st:
                                red[t] += red[t + st]
                            gpu.syncthreads()
                            st //= 2
                        if t == 0 and j0 + jj < dims:
                            T[j0 + jj] = red[0] / dinvs
                        gpu.syncthreads()

                if t == 0:
                    if nzeros == 0.0:
                        for j in range(dims):
                            y1[j] = T[j]
                    else:
                        r2 = nb_f64(0.0)
                        for j in range(dims):
                            rj = (T[j] - y[j]) * dinvs
                            r2 += rj * rj
                        r = math.sqrt(r2)
                        rinv = nzeros / r if r > _GEO_EPS else 0.0
                        a = 1.0 - rinv
                        if a < 0.0:
                            a = 0.0
                        b = rinv
                        if b > 1.0:
                            b = 1.0
                        for j in range(dims):
                            y1[j] = a * T[j] + b * y[j]
                    mv = nb_f64(0.0)
                    for j in range(dims):
                        diff = y[j] - y1[j]
                        mv += diff * diff
                    flag[0] = math.sqrt(mv)
                gpu.syncthreads()
                # Break before adopting y1. hdmedians and both CPU ports return
                # the previous estimate on the converging iteration.
                if flag[0] < _GEO_EPS:
                    break
                if t == 0:
                    for j in range(dims):
                        y[j] = y1[j]
                gpu.syncthreads()

        # distances to the median, then the within-group residuals.
        acc = nb_f64(0.0)
        for i in range(t, n, nt):
            if groupings[p, i] == g:
                d2 = nb_f64(0.0)
                for j in range(dims):
                    diff = samples[j, i] - y[j]
                    d2 += diff * diff
                acc += math.sqrt(d2)
        red[t] = acc
        gpu.syncthreads()
        st = nt // 2
        while st > 0:
            if t < st:
                red[t] += red[t + st]
            gpu.syncthreads()
            st //= 2
        sum_dist = red[0]
        gpu.syncthreads()
        mean_g = sum_dist / count

        acc = nb_f64(0.0)
        for i in range(t, n, nt):
            if groupings[p, i] == g:
                d2 = nb_f64(0.0)
                for j in range(dims):
                    diff = samples[j, i] - y[j]
                    d2 += diff * diff
                resid = math.sqrt(d2) - mean_g
                acc += resid * resid
        red[t] = acc
        gpu.syncthreads()
        st = nt // 2
        while st > 0:
            if t < st:
                red[t] += red[t + st]
            gpu.syncthreads()
            st //= 2

        if t == 0:
            out_count[p, g] = count
            out_sum[p, g] = sum_dist
            out_ss[p, g] = red[0]

    return _pd_median_bg


def _permutation_chunks(codes, permutations, seed, chunk):
    """Yield ``(start, block)`` covering the observed grouping then permutations.

    ``block[k]`` is the grouping for permutation ``start + k``; row 0 of the
    first block is the observed grouping. Permutations are drawn in the same
    order as :func:`._permdisp._run_permdisp_numba` whatever ``chunk`` is, so
    the p-value does not depend on the chunk size.
    """
    rng = get_rng(seed)
    n_total = permutations + 1
    start = 0
    while start < n_total:
        size = min(chunk, n_total - start)
        block = np.empty((size, codes.shape[0]), dtype=codes.dtype)
        for k in range(size):
            block[k] = codes if start + k == 0 else rng.permutation(codes)
        yield start, block
        start += size


def _assemble_stat_p(stats, permutations):
    """Observed statistic and Monte Carlo p-value from the per-permutation array."""
    stat = float(stats[0])
    if permutations == 0:
        return stat, np.nan
    p_value = (1.0 + np.sum(stats[1:] >= stats[0])) / (1.0 + permutations)
    return stat, float(p_value)


# Separate from _kernels: _get_kernel's cache key is only the backend name, so
# one shared dict would return whichever kernel was compiled first.
_median_kernels = {}


def _run_permdisp_gpu(
    gpu, builder, cache, sample_data, codes, num_groups, permutations, seed
):
    """Shared driver for the block-per-(permutation, group) kernels.

    ``sample_data`` is the host-side ordination, transposed to ``(dims, n)``
    and uploaded once. ``builder`` and ``cache`` select which statistic's
    kernel to compile and reuse.
    """
    n_samples, dims = sample_data.shape
    if num_groups > _MAX_GROUPS:
        raise ValueError(
            f"the permdisp GPU kernel supports at most {_MAX_GROUPS} groups; "
            f"got {num_groups}. Callers that dispatch automatically check "
            "_kernel_supports_shape first and take the host path instead."
        )
    if dims > _MAX_DIMS:
        raise ValueError(
            f"the permdisp GPU kernel supports at most {_MAX_DIMS} retained "
            f"dimensions; got {dims}. Callers that dispatch automatically "
            "check _kernel_supports_shape first and take the host path instead."
        )

    samples = np.ascontiguousarray(sample_data.T, dtype=np.float64)
    # Labels index groups and are never accumulated, and num_groups is already
    # bounded at _MAX_GROUPS above, but the kernel compares groupings against
    # gpu.blockIdx.y, an int32 on both cuda and hip, so the host array is int32.
    codes32 = np.ascontiguousarray(codes, dtype=np.int32)

    n_total = permutations + 1

    # No array-API namespace here (sample_data is host NumPy, see module
    # docstring), so the cache key is the Numba GPU module itself rather than
    # the array-API backend name _get_kernel's other callers use.
    kernel = _get_kernel(cache, builder, gpu, gpu.__name__)
    d_samples = gpu.to_device(samples)

    chunk = max(1, _PERM_CHUNK_BYTES // max(1, codes32.nbytes))
    stats = np.empty(n_total, dtype=np.float64)
    for start, block in _permutation_chunks(codes32, permutations, seed, chunk):
        size = block.shape[0]
        d_groupings = gpu.to_device(block)
        d_count = gpu.device_array((size, num_groups), dtype=np.float64)
        d_sum = gpu.device_array((size, num_groups), dtype=np.float64)
        d_ss = gpu.device_array((size, num_groups), dtype=np.float64)
        kernel[(size, num_groups), _TPB_BG](
            d_samples, n_samples, dims, d_groupings, d_count, d_sum, d_ss
        )
        gpu.synchronize()
        stats[start : start + size] = _assemble_f_from_partials(
            d_count.copy_to_host(),
            d_sum.copy_to_host(),
            d_ss.copy_to_host(),
            num_groups,
        )

    return _assemble_stat_p(stats, permutations)


def _run_permdisp_centroid_gpu(gpu, sample_data, codes, num_groups, permutations, seed):
    """Run permdisp's centroid test with the fused GPU kernel."""
    return _run_permdisp_gpu(
        gpu,
        _build_centroid_kernel_bg,
        _kernels,
        sample_data,
        codes,
        num_groups,
        permutations,
        seed,
    )


def _run_permdisp_median_gpu(gpu, sample_data, codes, num_groups, permutations, seed):
    """Run permdisp's median test with the fused GPU kernel."""
    return _run_permdisp_gpu(
        gpu,
        _build_median_kernel_bg,
        _median_kernels,
        sample_data,
        codes,
        num_groups,
        permutations,
        seed,
    )
