# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""GPU UniFrac (:mod:`skbio.diversity.beta._unifrac_gpu`).

Backend detection (``_detect_gpu_backend``/``_get_cuda_module``), the fused
Numba CUDA/HIP kernel shared by all five UniFrac methods, the runtime probe
that picks a safe tile configuration for it, and the ``*_gpu_or_xp``
wrappers that fall back to :mod:`skbio.diversity.beta._unifrac_xp` when the
fused kernel is unusable.
"""

import math
import threading
from warnings import warn

import numpy as np

UNWEIGHTED = 0
UNWEIGHTED_UNNORMALIZED = 1
WEIGHTED_NORMALIZED = 2
WEIGHTED_UNNORMALIZED = 3
GENERALIZED = 4

_backend_cache = None

# Compiled kernels, keyed by (id() of the cuda-API module, tile, node_chunk,
# variance_adjust). The module itself is stored alongside so the id stays valid.
_KERNEL_CACHE = {}

# Guards the check-compile-store sequence in `_make_unifrac_kernel`.
_KERNEL_CACHE_LOCK = threading.Lock()


def _detect_gpu_backend():
    """Detect which GPU backend, if any, is usable.

    Returns
    -------
    {'cuda', 'hip', None}
        'cuda' if an NVIDIA GPU and ``numba-cuda`` are both available,
        'hip' if no NVIDIA GPU is found but an AMD GPU and
        ``numba.hip`` are both available, otherwise ``None``.

    """
    global _backend_cache
    if _backend_cache is not None:
        return _backend_cache
    try:
        from numba import cuda

        if cuda.is_available():
            _backend_cache = "cuda"
            return _backend_cache
    except ImportError:
        pass
    try:
        import numba.hip as hip

        if hip.is_available():
            _backend_cache = "hip"
            return _backend_cache
    except ImportError:
        pass
    _backend_cache = None
    return _backend_cache


def _get_cuda_module():
    """Return a ``numba.cuda``-API-compatible module for the detected backend.

    Returns
    -------
    module
        ``numba.cuda`` itself on an NVIDIA backend, or ``numba.hip``
        after calling ``pose_as_cuda()`` on an AMD backend, so kernel
        code written against ``cuda.jit``/``cuda.device_array``/etc.
        runs unmodified on either.

    Raises
    ------
    ImportError
        If no usable GPU backend is detected.

    Notes
    -----
    NVIDIA support comes from the ``numba-cuda`` package, installable via
    scikit-bio's ``gpu-nvidia`` extra. AMD support requires ``numba-hip``,
    which is published only on test.pypi and must be pinned to match the
    installed ROCm version (for example ``numba-hip[rocm-7-0-0]==0.1.6`` for
    ROCm 7.0.0), so there is no single pip specification that works across
    ROCm releases and hence no ``gpu-amd`` extra. ``hip-python`` alone is not
    sufficient: it provides low-level HIP bindings, not
    ``numba.hip``/``pose_as_cuda()``.

    """
    backend = _detect_gpu_backend()
    if backend == "cuda":
        from numba import cuda

        return cuda
    if backend == "hip":
        import numba.hip as hip

        hip.pose_as_cuda()
        return hip
    raise ImportError(
        "No usable GPU backend found (checked numba-cuda and numba.hip). "
        "Install one of these optional dependencies to call this function "
        "directly, or use one of the `*_gpu_or_xp` wrappers (what "
        "engine='gpu' calls), which fall back to the array-API "
        "implementation automatically when no GPU backend is available."
    )


def _build_pair_index(n_samples):
    """Precompute condensed-index (i, j) pairs for i < j, scipy pdist order.

    ``np.triu_indices`` walks the upper triangle row by row, which is exactly
    ``scipy.spatial.distance.pdist``'s condensed order.
    """
    pair_i, pair_j = np.triu_indices(n_samples, k=1)
    return pair_i.astype(np.int32), pair_j.astype(np.int32)


def _build_block_pair_index(n_samples, tile):
    """Precompute (block_i, block_j) pairs of TILE x TILE sample blocks.

    Unlike `_build_pair_index`, diagonal blocks (``block_i == block_j``) are
    included (``k=0``, not ``k=1``): a diagonal block still contains valid
    ``i < j`` pairs within itself (just not the whole block), so it cannot be
    skipped the way a fully lower-triangle block can.
    """
    num_blocks = (n_samples + tile - 1) // tile
    block_i, block_j = np.triu_indices(num_blocks, k=0)
    return block_i.astype(np.int32), block_j.astype(np.int32)


def _make_unifrac_kernel(cuda, tile, node_chunk, variance_adjust):
    """Return the UniFrac pair kernel compiled against a cuda-API module.

    Shared by all 5 UniFrac methods, so only one kernel is compiled and
    launched regardless of which method a driver function dispatches.

    ``tile``/``node_chunk`` are closed over as compile-time constants because
    the shared-memory tiles are sized from them; which pair is safe differs
    per device, so it is discovered at runtime (see ``_get_tile_config``).
    ``variance_adjust`` is closed over too (unlike ``method``/``alpha``, which
    stay runtime arguments) so that the ``variance_adjust=False`` variant does
    not declare or fill the ``sh_cnt_i``/``sh_cnt_j`` tiles it never reads.

    ``cuda.jit`` pays full compilation cost (seconds) on every call, so the
    result is memoized in ``_KERNEL_CACHE``.

    """
    cache_key = (id(cuda), tile, node_chunk, variance_adjust)
    with _KERNEL_CACHE_LOCK:
        cached = _KERNEL_CACHE.get(cache_key)
        if cached is not None:
            return cached[1]

        # Imported here, not at module level, so this module stays importable
        # (and the array-API fallback reachable) with no numba installed.
        from numba import float64

        TILE = tile
        NODE_CHUNK = node_chunk
        # Must be a literal int closure variable: `numba.hip`'s
        # `cuda.shared.array` typing cannot constant-fold `NODE_CHUNK + 1`
        # inside the jitted function and raises a TypingError on it (the same
        # source compiles fine under numba-cuda).
        NODE_CHUNK_PADDED = node_chunk + 1
        VARIANCE_ADJUST = variance_adjust

        @cuda.jit
        def _unifrac_block_kernel(
            proportions,
            counts,
            sample_totals,
            branch_lengths,
            method,
            alpha,
            alpha_mode,
            block_i,
            block_j,
            n_samples,
            out,
        ):
            blk = cuda.blockIdx.x
            if blk >= block_i.shape[0]:
                return
            bi = block_i[blk]
            bj = block_j[blk]

            tx = cuda.threadIdx.x
            ty = cuda.threadIdx.y

            i = bi * TILE + tx
            j = bj * TILE + ty
            valid = i < n_samples and j < n_samples and i < j

            n_nodes = proportions.shape[1]

            # The i-tiles are read [tx, col] with tx fast-varying, so an even
            # column stride would hit shared-memory bank conflicts; the extra
            # padding column makes the stride odd. The j-tiles are read
            # [ty, col], warp-uniform, so they need no padding.
            sh_prop_i = cuda.shared.array((TILE, NODE_CHUNK_PADDED), dtype=float64)
            sh_prop_j = cuda.shared.array((TILE, NODE_CHUNK), dtype=float64)
            if VARIANCE_ADJUST:
                sh_cnt_i = cuda.shared.array((TILE, NODE_CHUNK_PADDED), dtype=float64)
                sh_cnt_j = cuda.shared.array((TILE, NODE_CHUNK), dtype=float64)
            sh_branch = cuda.shared.array(NODE_CHUNK, dtype=float64)

            lin_tid = ty * TILE + tx
            n_threads_blk = TILE * TILE
            tile_elems = TILE * NODE_CHUNK

            numerator = 0.0
            denominator = 0.0

            n_chunks = (n_nodes + NODE_CHUNK - 1) // NODE_CHUNK
            for c in range(n_chunks):
                chunk_start = c * NODE_CHUNK
                chunk_len = n_nodes - chunk_start
                if chunk_len > NODE_CHUNK:
                    chunk_len = NODE_CHUNK

                # Every thread in the block, valid or not, helps fill the
                # shared tiles for this node-chunk.
                for e in range(lin_tid, tile_elems, n_threads_blk):
                    row = e // NODE_CHUNK
                    col = e % NODE_CHUNK
                    node = chunk_start + col
                    samp_i = bi * TILE + row
                    samp_j = bj * TILE + row
                    if col < chunk_len and samp_i < n_samples:
                        sh_prop_i[row, col] = proportions[samp_i, node]
                        if VARIANCE_ADJUST:
                            sh_cnt_i[row, col] = counts[samp_i, node]
                    else:
                        sh_prop_i[row, col] = 0.0
                        if VARIANCE_ADJUST:
                            sh_cnt_i[row, col] = 0.0
                    if col < chunk_len and samp_j < n_samples:
                        sh_prop_j[row, col] = proportions[samp_j, node]
                        if VARIANCE_ADJUST:
                            sh_cnt_j[row, col] = counts[samp_j, node]
                    else:
                        sh_prop_j[row, col] = 0.0
                        if VARIANCE_ADJUST:
                            sh_cnt_j[row, col] = 0.0

                for e in range(lin_tid, NODE_CHUNK, n_threads_blk):
                    node = chunk_start + e
                    sh_branch[e] = branch_lengths[node] if e < chunk_len else 0.0

                cuda.syncthreads()

                if valid:
                    for col in range(chunk_len):
                        p_u = sh_prop_i[tx, col]
                        p_v = sh_prop_j[ty, col]
                        length = sh_branch[col]
                        s = p_u + p_v
                        d = abs(p_u - p_v)
                        if VARIANCE_ADJUST:
                            m = sample_totals[i] + sample_totals[j]
                            mi = sh_cnt_i[tx, col] + sh_cnt_j[ty, col]
                            vaw = math.sqrt(mi * (m - mi))
                            if vaw <= 0.0:
                                continue
                            s = s / vaw
                            d = d / vaw
                        if (
                            method == WEIGHTED_NORMALIZED
                            or method == WEIGHTED_UNNORMALIZED
                        ):
                            numerator += length * d
                            denominator += length * s
                        elif method == UNWEIGHTED or method == UNWEIGHTED_UNNORMALIZED:
                            if s <= 0.0:
                                continue
                            observed = p_u > 0.0 or p_v > 0.0
                            differs = (p_u > 0.0) != (p_v > 0.0)
                            if observed:
                                if VARIANCE_ADJUST:
                                    numerator += length / vaw if differs else 0.0
                                    denominator += length / vaw
                                else:
                                    numerator += length if differs else 0.0
                                    denominator += length
                        elif method == GENERALIZED:
                            if s == 0.0:
                                continue
                            # alpha is a runtime float, so `s**alpha` lowers
                            # to a pow() call per node. alpha_mode lets the
                            # three exact alpha values that have closed forms
                            # (1.0, 0.5, 0.0) skip it.
                            if alpha_mode == 1:
                                sum_pow = length * s
                            elif alpha_mode == 2:
                                sum_pow = length * math.sqrt(s)
                            elif alpha_mode == 3:
                                sum_pow = length
                            else:
                                sum_pow = length * s**alpha
                            numerator += sum_pow * (d / s)
                            denominator += sum_pow

                # Must finish before the next iteration overwrites the shared
                # tiles just read above.
                cuda.syncthreads()

            if valid:
                idx = i * n_samples - (i * (i + 1)) // 2 + (j - i - 1)
                if method == WEIGHTED_UNNORMALIZED or method == UNWEIGHTED_UNNORMALIZED:
                    out[idx] = numerator
                else:
                    out[idx] = 0.0 if denominator == 0.0 else numerator / denominator

        # The module reference keeps its id() from being reused by another
        # object while this entry lives.
        _KERNEL_CACHE[cache_key] = (cuda, _unifrac_block_kernel)
        return _unifrac_block_kernel


def _launch_unifrac_kernel(
    cuda,
    tile,
    node_chunk,
    proportions,
    counts_by_node,
    sample_totals,
    branch_lengths,
    method,
    alpha,
    variance_adjust,
):
    """Upload the inputs, launch the kernel, and return the condensed vector.

    Every array must already be contiguous float64. Shared by
    ``_run_unifrac_kernel`` and by the tile-config probe below, so the probe
    exercises the launch path production uses.
    """
    n_samples = proportions.shape[0]
    n_pairs = n_samples * (n_samples - 1) // 2
    block_i, block_j = _build_block_pair_index(n_samples, tile)

    # Bound to locals, not inlined into the launch below: the kernel launch is
    # asynchronous, so each device array must stay referenced until the
    # copy_to_host() that synchronizes on it.
    d_proportions = cuda.to_device(proportions)
    # Callers that need no division pass the same array twice; upload once.
    if counts_by_node is proportions:
        d_counts = d_proportions
    else:
        d_counts = cuda.to_device(counts_by_node)
    d_sample_totals = cuda.to_device(sample_totals)
    d_branch_lengths = cuda.to_device(branch_lengths)
    d_block_i = cuda.to_device(block_i)
    d_block_j = cuda.to_device(block_j)
    d_out = cuda.device_array(n_pairs, dtype=np.float64)

    # Only the GENERALIZED kernel branch reads `alpha_mode`.
    if alpha == 1.0:
        alpha_mode = 1
    elif alpha == 0.5:
        alpha_mode = 2
    elif alpha == 0.0:
        alpha_mode = 3
    else:
        alpha_mode = 0

    kernel = _make_unifrac_kernel(cuda, tile, node_chunk, variance_adjust)
    try:
        kernel[block_i.shape[0], (tile, tile)](
            d_proportions,
            d_counts,
            d_sample_totals,
            d_branch_lengths,
            method,
            alpha,
            alpha_mode,
            d_block_i,
            d_block_j,
            n_samples,
            d_out,
        )
        return d_out.copy_to_host()
    except Exception:
        # A kernel is cached as soon as it compiles, but a config can compile
        # and still fail at launch (NVIDIA's LAUNCH_OUT_OF_RESOURCES), which
        # would leave an unusable entry behind for the rest of the process.
        _evict_kernel_cache_entry(cuda, tile, node_chunk, variance_adjust)
        raise


def _evict_kernel_cache_entry(cuda, tile, node_chunk, variance_adjust):
    """Remove ``(id(cuda), tile, node_chunk, variance_adjust)`` from
    ``_KERNEL_CACHE``, if present.
    """
    cache_key = (id(cuda), tile, node_chunk, variance_adjust)
    with _KERNEL_CACHE_LOCK:
        _KERNEL_CACHE.pop(cache_key, None)


# (TILE, NODE_CHUNK) candidates for the kernel's block tiling, most
# aggressive first. Which one is safe depends on the device's register file
# and shared-memory budget, which varies across SKUs of the same vendor, so
# it cannot be read off a per-backend table and is probed at runtime instead
# (see `_probe_tile_config`). Measured notes:
#   (32, 32): fastest on AMD MI300A (gfx942), but fails at launch with
#       CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES on every NVIDIA card tested
#       (Turing, Ampere) -- 1024 threads/block needs more registers/thread
#       under NVVM than NVIDIA's per-SM budget allows.
#   (16, 32): the rung known good on those NVIDIA cards.
#   (8, 32): conservative floor for tighter per-block limits; untested.
_TILE_CANDIDATES = [
    (32, 32),
    (16, 64),
    (16, 32),
    (8, 32),
]

# Discovered-safe (TILE, NODE_CHUNK), keyed by backend name ('cuda'/'hip')
# rather than by device: a process talks to one vendor at a time in practice,
# and probing costs real compiles and launches.
_TILE_CONFIG_CACHE = {}

# Guards the check-probe-store sequence in `_get_tile_config`.
_TILE_CONFIG_LOCK = threading.Lock()


def _probe_tile_config_candidate(cuda, tile, node_chunk, variance_adjust):
    """Compile and launch the real kernel at ``(tile, node_chunk)`` against
    tiny synthetic dummy data, to see whether this configuration is safe on
    the GPU actually present.

    Goes through `_launch_unifrac_kernel` rather than a simplified stand-in,
    because the failures this guards against depend on the real kernel's
    register and shared-memory usage: an NVIDIA launch-time resource error,
    or an AMD code-generation error (numba.hip compiles lazily, on first
    invocation). The launch ends in a synchronizing ``copy_to_host()``, so a
    launch-time error surfaces here rather than later.

    Raises whatever exception the backend produces; the caller decides
    whether to try the next candidate.

    ``variance_adjust`` selects which compile-time kernel variant to probe.
    Register allocation is decided by NVVM/LLVM codegen over the whole
    function body and is not monotonic in the amount of source removed: on an
    NVIDIA A10 the ``variance_adjust=False`` variant hit
    ``CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES`` at a config where the
    ``variance_adjust=True`` variant launched fine. Each variant therefore
    needs its own probe.
    """
    # Two samples (one i<j pair) and two nodes are enough to exercise one
    # full block launch; the limits probed for do not depend on data size.
    n_samples = 2
    n_nodes = 2
    _launch_unifrac_kernel(
        cuda,
        tile,
        node_chunk,
        np.zeros((n_samples, n_nodes), dtype=np.float64),
        np.zeros((n_samples, n_nodes), dtype=np.float64),
        np.zeros(n_samples, dtype=np.float64),
        np.zeros(n_nodes, dtype=np.float64),
        UNWEIGHTED,
        1.0,
        variance_adjust,
    )


def _probe_tile_config(cuda, backend):
    """Find the first candidate in `_TILE_CANDIDATES` that works on the GPU
    behind ``cuda``, trying each in order (most aggressive first).

    Catches any exception, not just the two known failure shapes: this is a
    controlled probe against synthetic data, not user input.

    A candidate is accepted only if both kernel variants probe clean at it,
    since `_get_tile_config` caches one config per backend and hands it to
    whichever variant a real call needs.
    """
    failures = []
    for tile, node_chunk in _TILE_CANDIDATES:
        try:
            _probe_tile_config_candidate(cuda, tile, node_chunk, True)
            _probe_tile_config_candidate(cuda, tile, node_chunk, False)
            return tile, node_chunk
        except Exception as exc:
            failures.append(
                f"(TILE={tile}, NODE_CHUNK={node_chunk}): {type(exc).__name__}: {exc}"
            )

    detail = "\n  ".join(failures)
    raise RuntimeError(
        f"No safe (TILE, NODE_CHUNK) tile configuration could be found for "
        f"the '{backend}' GPU backend on this device; every candidate in "
        f"_TILE_CANDIDATES failed to compile or launch:\n  {detail}"
    )


def _get_tile_config(cuda, backend):
    """Return a (TILE, NODE_CHUNK) pair known to work on this device for
    ``backend``, probing and caching it on first use.
    """
    with _TILE_CONFIG_LOCK:
        cached = _TILE_CONFIG_CACHE.get(backend)
        if cached is None:
            cached = _probe_tile_config(cuda, backend)
            _TILE_CONFIG_CACHE[backend] = cached
        return cached


class _ValidationFailure(Exception):
    """Internal marker: wraps whatever ``_setup_multiple_unifrac`` raised for
    bad ``counts``/``taxa``/``tree`` input, before any GPU-specific code ran.

    ``_dispatch_gpu_or_xp`` catches this to re-raise the original exception
    without marking the GPU backend unavailable. The signal is *position*
    (raised only around the validation call, which never touches the GPU),
    not the exception's type: bad input can surface as ``ValueError``,
    ``TreeError``, or a plain ``KeyError``, and ``KeyError`` is also a
    plausible failure from numba/CUDA-driver internals.
    """

    def __init__(self, original):
        self.original = original
        super().__init__(str(original))


def _run_unifrac_kernel(
    counts, taxa, tree, method, *, alpha=1.0, variance_adjust=False, validate=True
):
    """Launch the fused kernel and return the condensed distance vector.

    Shared host-side body of the three ``*_unifrac_gpu`` drivers below: the
    five methods differ only in the ``method``/``alpha`` the one kernel is
    launched with, and in whether it reads node proportions or raw counts.
    Requires a GPU backend detected by ``_detect_gpu_backend``.
    """
    from skbio.diversity.beta._unifrac import _setup_multiple_unifrac, _get_tip_indices

    cuda = _get_cuda_module()
    # Validate before tile probing, which does a real compile and launch: a
    # bad `taxa`/`tree` should surface its own error rather than a probe
    # failure masking it. Wrapped so `_dispatch_gpu_or_xp` can tell input
    # errors from kernel errors by position (see `_ValidationFailure`).
    try:
        counts_by_node, tree_index, branch_lengths = _setup_multiple_unifrac(
            counts, taxa, tree, validate
        )
    except Exception as exc:
        raise _ValidationFailure(exc) from exc
    tile, node_chunk = _get_tile_config(cuda, _detect_gpu_backend())
    counts_by_node = np.ascontiguousarray(counts_by_node, dtype=np.float64)

    # Unweighted without variance adjustment reads sample_totals nowhere:
    # neither the kernel's `variance_adjust` branch nor the division below.
    needs_sample_totals = variance_adjust or method not in (
        UNWEIGHTED,
        UNWEIGHTED_UNNORMALIZED,
    )
    if needs_sample_totals:
        tip_indices = _get_tip_indices(tree_index)
        # A product against a 0/1 tip mask avoids materializing the
        # fancy-indexed temporary `counts_by_node[:, tip_indices]`.
        tip_mask = np.zeros(counts_by_node.shape[1], dtype=np.float64)
        tip_mask[tip_indices] = 1.0
        sample_totals = counts_by_node @ tip_mask
    else:
        sample_totals = np.zeros(counts_by_node.shape[0], dtype=np.float64)

    if method in (UNWEIGHTED, UNWEIGHTED_UNNORMALIZED):
        # The unweighted kernel branch only tests "proportions" for > 0, so
        # raw counts can stand in for them unnormalized.
        proportions = counts_by_node
    else:
        # Counts are non-negative, so sample_totals[i] == 0 implies that whole
        # row is zero; substituting 1.0 for the zero divisors leaves those
        # rows at zero.
        proportions = counts_by_node / np.where(
            sample_totals > 0, sample_totals, 1.0
        )[:, None]
        if not variance_adjust:
            # The kernel reads `counts` only when variance_adjust is True.
            # Aliasing lets `_launch_unifrac_kernel` skip the second upload.
            counts_by_node = proportions

    return _launch_unifrac_kernel(
        cuda,
        tile,
        node_chunk,
        proportions,
        counts_by_node,
        sample_totals,
        np.ascontiguousarray(branch_lengths, dtype=np.float64),
        method,
        alpha,
        variance_adjust,
    )


def _weighted_unifrac_gpu(
    counts, taxa, tree, normalized, variance_adjust=False, validate=True
):
    """Compute the condensed weighted UniFrac distance vector on a GPU.

    Host-side driver for the weighted UniFrac methods (``normalized`` selects
    between ``WEIGHTED_NORMALIZED`` and ``WEIGHTED_UNNORMALIZED``); see
    ``_run_unifrac_kernel``.

    """
    return _run_unifrac_kernel(
        counts,
        taxa,
        tree,
        WEIGHTED_NORMALIZED if normalized else WEIGHTED_UNNORMALIZED,
        variance_adjust=variance_adjust,
        validate=validate,
    )


def _unweighted_unifrac_gpu(
    counts, taxa, tree, normalized, variance_adjust=False, validate=True
):
    """Compute the condensed unweighted UniFrac distance vector on a GPU.

    Host-side driver for the unweighted UniFrac methods (``normalized``
    selects between ``UNWEIGHTED`` and ``UNWEIGHTED_UNNORMALIZED``); see
    ``_run_unifrac_kernel``.

    """
    return _run_unifrac_kernel(
        counts,
        taxa,
        tree,
        UNWEIGHTED if normalized else UNWEIGHTED_UNNORMALIZED,
        variance_adjust=variance_adjust,
        validate=validate,
    )


def _generalized_unifrac_gpu(
    counts, taxa, tree, alpha, variance_adjust=False, validate=True
):
    """Compute the condensed generalized UniFrac distance vector on a GPU.

    Host-side driver for the generalized UniFrac method (``GENERALIZED``);
    see ``_run_unifrac_kernel``.

    """
    return _run_unifrac_kernel(
        counts,
        taxa,
        tree,
        GENERALIZED,
        alpha=alpha,
        variance_adjust=variance_adjust,
        validate=validate,
    )


# Backends ('cuda'/'hip') whose fused kernel failed to build or run in this
# process. Populated by `_mark_backend_unavailable` so later calls skip
# straight to the array-API path instead of retrying a failing kernel.
_unavailable_backends = set()

# Guards the check-then-add in `_mark_backend_unavailable` and the
# check-then-attempt in `_dispatch_gpu_or_xp`.
_UNAVAILABLE_BACKENDS_LOCK = threading.Lock()


def _mark_backend_unavailable(backend):
    """Record that ``backend``'s fused UniFrac kernel cannot run this process.

    Called when the fused kernel raises after a GPU backend was detected
    (for example a numba-hip build that fails to compile on the running ROCm
    stack). Warns once per backend, then routes that backend to the
    array-API fallback from then on.
    """
    with _UNAVAILABLE_BACKENDS_LOCK:
        if backend not in _unavailable_backends:
            _unavailable_backends.add(backend)
            should_warn = True
        else:
            should_warn = False
    if should_warn:
        warn(
            f"The fused UniFrac GPU kernel could not be used for the "
            f"'{backend}' backend on this system; using the array-API "
            "fallback instead.",
            UserWarning,
        )


def _dispatch_gpu_or_xp(gpu_func, xp_func, *args, **kwargs):
    """Call ``gpu_func`` if the fused kernel is usable, else ``xp_func``.

    An exception out of the fused kernel marks its backend unavailable for
    the rest of the process and falls through to the array-API path, so one
    bad kernel build degrades performance rather than failing the call.

    Input-validation failures are excluded from that: they are the caller's
    fault, not the backend's, so re-raising them keeps one bad input from
    permanently demoting the process to the slow path. They arrive either
    wrapped in ``_ValidationFailure`` (from the real ``*_unifrac_gpu``
    drivers) or as a bare ``ValueError``/``skbio.tree.TreeError`` from any
    other ``gpu_func``.
    """
    from skbio.tree import TreeError

    backend = _detect_gpu_backend()
    with _UNAVAILABLE_BACKENDS_LOCK:
        backend_usable = backend is not None and backend not in _unavailable_backends
    if backend_usable:
        try:
            return gpu_func(*args, **kwargs)
        except _ValidationFailure as exc:
            raise exc.original from None
        except (ValueError, TreeError):
            raise
        except Exception:
            _mark_backend_unavailable(backend)
    return xp_func(*args, **kwargs)


def _weighted_unifrac_gpu_or_xp(
    counts, taxa, tree, normalized, variance_adjust=False, validate=True
):
    """Compute weighted UniFrac, preferring the fused GPU kernel when usable.

    Runs :func:`_weighted_unifrac_gpu` when a GPU backend is detected and its
    fused kernel builds/runs successfully; otherwise falls back to
    :func:`skbio.diversity.beta._unifrac_xp._weighted_unifrac_xp`, which is
    correct (if slower) on any array-API backend, including plain NumPy.

    """
    from skbio.diversity.beta._unifrac_xp import _weighted_unifrac_xp

    return _dispatch_gpu_or_xp(
        _weighted_unifrac_gpu,
        _weighted_unifrac_xp,
        counts,
        taxa,
        tree,
        normalized,
        variance_adjust=variance_adjust,
        validate=validate,
    )


def _unweighted_unifrac_gpu_or_xp(
    counts, taxa, tree, normalized, variance_adjust=False, validate=True
):
    """Compute unweighted UniFrac, preferring the fused GPU kernel when usable.

    See :func:`_weighted_unifrac_gpu_or_xp`; falls back to
    :func:`skbio.diversity.beta._unifrac_xp._unweighted_unifrac_xp`.

    """
    from skbio.diversity.beta._unifrac_xp import _unweighted_unifrac_xp

    return _dispatch_gpu_or_xp(
        _unweighted_unifrac_gpu,
        _unweighted_unifrac_xp,
        counts,
        taxa,
        tree,
        normalized,
        variance_adjust=variance_adjust,
        validate=validate,
    )


def _generalized_unifrac_gpu_or_xp(
    counts, taxa, tree, alpha, variance_adjust=False, validate=True
):
    """Compute generalized UniFrac, preferring the fused GPU kernel when usable.

    See :func:`_weighted_unifrac_gpu_or_xp`; falls back to
    :func:`skbio.diversity.beta._unifrac_xp._generalized_unifrac_xp`.

    """
    from skbio.diversity.beta._unifrac_xp import _generalized_unifrac_xp

    return _dispatch_gpu_or_xp(
        _generalized_unifrac_gpu,
        _generalized_unifrac_xp,
        counts,
        taxa,
        tree,
        alpha,
        variance_adjust=variance_adjust,
        validate=validate,
    )
