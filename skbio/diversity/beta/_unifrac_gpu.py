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

# Compiled kernels, keyed by id() of the cuda-API module they were built
# against. The module itself is stored alongside so the id stays valid.
#
# (`skbio.stats.distance._gpu._get_kernel` is the PERMANOVA/Mantel precedent
# for a "compile once per backend, cache it" helper, but its cache key is
# backend name alone; UniFrac's cache key needs an extra dimension -- `tile`,
# `node_chunk` -- since the tile-probe below tries several kernel builds per
# backend, and the probe also needs to evict one specific failed entry
# (`_evict_kernel_cache_entry`), which `_get_kernel`'s simpler contract has no
# equivalent of. So a small UniFrac-local cache is used here instead of
# reusing or extending that helper.)
_KERNEL_CACHE = {}

# Guards the check-compile-store sequence in `_make_unifrac_kernel` so
# concurrent calls from different threads cannot both miss the cache and
# redundantly recompile. Not a hot path (compilation is memoized after the
# first call), so a single coarse lock is sufficient.
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
    installed ROCm version (for example
    ``numba-hip[rocm-7-0-0]==0.1.6`` for ROCm 7.0.0), so there is no
    single pip specification that works across ROCm releases and hence no
    ``gpu-amd`` extra. Note that ``hip-python`` alone is not sufficient: it
    provides low-level HIP bindings, not ``numba.hip``/``pose_as_cuda()``.
    See ``docs/superpowers/specs/2026-10-02-ssu-unifrac-numba-phase1-design.md``
    for the full install path.

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


def _make_unifrac_kernel(cuda, tile, node_chunk):
    """Return the UniFrac pair kernel compiled against a cuda-API module.

    Shared by all 5 UniFrac methods (``UNWEIGHTED``, ``UNWEIGHTED_UNNORMALIZED``,
    ``WEIGHTED_NORMALIZED``, ``WEIGHTED_UNNORMALIZED``, ``GENERALIZED``) so that
    only one kernel is compiled and launched regardless of which method a given
    driver function dispatches.

    ``tile``/``node_chunk`` are closed over as compile-time constants, since
    the shared-memory tiles are sized from them; which pair is safe differs
    per device, so it is discovered at runtime (see ``_get_tile_config``).

    ``cuda.jit`` returns a fresh Dispatcher (and so pays full compilation
    cost, seconds) each time this runs, which would otherwise happen on every
    single ``*_unifrac_gpu`` call. The result is therefore memoized per
    backend module (and tile config) in ``_KERNEL_CACHE``; only the compiled
    kernel object is cached, never any per-call state.

    """
    cache_key = (id(cuda), tile, node_chunk)
    with _KERNEL_CACHE_LOCK:
        cached = _KERNEL_CACHE.get(cache_key)
        if cached is not None:
            return cached[1]

        # Imported here, not at module level: this module must stay
        # importable (and the array-API fallback reachable) on a system with
        # no numba installed at all. This function is only ever reached via
        # `_launch_unifrac_kernel`, after `_get_cuda_module()` has already
        # succeeded, which means some numba variant (numba-cuda or
        # numba.hip, both of which depend on core numba) is installed.
        from numba import float64

        TILE = tile
        NODE_CHUNK = node_chunk
        # A plain Python int, not an expression evaluated inside the jitted
        # function: `numba.hip`'s `cuda.shared.array` typing requires every
        # shape element to already be a literal int and cannot constant-fold
        # `NODE_CHUNK + 1` itself (confirmed on real AMD hardware -- that
        # form types as a non-literal int64 and fails with a TypingError on
        # the hip backend, even though the identical source compiles fine
        # under numba-cuda). Precomputing it here as its own closure
        # variable keeps both backends happy.
        NODE_CHUNK_PADDED = node_chunk + 1

        @cuda.jit
        def _unifrac_block_kernel(
            proportions,
            counts,
            sample_totals,
            branch_lengths,
            method,
            alpha,
            variance_adjust,
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

            # Offset within the i-block / j-block respectively.
            tx = cuda.threadIdx.x
            ty = cuda.threadIdx.y

            i = bi * TILE + tx
            j = bj * TILE + ty
            valid = i < n_samples and j < n_samples and i < j

            n_nodes = proportions.shape[1]

            # sh_prop_i/sh_cnt_i are indexed [tx, col] in the inner loop below,
            # with tx (a lane index) fast-varying and col loop-invariant
            # within a warp/wavefront -- a column stride of exactly
            # NODE_CHUNK, which is always a power of two in
            # `_TILE_CANDIDATES`, causing shared-memory bank conflicts
            # (conflict-free 8-byte access wants an odd element stride). Pad
            # by one column to break that. sh_prop_j/sh_cnt_j are read
            # [ty, col], a warp-uniform access (every lane in a
            # warp/wavefront reads the same address), which has no bank
            # conflict to fix, so they are left unpadded.
            sh_prop_i = cuda.shared.array((TILE, NODE_CHUNK_PADDED), dtype=float64)
            sh_prop_j = cuda.shared.array((TILE, NODE_CHUNK), dtype=float64)
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

                # Cooperative load: every thread in the block (valid or not)
                # helps fill the i-tile and j-tile shared buffers for this
                # node-chunk, strided by the number of threads in the block.
                for e in range(lin_tid, tile_elems, n_threads_blk):
                    row = e // NODE_CHUNK
                    col = e % NODE_CHUNK
                    node = chunk_start + col
                    samp_i = bi * TILE + row
                    samp_j = bj * TILE + row
                    if col < chunk_len and samp_i < n_samples:
                        sh_prop_i[row, col] = proportions[samp_i, node]
                        sh_cnt_i[row, col] = counts[samp_i, node]
                    else:
                        sh_prop_i[row, col] = 0.0
                        sh_cnt_i[row, col] = 0.0
                    if col < chunk_len and samp_j < n_samples:
                        sh_prop_j[row, col] = proportions[samp_j, node]
                        sh_cnt_j[row, col] = counts[samp_j, node]
                    else:
                        sh_prop_j[row, col] = 0.0
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
                        if variance_adjust:
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
                                if variance_adjust:
                                    numerator += length / vaw if differs else 0.0
                                    denominator += length / vaw
                                else:
                                    numerator += length if differs else 0.0
                                    denominator += length
                        elif method == GENERALIZED:
                            if s == 0.0:
                                continue
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

        # Keep a reference to the module so its id() cannot be reused by another
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
    exercises exactly the launch path production uses and cannot drift from
    it.
    """
    n_samples = proportions.shape[0]
    n_pairs = n_samples * (n_samples - 1) // 2
    block_i, block_j = _build_block_pair_index(n_samples, tile)

    # Bound to locals, not inlined into the launch below: the kernel launch is
    # asynchronous, so each device array must stay referenced until the
    # copy_to_host() that synchronizes on it.
    d_proportions = cuda.to_device(proportions)
    # For the unweighted methods, `proportions` *is* `counts_by_node` (the
    # caller needs no division), so reuse the one transfer already made for
    # `d_proportions` instead of uploading the identical host array twice.
    if counts_by_node is proportions:
        d_counts = d_proportions
    else:
        d_counts = cuda.to_device(counts_by_node)
    d_sample_totals = cuda.to_device(sample_totals)
    d_branch_lengths = cuda.to_device(branch_lengths)
    d_block_i = cuda.to_device(block_i)
    d_block_j = cuda.to_device(block_j)
    d_out = cuda.device_array(n_pairs, dtype=np.float64)

    kernel = _make_unifrac_kernel(cuda, tile, node_chunk)
    try:
        kernel[block_i.shape[0], (tile, tile)](
            d_proportions,
            d_counts,
            d_sample_totals,
            d_branch_lengths,
            method,
            alpha,
            variance_adjust,
            d_block_i,
            d_block_j,
            n_samples,
            d_out,
        )
        return d_out.copy_to_host()
    except Exception:
        # `_make_unifrac_kernel` caches the compiled kernel before it is
        # ever launched, so a candidate that compiles fine but fails here
        # (the NVIDIA tile-probe failure mode: launch-time
        # LAUNCH_OUT_OF_RESOURCES) would otherwise leave a cached-but-
        # unusable entry -- holding a reference to the CUDA module and a
        # broken Dispatcher -- for the rest of the process. Evict it so a
        # failed candidate does not linger; nothing else looks up this
        # specific (tile, node_chunk) key again once probing has moved on
        # (`_get_tile_config` only returns the config that actually worked).
        _evict_kernel_cache_entry(cuda, tile, node_chunk)
        raise


def _evict_kernel_cache_entry(cuda, tile, node_chunk):
    """Remove ``(id(cuda), tile, node_chunk)`` from ``_KERNEL_CACHE``, if present.

    See the ``except`` clause in ``_launch_unifrac_kernel`` for why this is
    needed: a kernel is cached as soon as it compiles, before it is ever
    launched.
    """
    cache_key = (id(cuda), tile, node_chunk)
    with _KERNEL_CACHE_LOCK:
        _KERNEL_CACHE.pop(cache_key, None)


# Block-tiling parameters for the 2D node-chunked kernel (see
# `_make_unifrac_kernel`). A block is (TILE, TILE) threads, i.e. TILE**2
# threads/block, using 2 shared (TILE, NODE_CHUNK + 1) float64 tiles
# (`sh_prop_i`/`sh_cnt_i`, padded by one column to avoid the shared-memory
# bank conflict their `[tx, col]` access pattern would otherwise cause --
# see the comment above their declaration in `_make_unifrac_kernel`), 2
# unpadded shared (TILE, NODE_CHUNK) tiles (`sh_prop_j`/`sh_cnt_j`, read
# warp-uniformly so padding buys nothing there), plus one (NODE_CHUNK,)
# tile. The safe choice is hardware-specific, not just backend-specific: it
# depends on the GPU's actual register file and shared-memory budget, which
# varies across SKUs of the *same* vendor, not only between vendors. There
# is no reliable way to pick it from a static per-backend table, so it is
# discovered at runtime instead (see `_get_tile_config`/`_probe_tile_config`
# below) by actually compiling and launching the real kernel against
# trivial dummy data.
#
# `_TILE_CANDIDATES`, most aggressive (most shared memory, most
# threads/block, fastest when it works) to least, are real data points from
# tuning this kernel on specific hardware, not invented, and double as the
# starting ladder for the runtime probe. Shared-memory totals below are
# post-padding (2*TILE*(NODE_CHUNK+1)*8 + 2*TILE*NODE_CHUNK*8 bytes for the
# 4 node tiles, plus NODE_CHUNK*8 for the branch-length tile):
#   (32, 32): 1024 threads/block, 2*32*33*8 + 2*32*32*8 = 33280 bytes
#       (~32.5KB) shared mem. Fastest config measured, on AMD MI300A
#       (gfx942) at n_samples=5000 (~12.5M pairs), median 0.952s (measured
#       before the A1 padding fix below; see that commit for the
#       post-padding re-measurement). Crashes at kernel *launch* on every
#       NVIDIA GPU tested (RTX 2080 Ti/Turing, A10/Ampere) with
#       CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES: this kernel's numba-cuda/NVVM
#       compilation uses more registers/thread at 1024 threads/block than
#       fit NVIDIA's per-SM register budget, even though the identical
#       source is fine on AMD's CU register file.
#   (16, 64): also 1024 threads/block but a narrower, taller shared tile
#       (2*16*65*8 + 2*16*64*8 = 33024 bytes, ~32.25KB). Measured on AMD
#       MI300A: works, but slower than (32, 32) (median 1.325s).
#   (16, 32): 256 threads/block, 2*16*33*8 + 2*16*32*8 = 16640 bytes
#       (~16.25KB) shared mem. The config known to launch and run correctly
#       on both NVIDIA cards tested (RTX 2080 Ti, A10); previously the
#       hardcoded NVIDIA default.
#   (8, 32): 64 threads/block, 2*8*33*8 + 2*8*32*8 = 8320 bytes (~8.125KB)
#       shared mem. A defensive floor rung below anything actually tuned,
#       for hardware with even tighter per-block limits than tested
#       (untested itself, but conservative enough that it should be broadly
#       safe).
# A real candidate that failed to compile during tuning: (16, 128) on AMD
# -- "local memory (66560) exceeds limit (65536)" (pre-padding: 4 shared
# tiles of 16*128*8 = 16KB each = 64KB, right at AMD's 64KB LDS/block limit
# and pushed over by other locals). Not included in the ladder since (32, 32)
# already dominates it on the one backend where it even compiles.
_TILE_CANDIDATES = [
    (32, 32),
    (16, 64),
    (16, 32),
    (8, 32),
]

# Discovered-safe (TILE, NODE_CHUNK) per backend, filled in by
# `_get_tile_config` the first time it is called for a given backend.
# Probing does real compiles/launches (not free), so the result is cached
# for the rest of the process rather than re-probed on every call. Keyed by
# backend name ('cuda'/'hip'), not by individual device: a process talks to
# one GPU backend/vendor at a time in practice, and re-probing per distinct
# device name would pay the probe cost again for every differently-named
# card with no real precedent (within this session's tuning) that it would
# ever choose a different answer on same-vendor hardware.
_TILE_CONFIG_CACHE = {}

# Guards the check-probe-store sequence in `_get_tile_config`, mirroring
# `_KERNEL_CACHE_LOCK` above for the same reason: probing is not a hot path
# (it runs at most once per backend per process), so a single coarse lock
# is sufficient to stop concurrent callers from redundantly re-probing.
_TILE_CONFIG_LOCK = threading.Lock()


def _probe_tile_config_candidate(cuda, tile, node_chunk):
    """Compile and launch the real kernel at ``(tile, node_chunk)`` against
    tiny synthetic dummy data, to see whether this configuration is safe on
    the GPU actually present.

    Deliberately goes through `_launch_unifrac_kernel`, i.e. a real block
    launch of the *actual* kernel rather than a simplified stand-in, since
    the failure modes this guards against are specific to that kernel's real
    register/shared-memory usage: an NVIDIA launch-time resource error, or
    an AMD compile-time (numba.hip compiles lazily, on first invocation)
    code-generation error. The launch ends in a ``copy_to_host()``, which
    synchronizes, so a launch-time error surfaces here rather than silently
    later.

    Raises whatever exception the backend produces; the caller decides
    whether to try the next candidate.
    """
    # Minimal data: 2 samples (one valid i<j pair) and 2 nodes (one real
    # node-chunk iteration), just enough to exercise one full (tile, tile)
    # block launch of the real kernel body -- cheap to construct, and the
    # resource limits this probes for (threads/block, shared-memory/block)
    # depend only on the launch configuration and compiled kernel, not on
    # the data size.
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
        False,
    )


def _probe_tile_config(cuda, backend):
    """Find the first candidate in `_TILE_CANDIDATES` that actually works
    on the GPU behind ``cuda``, trying each in order (most aggressive
    first) and catching both known failure shapes (NVIDIA launch-time,
    AMD compile-time) plus anything else, since this is a controlled probe
    against synthetic data, not user input -- contrast `_dispatch_gpu_or_xp`,
    whose ``except`` narrows deliberately because its inputs *are*
    user-controlled.
    """
    failures = []
    for tile, node_chunk in _TILE_CANDIDATES:
        try:
            _probe_tile_config_candidate(cuda, tile, node_chunk)
            return tile, node_chunk
        except Exception as exc:  # noqa: BLE001 -- see docstring above
            failures.append(
                f"(TILE={tile}, NODE_CHUNK={node_chunk}): {type(exc).__name__}: {exc}"
            )

    detail = "\n  ".join(failures)
    # If you land here debugging a real failure: every rung in
    # _TILE_CANDIDATES above failed on this device, which the ladder was
    # built to avoid (its smallest rung, (8, 32), is meant to be
    # conservative enough to not need this). What to do depends on the
    # error text in `detail` for the *smallest* candidate, (8, 32):
    #   - "LAUNCH_OUT_OF_RESOURCES" / a register-count complaint (NVIDIA
    #     shape): this device's register file is smaller than any tested
    #     so far. Add a new rung below (8, 32) -- e.g. (8, 16) or (4, 32)
    #     -- to _TILE_CANDIDATES above, in the same most-to-least-aggressive
    #     order, with a comment recording the device and the real error.
    #   - "local memory ... exceeds limit" / a shared-memory complaint (AMD
    #     shape): this device's LDS/shared-memory-per-block budget is
    #     smaller than 8*32*8*4 = 8KB. Same fix: add a smaller rung, sized
    #     under that device's real limit (check `detail` for the exact
    #     numbers the compiler reported).
    #   - Anything else (an import error, a missing symbol, a totally
    #     different exception shape): this probably is not a tile-size
    #     problem at all -- it's more likely a toolchain/driver issue on
    #     this specific machine. Don't just shrink the ladder in that case;
    #     investigate why `_get_cuda_module()` returned a module that can't
    #     actually compile this kernel.
    raise RuntimeError(
        f"No safe (TILE, NODE_CHUNK) tile configuration could be found for "
        f"the '{backend}' GPU backend on this device; every candidate in "
        f"_TILE_CANDIDATES failed to compile or launch:\n  {detail}"
    )


def _get_tile_config(cuda, backend):
    """Return a (TILE, NODE_CHUNK) pair known to work on this device for
    ``backend``, probing and caching it on first use.

    Not free the first time (it compiles and launches the real kernel,
    possibly more than once), so the result is memoized in
    `_TILE_CONFIG_CACHE` for the rest of the process; later calls for the
    same backend hit the cache and never re-probe. See `_probe_tile_config`
    for the probing strategy.
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
    immediately, without marking the GPU backend unavailable. The signal
    this relies on is *position* (raised only around the validation call,
    which never touches the GPU), not the original exception's type -- a
    type-based allowlist is not safe here, since the real exception for bad
    input varies (``ValueError``/``TreeError`` under ``validate=True``,
    plain ``KeyError`` from ``_nodes_by_counts`` under ``validate=False``
    for mismatched taxa/tree) and at least ``KeyError`` is also a plausible
    failure mode from unrelated numba/CUDA-driver internals during genuine
    kernel compilation or launch.
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
    # Validate before any GPU-specific code, including tile-config probing:
    # probing now does a real compile+launch against this device (unlike
    # the old fixed-dict lookup), so it must not run ahead of input
    # validation -- a bad `taxa`/`tree` should surface its own error
    # immediately, not after paying for a wasted probe (or worse, a probe
    # failure masking the real validation error). `_setup_multiple_unifrac`
    # never touches the GPU, so anything it raises here -- a `ValueError`/
    # `TreeError` under `validate=True`, or a plain `KeyError` from
    # `_nodes_by_counts` for mismatched taxa under `validate=False` -- is by
    # construction an input problem, not a kernel/hardware one. Wrap it so
    # `_dispatch_gpu_or_xp` can tell the two apart by *where* the exception
    # came from rather than by its type (see `_ValidationFailure` and
    # `_dispatch_gpu_or_xp`'s docstring: exception type alone is not a safe
    # signal, since e.g. `KeyError` could also plausibly come from the real
    # GPU kernel/launch machinery for an unrelated reason).
    try:
        counts_by_node, tree_index, branch_lengths = _setup_multiple_unifrac(
            counts, taxa, tree, validate
        )
    except Exception as exc:
        raise _ValidationFailure(exc) from exc
    tile, node_chunk = _get_tile_config(cuda, _detect_gpu_backend())
    counts_by_node = np.ascontiguousarray(counts_by_node, dtype=np.float64)
    tip_indices = _get_tip_indices(tree_index)
    sample_totals = counts_by_node[:, tip_indices].sum(axis=1)
    if method in (UNWEIGHTED, UNWEIGHTED_UNNORMALIZED):
        # The unweighted kernel branch only tests "proportions" for > 0
        # (presence/absence), so raw counts can be passed directly in place of
        # proportions, no division needed.
        proportions = counts_by_node
    else:
        proportions = np.divide(
            counts_by_node,
            sample_totals[:, None],
            out=np.zeros_like(counts_by_node),
            where=sample_totals[:, None] > 0,
        )
    return _launch_unifrac_kernel(
        cuda,
        tile,
        node_chunk,
        proportions,
        counts_by_node,
        sample_totals,
        branch_lengths.astype(np.float64),
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


# -----------------------------------------------------------------------------
# Two-tier dispatch: fused GPU kernel when usable, array-API fallback otherwise
# -----------------------------------------------------------------------------
#
# The three driver functions above require a real, detected GPU backend
# (``_get_cuda_module`` raises ``ImportError`` otherwise): they already run
# correctly on real NVIDIA/AMD hardware, including proactively moving plain
# NumPy input onto the device for the fused-kernel speedup. The wrappers below
# add the fallback: when no GPU backend is detected, or when the fused kernel
# fails to build/run on the running stack, they run the array-API-generic
# implementation in :mod:`skbio.diversity.beta._unifrac_xp` instead, which is
# correct (if slower) on any array-API backend, including plain NumPy on CPU.

# Backends ('cuda'/'hip') whose fused kernel failed to build or run in this
# process. Populated by `_mark_backend_unavailable` so later calls skip
# straight to the array-API path instead of retrying a failing kernel.
# (`skbio.stats.distance._gpu._mark_gpu_unavailable` is the PERMANOVA/Mantel
# precedent for this, but it keys off an array-API array's device/backend;
# UniFrac's GPU dispatch is keyed off `_detect_gpu_backend()`'s backend name
# instead -- plain NumPy input has no device of its own -- so a small
# UniFrac-local equivalent is used here rather than reusing that helper.)
_unavailable_backends = set()

# Guards the check-then-add in `_mark_backend_unavailable` and the
# check-then-attempt in `_dispatch_gpu_or_xp`, mirroring `_KERNEL_CACHE_LOCK`/
# `_TILE_CONFIG_LOCK` above for the same reason: without it, two concurrent
# callers could both pass the `not in` check and both attempt the failing
# GPU path (or both emit the "kernel could not be used" warning).
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
    bad kernel build degrades performance rather than failing the call --
    *except* for input-validation failures. There are two ways ``gpu_func``
    surfaces one of those:

    - Wrapped in ``_ValidationFailure``, when ``gpu_func`` is one of the real
      ``*_unifrac_gpu`` drivers: ``_run_unifrac_kernel`` wraps whatever its
      own ``_setup_multiple_unifrac`` call raises before any GPU-specific
      code runs. See ``_ValidationFailure`` for why this is keyed off
      *where* the exception came from rather than its type (a type-based
      allowlist is not safe for every shape a validation error can take,
      e.g. plain ``KeyError`` under ``validate=False``).
    - A bare ``ValueError``/``skbio.tree.TreeError``, for any other
      ``gpu_func`` (e.g. in tests exercising this dispatcher directly) that
      raises one of those two types itself.

    Either way, these are the caller's fault, not the backend's, so the
    original exception is re-raised immediately rather than being
    misclassified as "this backend's kernel can't run" -- which would both
    hide the real error behind a fallback and permanently (for the rest of
    the process) degrade every later call to the slow array-API path over
    one bad input.
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
