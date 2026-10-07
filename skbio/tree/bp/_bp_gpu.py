# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""GPU kernels of :class:`skbio.tree.BPTree` (Numba CUDA / HIP).

The kernels are written once and compiled for the Numba GPU module of the
tree's device (``numba_cuda_mlir.cuda`` or ``numba.cuda`` on NVIDIA,
``numba.hip`` on AMD; see :func:`._gpu._numba_gpu_module_for`), on first
use, from the device-function build of the navigation primitives
(:func:`._bp_numba.gpu_primitives`). They are reached through ``BPTree`` with
``engine="numba"`` on a tree whose ``data`` lives on a CUDA or ROCm device;
every other call runs on the CPU engines.

Nothing here imports a GPU module: the caller passes it in, so importing
scikit-bio never initializes a GPU.

Each function takes the module ``gpu``, the backend name (the kernel cache
key) and its arrays either on the device (anything exposing
``__cuda_array_interface__``, read in place) or on the host (uploaded). A
result is written into ``out`` when it is given, a device array of the
caller's backend, and is otherwise returned as a new NumPy array.
"""

import warnings

import numpy as np

from . import _bp_numba
from ._gpu import _TPB, _get_kernel, _sync_stream

# Compiled kernels, by backend name (see ``_get_kernel``).
_kernels = {}


def _build_kernels(gpu):
    """Compile the kernels for a Numba GPU module.

    ``gpu.jit`` compiles each kernel at its first launch, for the argument
    types of that launch.
    """
    G = _bp_numba.gpu_primitives(gpu)
    # device functions are bound to plain names: a kernel calls a dispatcher,
    # not an attribute of a tuple of them
    close, parent, lca, level_ancestor = G.close, G.parent, G.lca, G.level_ancestor

    # -- batch navigation: one thread per query ------------------------------

    @gpu.jit
    def close_batch(T, idx, out):  # pragma: no cover (runs on the device)
        t = gpu.grid(1)
        if t < out.shape[0]:
            out[t] = close(T, idx[t])

    @gpu.jit
    def parent_batch(T, idx, out):  # pragma: no cover
        t = gpu.grid(1)
        if t < out.shape[0]:
            out[t] = parent(T, idx[t])

    @gpu.jit
    def lca_batch(T, i, j, out):  # pragma: no cover
        t = gpu.grid(1)
        if t < out.shape[0]:
            out[t] = lca(T, i[t], j[t])

    @gpu.jit
    def level_ancestor_batch(T, idx, d, out):  # pragma: no cover
        t = gpu.grid(1)
        if t < out.shape[0]:
            out[t] = level_ancestor(T, idx[t], d[t])

    # -- cophenetic distances ----------------------------------------------------
    #
    # Every element of the square matrix is written once, in rows, as
    # (d[a] - d[v]) + (d[c] - d[v]) for tips a, c with lowest common ancestor
    # v: the CPU kernels' formula, whose two terms only swap between (a, c) and
    # (c, a). Floating-point addition is commutative, so the matrix is
    # bit-identical to the CPU engines' and symmetric, and with no multiply
    # there is no fused multiply-add to change a rounding.

    @gpu.jit
    def tip_distance_rows(
        tips, slot, parent, begin, end, dist, out
    ):  # pragma: no cover
        # one block per row: its threads walk the row's ancestors in step (the
        # same reads for all of them) and share the columns of each ancestor.
        # The tips of ancestor v are the sorted tips [begin[v], end[v]); those
        # outside the subtree already covered have v as their lca.
        x = gpu.blockIdx.x
        t = gpu.threadIdx.x
        step = gpu.blockDim.x
        n = tips.shape[0]
        r = slot[x]
        a = tips[x]
        da = dist[a]
        if t == 0:
            out[r, r] = 0.0
        lo = x
        hi = x + 1
        v = parent[a]
        while (lo > 0 or hi < n) and v >= 0:
            dv = dist[v]
            v_lo = begin[v]
            v_hi = end[v]
            for col in range(v_lo + t, lo, step):
                out[r, slot[col]] = (da - dv) + (dist[tips[col]] - dv)
            for col in range(hi + t, v_hi, step):
                out[r, slot[col]] = (da - dv) + (dist[tips[col]] - dv)
            lo = v_lo
            hi = v_hi
            v = parent[v]

    return {
        "close_batch": close_batch,
        "parent_batch": parent_batch,
        "lca_batch": lca_batch,
        "level_ancestor_batch": level_ancestor_batch,
        "tip_distance_rows": tip_distance_rows,
    }


def _launch(kernel, grid, block, *args):
    """Launch a kernel and wait for it.

    Numba warns about the low occupancy of a small grid on every launch, but
    a small batch or tree is a legitimate input here, not a misconfiguration.
    """
    with warnings.catch_warnings():
        # matched by message: each Numba GPU extension raises its own class of it
        warnings.filterwarnings("ignore", message="Grid size .* under-utilization")
        kernel[grid, block](*args)


def _on_device(gpu, arr):
    """A Numba device array of ``arr``: in place if on the device, else a copy.

    An array of a GPU backend is read in place once the work queued for it on
    that backend's stream is done (see :func:`._gpu._sync_stream`).
    """
    if hasattr(arr, "__cuda_array_interface__"):
        _sync_stream(arr)
        return gpu.as_cuda_array(arr)
    return gpu.to_device(np.ascontiguousarray(arr))


def tree_arrays(gpu, tree):
    """The :class:`._bp_numba.BPArrays` of a tree on the device.

    Uploaded from the tree's host index, which the device index would equal,
    and cached on the tree per GPU module.
    """
    cache = tree._gpu_arrays
    name = gpu.__name__
    if name not in cache:
        cache[name] = _bp_numba.bp_arrays(tree, asarray=gpu.to_device)
    return cache[name]


def run_batch(gpu, backend, T, name, args, out=None):
    """Run a batch navigation kernel, one thread per query.

    Parameters
    ----------
    gpu : module
        The Numba GPU module.
    backend : str
        The array backend name, keying the kernel cache.
    T : BPArrays
        The tree's arrays on the device (:func:`tree_arrays`).
    name : str
        The kernel: ``"close_batch"``, ``"parent_batch"``, ``"lca_batch"`` or
        ``"level_ancestor_batch"``.
    args : sequence of array of int64
        The kernel's flat query arrays, all of one length.
    out : array of int64, optional
        A device array of that length to write the result into.

    Returns
    -------
    array
        ``out`` if given, else a NumPy array of intp.
    """
    kernel = _get_kernel(_kernels, _build_kernels, gpu, backend)[name]
    n = args[0].shape[0]
    d_out = gpu.device_array(n, dtype=np.intp) if out is None else _on_device(gpu, out)
    if n:
        d_args = [_on_device(gpu, a) for a in args]
        _launch(kernel, -(-n // _TPB), _TPB, T, *d_args, d_out)
        gpu.synchronize()
    return d_out.copy_to_host() if out is None else out


def tip_distances(gpu, backend, tips, slot, parent, begin, end, dist, out=None):
    """Pairwise path distances between tips, as a square matrix.

    The GPU counterpart of ``_bp_cy.tip_distances``, with ``begin`` (the
    number of the sorted tips before each node) alongside ``end``.

    Returns
    -------
    array
        ``out`` if given, else a NumPy array of float64.
    """
    kernel = _get_kernel(_kernels, _build_kernels, gpu, backend)["tip_distance_rows"]
    n = tips.shape[0]
    d_out = (
        gpu.device_array((n, n), dtype=np.float64)
        if out is None
        else _on_device(gpu, out)
    )
    if n:
        args = [_on_device(gpu, a) for a in (tips, slot, parent, begin, end, dist)]
        # one block per row
        _launch(kernel, n, _TPB, *args, d_out)
        gpu.synchronize()
    return d_out.copy_to_host() if out is None else out
