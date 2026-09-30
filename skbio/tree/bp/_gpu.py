# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Shared helpers for the single-source Numba CUDA/HIP kernels of ``BPTree``.

The batch operations and ``cophenet`` of :class:`skbio.tree.BPTree` each have a
kernel (:mod:`._bp_gpu`) that compiles through ``numba_cuda_mlir.cuda`` or
``numba.cuda`` on NVIDIA and ``numba.hip`` on AMD. This module holds the piece
they share: choosing the Numba GPU module for a tree's device-resident data,
and a correctness-first fallback to the CPU engines whenever a kernel is
unavailable or fails to build on the running stack.

It follows :mod:`skbio.stats.distance._gpu`, the helpers of the permutation
tests, but keeps its own record of failed backends, so that a BPTree kernel
failing on a backend does not disable the PERMANOVA and Mantel kernels there,
or the reverse.
"""

from warnings import warn

import array_api_compat as _aac

from skbio.util._array import _get_backend_name

# Threads per block of the kernels.
_TPB = 128

# Backend names whose BPTree GPU kernels failed to build or run in this
# process. Populated by ``_mark_gpu_unavailable`` so that later calls skip the
# kernels and take the CPU engines instead of retrying a failing compilation.
_unavailable = set()


def _numba_gpu_module_for(arr):
    """Return the Numba GPU module for ``arr``'s device, or None.

    CUDA-built CuPy and PyTorch map to ``numba_cuda_mlir.cuda`` if it is installed,
    and to ``numba.cuda`` otherwise; ROCm-built CuPy and PyTorch map to
    ``numba.hip`` (both report the same array-API namespace, so the build's own
    flag disambiguates them). Returns None for any other namespace (e.g. JAX,
    Dask), an unavailable backend, when Numba GPU support is not installed, or when
    this backend's kernels have already failed once in this process (see
    :func:`_mark_gpu_unavailable`); the caller then takes the CPU engines.

    Parameters
    ----------
    arr : array
        A non-NumPy, array-API-compatible buffer (the tree's data).

    Returns
    -------
    module or None
        The Numba GPU module usable on this array's device, else None.
    """
    name = _get_backend_name(_aac.array_namespace(arr))
    if name in _unavailable:
        return None
    if name == "cupy":
        import cupy

        # cuPy reports the same namespace for CUDA and ROCm builds; the build's
        # own is_hip flag disambiguates them (ROCm -> numba.hip).
        want = "hip" if getattr(cupy.cuda.runtime, "is_hip", False) else "cuda"
    elif name == "torch":
        import torch

        # torch reports the same namespace for CUDA and ROCm builds; the build's
        # own hip version tag disambiguates them (ROCm -> numba.hip).
        want = "hip" if torch.version.hip is not None else "cuda"
    else:
        return None

    if want == "cuda":
        # numba-cuda is in maintenance mode and does not support NumPy 2.5, so its
        # successor, numba-cuda-mlir, is preferred when it is installed.
        mod = _available_module("numba_cuda_mlir", "cuda")
        if mod is not None:
            return mod
    return _available_module("numba", want)


def _available_module(package, name):
    """Return ``package.name`` if it imports and reports a usable GPU, else None."""
    try:
        mod = getattr(__import__(package, fromlist=[name]), name)
    except Exception:
        return None
    try:
        return mod if mod.is_available() else None
    except Exception:
        return None


def _mark_gpu_unavailable(arr):
    """Record that ``arr``'s backend cannot run the BPTree kernels this process.

    Called by ``BPTree`` when a kernel raises (for example a Numba GPU module
    that fails to compile on the running stack). Warns once per backend, then
    routes that backend to the CPU engines from then on. The warning begins as
    the permutation tests' does, so a filter on it covers both.
    """
    name = _get_backend_name(_aac.array_namespace(arr))
    if name not in _unavailable:
        _unavailable.add(name)
        warn(
            f"The Numba GPU kernel could not be used for the '{name}' backend on "
            "this system; using the CPU engines instead.",
            UserWarning,
        )


def _get_kernel(cache, builder, gpu, backend_name):
    """Return the compiled kernels for a Numba GPU module, building them once.

    ``cache`` is a dict (backend name -> compiled kernels); ``builder`` compiles
    the kernels for a given Numba GPU module.
    """
    if backend_name not in cache:
        cache[backend_name] = builder(gpu)
    return cache[backend_name]
