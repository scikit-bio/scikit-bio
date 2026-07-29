# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Shared helpers for the single-source Numba CUDA/HIP GPU backends.

The permutation-test statistics (PERMANOVA, Mantel) each ship a fused kernel that
compiles through ``numba.cuda`` on NVIDIA and ``numba.hip`` on AMD. This module
holds the piece they share: choosing the Numba GPU module for a given
device-resident array (``numba.cuda`` for CuPy / CUDA-PyTorch, ``numba.hip`` for
ROCm-PyTorch), and a correctness-first fallback to the array-API path whenever the
fused kernel is unavailable or fails to build on the running stack.
"""

from warnings import warn

import array_api_compat as _aac

from skbio.util._array import _get_backend_name

# One thread block per permutation; _TPB threads walk the columns (the tile width
# used by scikit-bio-binaries).
_TPB = 128

# Backend names whose fused Numba GPU kernel failed to build or run in this
# process. Populated by ``_mark_gpu_unavailable`` so that later calls skip the
# kernel and take the array-API path instead of retrying a failing compilation.
_unavailable = set()


def _numba_gpu_module_for(arr):
    """Return the Numba GPU module for ``arr``'s device, or None.

    CuPy and CUDA-built PyTorch map to ``numba.cuda``; ROCm-built PyTorch maps to
    ``numba.hip``. Returns None for any other namespace (e.g. JAX, Dask), an
    unavailable backend, when Numba GPU support is not installed, or when this
    backend's fused kernel has already failed once in this process (see
    :func:`_mark_gpu_unavailable`); the caller then takes the array-API path. A
    CuPy array is routed to ``numba.cuda``, never ``numba.hip`` (the two target
    different GPU vendors).

    Parameters
    ----------
    arr : array
        A non-NumPy, array-API-compatible buffer (the DistanceMatrix data).

    Returns
    -------
    module or None
        The Numba GPU module usable on this array's device, else None.
    """
    name = _get_backend_name(_aac.array_namespace(arr))
    if name in _unavailable:
        return None
    if name == "cupy":
        want = "cuda"
    elif name == "torch":
        import torch

        # torch reports the same namespace for CUDA and ROCm builds; the build's
        # own hip version tag disambiguates them (ROCm -> numba.hip).
        want = "hip" if torch.version.hip is not None else "cuda"
    else:
        return None

    try:
        mod = getattr(__import__("numba", fromlist=[want]), want)
    except Exception:
        return None
    try:
        return mod if mod.is_available() else None
    except Exception:
        return None


def _mark_gpu_unavailable(arr):
    """Record that ``arr``'s backend cannot run the fused kernel this process.

    Called by the statistic dispatchers when the fused kernel raises (for example
    a numba-hip build that fails to compile on the running ROCm stack). Warns once
    per backend, then routes that backend to the array-API path from then on.
    """
    name = _get_backend_name(_aac.array_namespace(arr))
    if name not in _unavailable:
        _unavailable.add(name)
        warn(
            f"The Numba GPU kernel could not be used for the '{name}' backend on "
            "this system; using the array-API fallback instead.",
            UserWarning,
        )


def _get_kernel(cache, builder, gpu, backend_name):
    """Return the compiled kernel for a Numba GPU module, building it once.

    ``cache`` is a per-statistic dict (backend name -> compiled kernel); ``builder``
    compiles the kernel for a given Numba GPU module.
    """
    if backend_name not in cache:
        cache[backend_name] = builder(gpu)
    return cache[backend_name]
