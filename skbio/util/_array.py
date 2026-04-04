# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Array-like objects and namespaces."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import array_api_compat as aac


if TYPE_CHECKING:  # pragma: no cover
    from types import ModuleType
    from ._typing import ArrayLike, StdArray


# -------------------------------------------------------------------------------------
# Important note: This script provides support for array-like objects. Two notions are
# currently used interchangeably in scikit-bio:
#
# 1. `np.typing.ArrayLike`: Used in most legacy code. It is usually converted into a
#    NumPy array via `np.asarray`.
#
# 2. `skbio.util._array.StdArray`: Used in some new code that focuses on GPU support.
#    It is essentially a union of NumPy array-like and array objects compliant with the
#    Python array API standard. It is parsed using the `_get_array` function below.
#
# The two types and the two functions can be used interchangeably, but with certain
# considerations, as detailed below.
# -------------------------------------------------------------------------------------


def _get_array(arr: ArrayLike, /, *, to_numpy: bool = False) -> StdArray:
    r"""Convert an array-like variable into an array object.

    Parameters
    ----------
    arr : array_like
        An array-like variable.
    to_numpy : bool, optional
        If True, make sure the returned array is a NumPy array.

    Returns
    -------
    arr : ndarray
        An array object.

    Notes
    -----
    When `to_numpy=True`, this function's behavior is very similar to `np.asarray`.
    Therefore, in most legacy code, a simple `np.asarray` call can replace the
    current function. However, there is a subtle difference that makes the current
    function more suitable for GPU-oriented code: `np.asarray` forces a device-to-
    host copy, whereas `np.from_dlpack` attempts to share the buffer without making
    a copy, thereby improving performance.

    """
    # Cast a non-array object into a NumPy array.
    if not aac.is_array_api_obj(arr):
        arr = np.asarray(arr)

    # Cast a non-NumPy array into a NumPy array.
    elif to_numpy and not aac.is_numpy_array(arr):
        # The Python array API standard requires any compliant array objects have a
        # `__dlpack__` protocol, which guarantees that they can be converted into a
        # NumPy array via `np.from_dlpack`.
        # See: https://numpy.org/devdocs//user/basics.interoperability.html
        try:
            arr = np.from_dlpack(arr)

        # For array objects that doesn't have the `__dlpack__` protocol (e.g., Dask,
        # at the time of coding), or array objects that are not on the same device
        # (e.g., cuda), `np.asarray` is called as a fallback, which should work for
        # objects that have an `__array__` protocol.
        except (AttributeError, TypeError, RuntimeError):
            arr = np.asarray(arr)

    return arr


# -------------------------------------------------------------------------------------
# Array API backend helpers
# -------------------------------------------------------------------------------------

_BACKEND_CHECKERS = {
    "numpy": aac.is_numpy_namespace,
    "cupy": aac.is_cupy_namespace,
    "torch": aac.is_torch_namespace,
    "jax": aac.is_jax_namespace,
    "dask": aac.is_dask_namespace,
}


def _get_namespace_from_args(args, kwargs):
    r"""Return the array namespace if any argument is an array object.

    Parameters
    ----------
    args : tuple
        Positional arguments.
    kwargs : dict
        Keyword arguments.

    Returns
    -------
    namespace or None
        A uniform array namespace for all detected array objects, or `None` if no array
        object is found.

    Raises
    ------
    TypeError
        If `args` or `kwargs` contain arrays from different array libraries.

    """
    arrays = [
        obj for obj in list(args) + list(kwargs.values()) if aac.is_array_api_obj(obj)
    ]
    if not arrays:
        return None
    return aac.array_namespace(*arrays)


def _check_array_api_backend(xp, backends, func_name):
    r"""Raise TypeError if xp is not among the declared supported backends.

    Parameters
    ----------
    xp : namespace
        Array namespace to check.
    backends : list of str
        Allowed backend names (e.g. `['numpy', 'torch']`).
    func_name : str
        Name of the calling function, used in the error message.

    Raises
    ------
    TypeError
        If `xp` is not a supported backend.

    """
    for name in backends:
        checker = _BACKEND_CHECKERS.get(name)
        if checker is not None and checker(xp):
            return
    supported = ", ".join(backends)
    raise TypeError(
        f"{func_name}() received an array from an unsupported backend. "
        f"Supported backends are: {supported}."
    )


# -------------------------------------------------------------------------------------
# Conversion and device helpers (used primarily in testing)
# -------------------------------------------------------------------------------------


def _to_numpy(arr):
    r"""Convert any array-API-compatible array to a NumPy array on CPU.

    Unlike `_get_array(arr, to_numpy=True)`, this handles GPU arrays
    explicitly without relying on `__dlpack__` or `__array__`, which can
    fail for device arrays (e.g., `np.asarray` raises `RuntimeError` on a
    PyTorch CUDA tensor).

    Parameters
    ----------
    arr : array
        Any array-API-compatible array object, or a NumPy array.

    Returns
    -------
    numpy.ndarray
        The data as a NumPy array on CPU.

    Notes
    -----
    This is intended for use in **test assertions**, not in library functions.
    Library functions should stay on the input device and use the `xp`
    namespace throughout.

    """
    if isinstance(arr, np.ndarray):
        return arr

    # PyTorch: must detach from computation graph and move to CPU first.
    # Without .detach().cpu(), np.asarray raises RuntimeError for CUDA tensors
    # and tensors that require grad.
    if hasattr(arr, "detach"):
        return arr.detach().cpu().numpy()

    # CuPy: .get() transfers GPU → CPU and returns a NumPy array.
    if hasattr(arr, "get"):
        return arr.get()

    # JAX, Dask, and others: np.asarray works (JAX arrays are on-host by
    # default; Dask triggers compute).
    return np.asarray(arr)


def _move_to_device(arr, xp, device):
    r"""Move an array to the specified device.

    Parameters
    ----------
    arr : array
        An array-API-compatible array.
    xp : module
        The array namespace (e.g. `numpy`, `torch`, `jax.numpy`).
    device : str or None
        Target device string (`'cpu'`, `'cuda'`, `'gpu'`, etc.).
        If `None` or `'cpu'`, the array is returned unchanged.

    Returns
    -------
    array
        The array on the requested device.

    Notes
    -----
    This is intended for use in **test setup**, to place test data on the
    correct device before calling the function under test.

    """
    if device is None or device == "cpu":
        return arr

    xp_name = xp.__name__.split(".")[0]

    if xp_name == "torch":
        return arr.to(device)

    if xp_name == "jax":
        import jax

        dev = jax.devices(device)[0]
        return jax.device_put(arr, dev)

    # CuPy arrays are always on GPU; numpy arrays are always on CPU.
    return arr


def _get_backend_name(xp):
    r"""Return a short backend name string for an array namespace.

    Parameters
    ----------
    xp : module
        An array namespace.

    Returns
    -------
    str
        One of `'numpy'`, `'torch'`, `'jax'`, `'cupy'`, `'dask'`,
        or the raw module name if unrecognized.

    """
    for name, checker in _BACKEND_CHECKERS.items():
        if checker(xp):
            return name
    return getattr(xp, "__name__", str(xp))


# -------------------------------------------------------------------------------------
# Main public helper
# -------------------------------------------------------------------------------------


def ingest_array(
    *arrays: ArrayLike, to_numpy: bool = False
) -> tuple[ModuleType, StdArray | ArrayLike]:
    r"""Convert array-like variables into array objects and their shared namespace.

    Parameters
    ----------
    *arrays : array_like
        One or more array-like variables.
    to_numpy : bool, optional
        If True, make sure the returned arrays are NumPy arrays.

    Returns
    -------
    xp : namespace
        The array API compatible namespace corresponding `arrays`.
    *arrays : object
        One or more array objects that are consistent with `xp`.

    See Also
    --------
    _get_array
    skbio.util._typing.ArrayLike
    array_api_compat.array_namespace
    numpy.typing.ArrayLike
    numpy.from_dlpack

    Notes
    -----
    A compatible array-like object must be either:

    1. An array object that is compatible with the Python array API standard [1]_.

        * Examples are numpy.ndarray, cupy.ndarray, torch.Tensor, jax.Array,
          dask.array.Array, and sparse.SparseArray.

    2. An object that can be casted into a NumPy array [2]_.

        * Examples are numpy.ndarray, Python list, tuple, array, and scalar.

    `xp` is the namespace wrapper for different array libraries.

    If the input `arr` is already a compatible object, it will be returned as-is.
    Otherwise, it will be converted into a NumPy array.

    References
    ----------
    .. [1] https://data-apis.org/array-api/latest/

    .. [2] https://numpy.org/devdocs/glossary.html#term-array_like

    """
    arrays = tuple(_get_array(x, to_numpy=to_numpy) for x in arrays)
    return aac.array_namespace(*arrays), *arrays
