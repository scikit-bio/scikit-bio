# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Array-like objects and namespaces."""

from typing import Union, Tuple, TYPE_CHECKING

import numpy as np
import array_api_compat as aac


if TYPE_CHECKING:  # pragma: no cover
    from types import ModuleType
    from ._typing import ArrayLike, StdArray


def _ingest_array(arr: "ArrayLike", /, *, to_numpy: bool = False) -> "StdArray":
    r"""Convert an array-like variable into an array object."""
    # Cast a non-array object into a NumPy array.
    if not aac.is_array_api_obj(arr):
        arr = np.asarray(arr)

    # Cast a non-NumPy array into a NumPy array.
    elif to_numpy and not aac.is_numpy_array(arr):
        # The Python array API standard requires any compliant array objects have a
        # `__dlpack__` protocol, which guarantees that they can be converted into a
        # NumPy array via `np.from_dlpack`.
        try:
            arr = np.from_dlpack(arr)

        # For array objects that doesn't have the `__dlpack__` protocol (e.g., Dask,
        # at the time of coding), `np.asarray` is called as a fallback, which should
        # work for objects that have an `__array__` protocol.
        except (AttributeError, TypeError):
            arr = np.asarray(arr)

    return arr


def ingest_array(
    *arrays: "ArrayLike", to_numpy: bool = False
) -> Tuple["ModuleType", Union["StdArray", "ArrayLike"]]:
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
        The array API compatible namespace corresponding ``arrays``.
    *arrays : object
        One or more array objects that are consistent with ``xp``.

    See Also
    --------
    _ingest_array
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

    ``xp`` is the namespace wrapper for different array libraries.

    If the input ``arr`` is already a compatible object, it will be returned as-is.
    Otherwise, it will be converted into a NumPy array.

    References
    ----------
    .. [1] https://data-apis.org/array-api/latest/

    .. [2] https://numpy.org/devdocs/glossary.html#term-array_like

    """
    arrays = tuple(_ingest_array(x, to_numpy=to_numpy) for x in arrays)
    return aac.array_namespace(*arrays), *arrays
