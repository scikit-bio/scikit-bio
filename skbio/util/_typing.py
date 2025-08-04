# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, Union, Protocol, runtime_checkable, Tuple, Any

import pandas as pd
import numpy as np
from numpy.typing import ArrayLike as NPArrayLike

from skbio.table import Table
from skbio.util import get_package


# ------------------------------------------------
# ArrayLike (see: skbio.util._array)
# ------------------------------------------------


@runtime_checkable
class StdArray(Protocol):
    r"""Any object compliant with the Python array API standard [1]_.

    Examples are numpy.ndarray, cupy.ndarray, torch.Tensor, jax.Array,
    dask.array.Array, and sparse.SparseArray.

    See Also
    --------
    ._array._get_array

    References
    ----------
    .. [1] https://data-apis.org/array-api/latest/

    """

    def __array_namespace__(self, api_version: Optional[str] = None): ...

    @property
    def shape(self) -> Tuple[int, ...]: ...

    @property
    def ndim(self) -> int: ...

    @property
    def dtype(self) -> Any: ...

    @property
    def device(self) -> Any: ...


ArrayLike = Union[NPArrayLike, StdArray]


# ------------------------------------------------
# TableLike (see: skbio.table._tabular)
# ------------------------------------------------

# Base types which are always available
TableLike = Union[pd.DataFrame, np.ndarray, Table]

# add other types depending on availability
pl = get_package("polars", raise_error=False)
if pl is not None:
    TableLike = Union[TableLike, pl.DataFrame]  # type: ignore[misc]

adt = get_package("anndata", raise_error=False)
if adt is not None:
    TableLike = Union[TableLike, adt.AnnData]  # type: ignore[misc]


# ------------------------------------------------
# SeedLike (see: skbio.util.get_rng)
# ------------------------------------------------

SeedLike = Union[int, np.random.Generator, np.random.RandomState]
