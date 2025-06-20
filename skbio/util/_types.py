# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Union

import pandas as pd
import numpy as np

from skbio.table import Table
from skbio.util import get_package

# Base types which are always available
DataTable = Union[pd.DataFrame, np.ndarray, Table]

# add other types depending on availability
pl = get_package("polars", raise_error=False)
if pl is not None:
    DataTable = Union[DataTable, pl.DataFrame]  # type: ignore[misc]

adt = get_package("anndata", raise_error=False)
if adt is not None:
    DataTable = Union[DataTable, adt.AnnData]  # type: ignore[misc]
