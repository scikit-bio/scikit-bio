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
from skbio.util.config import _get_package

# Base types which are always available
DataTable = Union[pd.DataFrame, np.ndarray, Table]

# add other types depending on availability
pl, has_polars = _get_package("polars", raise_error=False, no_bool=False)
if has_polars:
    DataTable = Union[DataTable, pl.DataFrame]

adt, has_anndata = _get_package("anndata", raise_error=False, no_bool=False)
if has_anndata:
    DataTable = Union[DataTable, adt.AnnData]
