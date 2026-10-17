# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Shared QIIME 1.9.1 tiny-test fixture loading for unifrac tests."""

from unittest.mock import patch

import numpy as np
import pandas as pd

from skbio import TreeNode, DistanceMatrix
from skbio.util._testing import get_data_path


def patch_gpu_backend(backend):
    """Patch `_detect_gpu_backend` to report ``backend``.

    ``None`` forces the no-usable-backend branch of the `*_gpu_or_xp`
    dispatch wrappers (i.e. the array-API fallback), and 'cuda'/'hip' force
    the GPU branch, deterministically and independently of what hardware the
    machine running the tests actually has.
    """
    return patch(
        'skbio.diversity.beta._unifrac_gpu._detect_gpu_backend', return_value=backend
    )


class QiimeTinyTestMixin:
    """Mixin providing access to the qiime-191-tt fixture directory."""

    def _load_qiime_191_tt(self):
        base = get_data_path('qiime-191-tt', subfolder='data')
        df = pd.read_csv(f'{base}/otu-table.tsv', sep='\t', skiprows=1, index_col=0)
        taxa = df.index.astype(str).tolist()
        sample_ids = df.columns.tolist()
        table = df.T.values
        tree = TreeNode.read(f'{base}/tree.nwk')
        return table, taxa, tree, sample_ids

    def _load_dm_fixture(self, filename):
        base = get_data_path('qiime-191-tt', subfolder='data')
        return DistanceMatrix.read(f'{base}/{filename}')

    def _assert_matches_dm_fixture(self, condensed, sample_ids, filename, atol):
        """Assert a condensed vector matches a qiime-191-tt fixture matrix.

        Both matrices are filtered to ``sample_ids`` so the comparison is
        order-independent (the fixture files are not necessarily in the OTU
        table's sample order).
        """
        obs = DistanceMatrix(condensed, sample_ids)
        expected = self._load_dm_fixture(filename)
        np.testing.assert_allclose(
            obs.filter(sample_ids).data, expected.filter(sample_ids).data,
            rtol=0, atol=atol)
