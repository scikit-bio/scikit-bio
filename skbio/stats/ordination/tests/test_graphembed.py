
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pytest
import numpy as np
from scipy.sparse import csr_matrix

from skbio.stats.ordination import graphembed


class TestGraphEmbed(unittest.TestCase):
    def test_graphembed_missing_dependency(self):
        try:
            import graphembed_rs  # noqa: F401
            has_ge = True
        except ImportError:
            has_ge = False

        adj = np.array([[0, 1], [1, 0]])

        if not has_ge:
            with self.assertRaises(ImportError):
                graphembed(adj, method="sketching")
        else:
            res = graphembed(adj, method="sketching", dimensions=2, nbiter=1)
            self.assertEqual(res.samples.shape, (2, 2))


@pytest.mark.skipif(True, reason="Run manually or dynamically check dependency")
def test_graphembed_sparse():
    import graphembed_rs  # noqa: F401
    
    adj = csr_matrix([[0, 1.5], [1.5, 0]])
    res = graphembed(adj, method="hope", dimensions=2)
    assert res.samples.shape == (2, 2)


if __name__ == "__main__":
    unittest.main()
