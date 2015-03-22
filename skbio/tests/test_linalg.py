# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main
import numpy.testing as npt

import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg

from skbio.linalg import ssvd

class SSVDTests(TestCase):

    def setUp(self):
        np.random.seed(0)
        self.mat1 = np.random.randn(9, 6)
        self.mat2 = np.dot(self.mat1.T, self.mat1)

    def test_ssvd(self):
        np.random.seed(0)
        test_s, test_U, test_V = ssvd(mat1, k=3)
        actual_U, actual_s, actual_V = scipy.sparse.linalg.svds(mat1, k=3)
        npt.assert_allclose(-1*test_U, actual_U, rtol=4)
        npt.assert_allclose(test_s, actual_s, rtol=4)
        npt.assert_allclose(-1*test_V, actual_V, rtol=4)

    def test_ssvd_eig(self):
        np.random.seed(0)
        test_s, test_U = ssvd(mat2, k=5, compute_v=False)
        actual_s, actual_U = scipy.sparse.linalg.eigsh(mat2, k=5)
        npt.assert_allclose(-1*test_U, actual_U, rtol=4)
        npt.assert_allclose(test_s, actual_s, rtol=4)


if __name__=="__main__":
    main()
