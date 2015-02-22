from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np

from skbio.stats.composition import (closure, multiplicative_replacement,
                                     perturb, perturb_inv, power,
                                     clr, centralize)

class CompositionTests(TestCase):


    def setUp(self):
        self.data1 = np.vstack((
            np.array([2, 2, 6]),
            np.array([4, 4, 2])))
        self.data2 = np.array([2, 2, 6])
        
    def test_closure(self):        
        np.testing.assert_array_almost_equal(closure(self.data1),
                                          np.vstack((
                                              np.array([.2, .2, .6]),
                                              np.array([.4, .4, .2]))))
        np.testing.assert_array_almost_equal(closure(self.data2),
                                             np.array([.2, .2, .6]))

    def test_perturb(self):
        pmat = perturb(closure(self.data1), np.array([1, 1, 1]))
        np.testing.assert_array_almost_equal(pmat,
                                             np.vstack((
                                                 np.array([.2, .2, .6]),
                                                 np.array([.4, .4, .2]))))
        
        pmat = perturb(closure(self.data1), np.array([10, 10, 20]))
        np.testing.assert_array_almost_equal(pmat,
                                             np.vstack((
                                                 np.array([.125, .125, .75]),
                                                 np.array([1./3, 1./3, 1./3]))))

        pmat = perturb(closure(self.data1), np.array([10, 10, 20]))
        np.testing.assert_array_almost_equal(pmat,
                                             np.vstack((
                                                 np.array([.125, .125, .75]),
                                                 np.array([1./3, 1./3, 1./3]))))

        pmat = perturb(closure(self.data2), np.array([1, 2, 1]))
        np.testing.assert_array_almost_equal(pmat,np.array([1./6, 2./6, 3./6]))
                                                 

if __name__ == "__main__":
    unittest.main()
