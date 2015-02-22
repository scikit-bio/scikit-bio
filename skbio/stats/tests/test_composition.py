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

from skbio.stats.composition import (closure, zero_imputation,
                                     perturb, perturb_inv, power,
                                     clr, centralize)

class CompositionTests(TestCase):


    def setUp(self):
        self.data1 = np.vstack((
            np.array([2, 2, 6]),
            np.array([4, 4, 2])))
        self.data2 = np.array([2, 2, 6])
        D = 5
        self.data3 = np.vstack((
            np.array(range(1,4) + [0]*1 + [5]),
            np.array(range(1,2) + [0]*2 + range(4,6)),
            np.array(range(1,D+1))))
        
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
        
    def test_power(self):
        pmat = power(closure(self.data1), 2)
        np.testing.assert_array_almost_equal(pmat,
                                             np.vstack((
                                                 np.array([.04, .04, .36])/.44,
                                                 np.array([.16, .16, .04])/.36)))


    def test_perturb_inv(self):
        pmat = perturb_inv(closure(self.data1), np.array([.1, .1, .1]))
        imat = perturb(closure(self.data1), np.array([10, 10, 10]))        
        np.testing.assert_array_almost_equal(pmat, imat)


        pmat = perturb_inv(closure(self.data1), np.array([1, 1, 1]))
        np.testing.assert_array_almost_equal(pmat,
                                             np.vstack((
                                                 np.array([.2, .2, .6]),
                                                 np.array([.4, .4, .2]))))
    
    def test_multiplicative_replacement(self):
        amat = zero_imputation(self.data3)
        np.testing.assert_array_almost_equal(amat,
                np.array([[ 0.09056604,  0.18113208,  0.27169811,  0.00377358,  0.45283019],
                          [ 0.09913793,  0.00431034,  0.00431034,  0.39655172,  0.49568966],
                          [ 0.06666667,  0.13333333,  0.2       ,  0.26666667,  0.33333333]]))


    def test_clr(self):
        cmat = clr(closure(self.data1))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])

        np.testing.assert_array_almost_equal(cmat,
                                             np.vstack((
                                                 np.log(A / np.exp(np.log(A).mean())) ,
                                                 np.log(B / np.exp(np.log(B).mean())) )))

    
if __name__ == "__main__":
    unittest.main()
