#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from unittest import main, TestCase as PyTestCase
from numpy import asarray, isfinite, zeros, ravel
from scipy.stats.stats import ttest_ind, ttest_1samp
from scipy.stats.contingency import chi2_contingency 

class TestCase(PyTestCase):
    _suite_pvalue = None
    def _set_suite_pvalue(self, pvalue):
        """Sets the test suite pvalue to be used in similarity tests

        This value is by default None. The pvalue used in this case is
        specified in the test module itself. The purpose of this method is to
        set the pvalue to be used when running a massive test suite
        """
        self._suite_pvalue = pvalue

    def assertSimilarMeans(self, observed, expected, pvalue=0.01, msg=None):
        """Fail if observed p is lower than pvalue"""
        if self._suite_pvalue:
            pvalue = self._suite_pvalue

        observed, expected = asarray(observed), asarray(expected)
        
        if observed.size == 1:
            t, p = ttest_1samp(expected, observed)
        elif expected.size == 1:
            t, p = ttest_1samp(observed, expected)
        elif observed.size == 1 and expected.size == 1:
            t, p = ttest_1samp(expected, observed)
        else:
            t, p = ttest_ind(observed, expected)
        
        if p > pvalue:
            return
        elif p is None or not isfinite(p): 
            #handle case where all elements were the same
            if not observed.shape:
                observed = observed.reshape((1,))
            if not expected.shape:
                expected = expected.reshape((1,))
            if observed[0] == expected[0]:
                return
        else:
            raise self.failureException(msg or 'p-value %s, t-test p %s' % \
                                        (`pvalue`, `p`))

    def assertSimilarFreqs(self, observed, expected, pvalue=0.01, msg=None):
        """Fail if observed p is lower than pvalue"""
        if self._suite_pvalue:
            pvalue = self._suite_pvalue

        obs_ravel = ravel(asarray(observed))
        exp_ravel = ravel(asarray(expected))

        m = zeros((2,len(obs_ravel)))
        m[0,:] = obs_ravel
        m[1,:] = exp_ravel

        G, p, _, _= chi2_contingency(m, lambda_="log-likelihood")

        if p > pvalue:
            return
        else:
            raise self.failureException(msg or 'p-value %s, G-test p %s' % \
                                        (`pvalue`, `p`))

