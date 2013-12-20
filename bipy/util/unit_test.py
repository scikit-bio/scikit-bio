#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from unittest import main, TestCase as PyTestCase
from math import log10
from numpy import asarray, isfinite, zeros, ravel
from numpy.testing import assert_almost_equal, assert_allclose
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
            # handle case where all elements were the same
            if not observed.shape:
                observed = observed.reshape((1,))
            if not expected.shape:
                expected = expected.reshape((1,))
            if observed[0] == expected[0]:
                return
        else:
            raise self.failureException(msg or 'p-value %s, t-test p %s' %
                                        (repr(pvalue), repr(p)))

    def assertSimilarFreqs(self, observed, expected, pvalue=0.01, msg=None):
        """Fail if observed p is lower than pvalue"""
        if self._suite_pvalue:
            pvalue = self._suite_pvalue

        obs_ravel = ravel(asarray(observed))
        exp_ravel = ravel(asarray(expected))

        m = zeros((2, len(obs_ravel)))
        m[0, :] = obs_ravel
        m[1, :] = exp_ravel

        G, p, _, _ = chi2_contingency(m, lambda_="log-likelihood")

        if p > pvalue:
            return
        else:
            raise self.failureException(msg or 'p-value %s, G-test p %s' %
                                        (repr(pvalue), repr(p)))

    def assertIsProb(self, observed, msg=None):
        """Fail is observed is not between 0.0 and 1.0"""
        try:
            if observed is None:
                raise ValueError
            if (asarray(observed) >= 0.0).all() and \
               (asarray(observed) <= 1.0).all():
                return
        except:
            pass
        raise self.failureException(msg or 'Observed %s has elements that are '
                                    'not probs' % (repr(observed)))

    def assertFloatEqual(self, observed, expected, eps=1e-6):
        """Tests whether two floating point numbers are approximately equal.

        If one of the arguments is zero, tests the absolute magnitude of the
        difference; otherwise, tests the relative magnitude.

        Use this method as a reasonable default.

        Conveniently wraps numpy.testing.assert_almost_equal
        """
        # don't change the eps interface that assertFloatEqual provides and
        # calculate the number of decimal digits that the values are compared
        # to
        assert_almost_equal(observed, expected, decimal=int(abs(log10(eps))))

    def assertFloatEqualRel(self, observed, expected, eps=1e-6):
        """Tests whether two floating point numbers/arrays are approx. equal.

        Checks whether the distance is within epsilon relative to the value
        of the sum of observed and expected. Use this method when you expect
        the difference to be small relative to the magnitudes of the observed
        and expected values.

        Conveniently wraps numpy.testing.assert_allclose with the rtol argument
        """
        assert_allclose(observed, expected, rtol=eps)

    def assertFloatEqualAbs(self, observed, expected, eps=1e-6):
        """
        Tests whether two floating point numbers are approximately equal.

        Checks whether the absolute value of (a - b) is within epsilon. Use
        this method when you expect that one of the values should be very
        small, and the other should be zero.

        Conveniently wraps numpy.testing.assert_allclose with the atol argument
        """
        assert_allclose(observed, expected, atol=eps)

    def assertEqualItems(self, observed, expected, msg=None):
        """Fail if the two items contain unequal elements"""
        obs_items = list(observed)
        exp_items = list(expected)
        if len(obs_items) != len(exp_items):
            raise self.failureException(
                msg or 'Observed and expected are different lengths: %s and %s'
                % (len(obs_items), len(exp_items)))

        obs_items.sort()
        exp_items.sort()
        for index, (obs, exp) in enumerate(zip(obs_items, exp_items)):
            if obs != exp:
                raise self.failureException(
                    msg or 'Observed %s and expected %s at sorted index %s' %
                    (obs, exp, index))

    def assertNotEqualItems(self, observed, expected, msg=None):
        """Fail if the two items contain only equal elements when sorted"""
        try:
            self.assertEqualItems(observed, expected, msg)
        except:
            pass
        else:
            raise self.failureException(
                msg or 'Observed %s has same items as %s' %
                (repr(observed), repr(expected)))
