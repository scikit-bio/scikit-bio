#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from bipy.util.unit_test import TestCase, main
from numpy import array

### ROB, GAVIN??, DANIEL, SANDRA??

class TestCaseTests(TestCase):
    """Tests for extension of the built-in unittest framework.

    For each test, includes an example of success and failure.
    """
    def test_assertSimilarMeans_one_obs_true(self):
        """assertSimilarMeans should raise on a single obs"""
        obs = [5]
        expected = [1,2,3,4,5,6,7,8,9,10,11]
        self.assertSimilarMeans(obs, expected)
        self.assertSimilarMeans(obs, expected, pvalue=0.25)
        self._set_suite_pvalue(0.10)
        self.assertSimilarMeans(obs, expected)

    def test_assertSimilarMeans_one_obs_fail(self):
        """assertSimilarMeans should raise on a single obs"""
        obs = [5]
        expected = [.001,.009,.00012, .000111, .002]

        with self.assertRaises(AssertionError):
            self.assertSimilarMeans(obs, expected)
        
        with self.assertRaises(AssertionError): 
            self.assertSimilarMeans(obs, expected, 0.1)
        
        self._set_suite_pvalue(0.001)
        with self.assertRaises(AssertionError): 
            self.assertSimilarMeans(obs, expected)

    def test_assertSimilarMeans_twosample_true(self):
        """assertSimilarMeans should pass when p > pvalue"""
        obs = [4,5,6]
        expected = [1,2,3,4,5,6,7,8,9]
        self.assertSimilarMeans(obs, expected)
        self.assertSimilarMeans(obs, expected, pvalue=0.25)
        self._set_suite_pvalue(0.10)
        self.assertSimilarMeans(obs, expected)

    def test_assertSimilarMeans_twosample_false(self):
        """assertSimilarMeans should raise when p < pvalue"""
        obs = [1,2,3]
        expected = [6,7,8,9,10,11,12,13,14]
        with self.assertRaises(AssertionError): 
            self.assertSimilarMeans(obs, expected)
        
        with self.assertRaises(AssertionError):
            self.assertSimilarMeans(obs, expected, 0.1)
        
        self._set_suite_pvalue(0.001)
        with self.assertRaises(AssertionError):
            self.assertSimilarMeans(obs, expected)

    def test_assertSimilarFreqs_true(self):
        """assertSimilarFreqs should pass when p > pvalue"""
        observed = [2,2,3,2,1,2,2,2,2]
        expected = [2,2,2,2,2,2,2,2,2]
        self.assertSimilarFreqs(observed, expected)
        self.assertSimilarFreqs(observed, expected, pvalue=0.25)
        self._set_suite_pvalue(0.10)
        self.assertSimilarFreqs(observed, expected)

    def test_assertSimilarFreqs_false(self):
        """assertSimilarFreqs should raise when p < pvalue"""
        observed = [10,15,20,10,12,12,13]
        expected = [100,50,10,20,700,2,100]
        with self.assertRaises(AssertionError):
            self.assertSimilarFreqs(observed, expected)
        
        with self.assertRaises(AssertionError): 
            self.assertSimilarFreqs(observed, expected, 0.2)
        
        self._set_suite_pvalue(0.001)
        with self.assertRaises(AssertionError):
            self.assertSimilarFreqs(observed, expected)

    def test_assertSimilarFreqs_numpy_array_true(self):
        """assertSimilarFreqs should pass when p > pvalue"""
        observed = array([2,2,3,2,1,2,2,2,2])
        expected = array([2,2,2,2,2,2,2,2,2])
        self.assertSimilarFreqs(observed, expected)
        self.assertSimilarFreqs(observed, expected, pvalue=0.25)
        self._set_suite_pvalue(0.10)
        self.assertSimilarFreqs(observed, expected)

    def test_assertSimilarFreqs_numpy_array_false(self):
        """assertSimilarFreqs should raise when p < pvalue"""
        observed = array([10,15,20,10,12,12,13])
        expected = array([100,50,10,20,700,2,100])
        with self.assertRaises(AssertionError):
            self.assertSimilarFreqs(observed, expected)
        
        with self.assertRaises(AssertionError): 
            self.assertSimilarFreqs(observed, expected, 0.2)
        
        self._set_suite_pvalue(0.001)
        with self.assertRaises(AssertionError):
            self.assertSimilarFreqs(observed, expected)

    def test_set_suite_pvalue(self):
        """Should set the suite pvalue"""
        # force stats to fail
        self._set_suite_pvalue(0.99)
        obs = [2,5,6]
        exp = [1,2,3,4,5,6,7,8,9]
        with self.assertRaises(AssertionError):
            self.assertSimilarMeans(obs, exp)

        # force stats to pass
        self._set_suite_pvalue(0.01)
        self.assertSimilarMeans(obs, exp)

if __name__ == '__main__':
    main()
