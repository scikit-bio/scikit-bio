#!/usr/bin/env python
"""Unit tests for fit function.
"""

from unittest import TestCase, main

import numpy as np
from numpy.random import rand
from numpy.testing import assert_allclose
from skbio.maths.fit_function import fit_function


class FitFunctionTests(TestCase):
    """Tests of top-level fit functions."""

    def test_constant(self):
        """test constant approximation"""
        # defining our fitting function
        def f(x,a):
            return a[0]

        exp_params = [2]
        x = np.arange(-1,1,.01)
        y = f(x, exp_params)
        y_noise = y + rand(len(x))

        params = fit_function(x, y_noise, f, 1, 5)
        assert_allclose(params, exp_params, 0.5)

    def test_linear(self):
        """test linear approximation"""
        # defining our fitting function
        def f(x,a):
            return (a[0]+x*a[1])

        exp_params = [2, 10]
        x = np.arange(-1,1,.01)
        y = f(x, exp_params)
        y_noise = y + rand(len(y))

        params = fit_function(x, y_noise, f, 2, 5)
        assert_allclose(params, exp_params, 0.5)

    def test_exponential(self):
        """test exponential approximation"""
        # defining our fitting function
        def f(x,a):
            return np.exp(a[0]+x*a[1])

        exp_params = [2, 10]
        x = np.arange(-1,1,.01)
        y = f(x, exp_params)
        y_noise = y + rand(len(y))
        
        params = fit_function(x, y_noise, f, 2, 5)
        assert_allclose(params, exp_params, 0.5)


if __name__ == '__main__':
    main()
