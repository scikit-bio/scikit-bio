#!/usr/bin/env python
"""
Fit a Function to a Model (:mod:`skbio.math.fit_function`)
==========================================================

.. currentmodule:: skbio.maths.fit_function

This module provides a function to fit ``x`` and ``y`` samples to a model.

Functions
---------

.. autosummary::
   :toctree: generated/

   fit_function

"""
from __future__ import division

import numpy as np
from scipy.optimize import fmin


def fit_function(x_vals, y_vals, func, n_params, iterations=2):
    """Fit any function to any array of values of x and y.

    Parameters
    ----------

    x_vals : array_like
        Values for x to fit the function func.
    y_vals : array_like
        Values for y to fit the function func.
    func : callable ``f(x, a)``
        Objective function (model) to be fitted to the data. This function
        should return either an array for models that are not a constant,
        i.e. ``f(x)=exp(a[0]+x*a[1])``, or a single value for models that are a
        cosntant, i.e. ``f(x)=a[0]``
    n_params : int
        Number of parameters to fit in func
    iterations : int, optional
        Number of iterations to fit func

    Returns
    -------

    param_guess : array
        Values for each of the arguments to fit func to x_vals and y_vals

    Notes
    -----

    Fit a function to a given array of values x and y using simplex to minimize
    the error.

    Examples
    --------

    >>> from skbio.maths.fit_function import fit_function
    >>> from numpy import array
    >>> def f(x,a):
    ...    return a[0]
    >>> fit_function(array([1, 2, 4, 6]), array([0.5, 22, 5]), f, 1, 6)
    array([ 9.1666875])

   """

    # internal function to minimize the error
    def f2min(a):
        # sum square deviation
        return ((func(x_vals, a) - y_vals) ** 2).sum()

    param_guess = np.array(range(n_params))
    for i in range(iterations):
        xopt = fmin(f2min, param_guess, disp=0)
        param_guess = xopt

    return xopt
