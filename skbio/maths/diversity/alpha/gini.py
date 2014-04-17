#!/usr/bin/env python
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from .base import _validate


def gini_index(data, method='rectangles'):
    """Calculates the gini index of data.

    Notes:
     formula is G = A/(A+B) where A is the area between y=x and the Lorenz curve
     and B is the area under the Lorenz curve. Simplifies to 1-2B since A+B=.5
     Formula available on wikipedia.
    Inputs:
     data - list or 1d arr, counts/abundances/proportions etc. All entries must
     be non-negative.
     method - str, either 'rectangles' or 'trapezoids'. see
     lorenz_curve_integrator for details.

    """
    # Suppress cast to int because this method supports ints and floats.
    data = _validate(data, suppress_cast=True)
    lorenz_points = _lorenz_curve(data)
    B = _lorenz_curve_integrator(lorenz_points, method)
    return 1 - 2 * B


def _lorenz_curve(data):
    """Calculates the Lorenz curve for input data.

    Notes:
     Formula available on wikipedia.
    Inputs:
     data - list or 1d arr, counts/abundances/proportions etc. All entries must
     be non-negative.

    """
    if any(np.array(data) < 0):
        raise ValueError("Lorenz curves aren't meaningful for non-positive "
                         "data.")

    # dont wan't to change input, copy and sort
    sdata = np.array(sorted((data[:])))
    n = float(len(sdata))
    Sn = sdata.sum()
    # ind+1 because must sum first point, eg. x[:0] = []
    lorenz_points = [((ind + 1) / n, sdata[:ind + 1].sum() / Sn)
                     for ind in range(int(n))]
    return lorenz_points


def _lorenz_curve_integrator(lc_pts, method):
    """Calculates the area under a lorenz curve.

    Notes:
     Could be utilized for integrating other simple, non-pathological
     'functions' where width of the trapezoids is constant.
     Two methods are available.
     1. Trapezoids, connecting the lc_pts by linear segments between them.
        Basically assumes that given sampling is accurate and that more features
        of given data would fall on linear gradients between the values of this
        data. formula is: dx[(h_0+h_n)/2 + sum(i=1,i=n-1,h_i)]
     2. Rectangles, connecting lc_pts by lines parallel to x axis. This is the
        correct method in my opinion though trapezoids might be desirable in
        some circumstances. forumla is : dx(sum(i=1,i=n,h_i))
    Inputs:
     lc_pts - list of tuples, output of lorenz_curve.
     method - str, either 'rectangles' or 'trapezoids'

    """
    if method is 'trapezoids':
        dx = 1. / len(lc_pts)  # each point differs by 1/n
        h_0 = 0.0  # 0 percent of the population has zero percent of the goods
        h_n = lc_pts[-1][1]
        sum_hs = sum([pt[1] for pt in lc_pts[:-1]])  # the 0th entry is at x=
        # 1/n
        return dx * ((h_0 + h_n) / 2. + sum_hs)
    elif method is 'rectangles':
        dx = 1. / len(lc_pts)  # each point differs by 1/n
        return dx * sum([pt[1] for pt in lc_pts])
    else:
        raise ValueError("Method '%s' not implemented. Available methods: "
                         "'rectangles', 'trapezoids'." % method)
