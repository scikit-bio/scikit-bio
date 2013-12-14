#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""Translations of functions from Release 2.3 of the Cephes Math Library, 
which is (c) Stephen L. Moshier 1984, 1995.
"""
from __future__ import division

from bipy.maths.stats.special import (fix_rounding_error, expm1, log1p, betai,
                                      igamc, erf, erfc, GB, SQRTH, LP, LQ, EQ)

def chi_high(x, df):
    """Returns right-hand tail of chi-square distribution (x to infinity).
    
    df, the degrees of freedom, ranges from 1 to infinity (assume integers).
    Typically, df is (r-1)*(c-1) for a r by c table.
    
    Result ranges from 0 to 1.
    
    See Cephes docs for details.
    """
    x = fix_rounding_error(x)
    
    if x < 0:
        raise ValueError, "chi_high: x must be >= 0 (got %s)." % x
    if df < 1:
        raise ValueError, "chi_high: df must be >= 1 (got %s)." % df
    return igamc(df/2, x/2)

def z_high(x):
    """Returns right-hand tail of z distribution (0 to x). 
    
    x ranges from -infinity to +infinity; result ranges from 0 to 1
    
    See Cephes docs for details."""
    y = x * SQRTH
    z = abs(y)
    if z < SQRTH:
        return 0.5 - 0.5 * erf(y)
    else:
        if x < 0:
            return 1 - 0.5 * erfc(z)
        else:
            return 0.5 * erfc(z)

def zprob(x):
    """Returns both tails of z distribution (-inf to -x, inf to x)."""
    return 2 * z_high(abs(x))

def f_high(df1, df2, x):
    """Returns right-hand tail of f distribution (x to infinity).
    
    Result ranges from 0 to 1.
    
    See Cephes docs for details.
    """
    return fdtrc(df1, df2, x)

def fdtrc(a, b, x):
    """Returns right tail of F distribution, x to infinity.

    See Cephes docs for details.
    """
    if min(a, b) < 1:
        raise ValueError, "F a and b (degrees of freedom) must both be >= 1."
    if x < 0:
        raise ValueError, "F distribution value of f must be >= 0."
    w = float(b) / (b + a * x)
    return betai(0.5 * b, 0.5 * a, w)

def binomial_high(successes, trials, prob):
    """Returns right-hand binomial tail (X > successes) given prob(success)."""
    if -1 <= successes < 0:
        return 1
    return bdtrc(successes, trials, prob)

def bdtrc(k, n, p):
    """Complement of binomial distribution, k+1 through n.

    Uses formula bdtrc(k, n, p) = betai(k+1, n-k, p)

    See Cephes docs for details.
    """
    p = fix_rounding_error(p)
    if (p < 0) or (p > 1):
        raise ValueError, "Binomial p must be between 0 and 1."
    if (k < 0) or (n < k):
        raise ValueError, "Binomial k must be between 0 and n."
    if k == n:
        return 0
    dn = n - k
    if k == 0:
        if p < .01:
            dk = -expm1(dn * log1p(-p))
        else:
            dk = 1 - pow(1.0-p, dn)
    else:
        dk = k + 1
        dk = betai(dk, dn, p)
    return dk

def pseries(a, b, x):
    """Power series for incomplete beta integral.

    Use when b * x is small and x not too close to 1.

    See Cephes docs for details.
    """
    ai = 1 / a
    u = (1-b) * x
    v = u / (a + 1)
    t1 = v
    t = u
    n = 2
    s = 0
    z = MACHEP * ai
    while abs(v) > z:
        u = (n - b) * x / n
        t *= u
        v = t / (a + n)
        s += v
        n += 1
    s += t1
    s += ai

    u = a * log(x)
    if ((a + b) < MAXGAM) and (abs(u) < MAXLOG):
        t = Gamma(a+b)/(Gamma(a)*Gamma(b))
        s = s * t * pow(x, a)
    else:
        t = lgam(a+b) - lgam(a) - lgam(b) + u + log(s)
        if t < MINLOG:
            s = 0
        else:
            s = exp(t)
    return(s)
