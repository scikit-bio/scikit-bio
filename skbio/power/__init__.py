r"""Empirical Power Estimation (:mod:`skbio.power`)
=====================================================

.. currentmodule:: skbio.power

The purpose of this module is to provide empirical, post-hoc power estimation
of normally and non-normally distributed data. It also provides support to
subsample data to facilitate this analysis.

The underlying principle is based on subsampling and Monte Carlo simulation.
Assume that there is some set of populations, :math:`K_{1}, K_{2}, ... K_{n}`
which have some property, :math:`\mu` such that :math:`\mu_{1} \neq \mu_{2}
\neq ... \neq \mu_{n}`. For each of the populations, a sample, :math:`S` can be
drawn, with a parameter, :math:`x` where :math:`x \approx \mu` and for the
samples, we can use a test, :math:`f`, to show that :math:`x_{1} \neq x_{2}
\neq ... \neq x_{n}`.

Since we know that :math:`\mu_{1} \neq \mu_{2} \neq ... \neq \mu_{n}`,
we know we should reject the null hypothesis. If we fail to reject the null
hypothesis, we have committed a Type II error and our result is a false
negative. We can estimate the frequency of Type II errors at various sampling
depths by repeatedly subsampling the populations and observing how often we
see a false negative. If we repeat this several times for each subsampling
depth, and vary the depths we use, we can start to approximate a relationship
between the number of samples we use and the rate of false negatives, also
called the statistical power of the test.

To generate complete power curves from data which appears underpowered, the
`statsmodels.stats.power` package can be used to solve for an effect size. The
effect size can be used to extrapolate a power curve for the data.

Most functions in this module accept a statistical test function which takes a
list of samples and returns a p value. The test is then evaluated over a series
of subsamples.

Sampling may be handled in two ways. For any set of samples, we may simply
choose to draw :math:`n` observations at random for each sample. Alternatively,
if metadata is available, samples can be matched based on a set of control
categories so that paired samples are drawn at random from the set of available
matches.

Functions
---------

.. autosummary::
    :toctree: generated/

    subsample_power
    subsample_paired_power
    confidence_bound
    paired_subsamples

Examples
--------
Suppose we wanted to test that there's a relationship between two random
variables, `ind` and `dep`. Let's use random subsampling to estimate the
statistical power of our test with an alpha of 0.1, 0.01, and 0.001.

To control for the pseudo-random number generation, we will use a seed.
When using these functions with your own data, you don't need to include the
step.

>>> import numpy as np
>>> np.random.seed(20)
>>> ind = np.random.randint(0, 20, 15)
>>> ind
array([ 3, 15,  9, 11,  7,  2,  0,  8, 19, 16,  6,  6, 16,  9,  5])
>>> dep = (3 * ind + 5 + np.random.randn(15) * 5).round(3)
>>> dep
array([ 15.617,  47.533,  28.04 ,  33.788,  19.602,  12.229,   4.779,
        36.838,  67.256,  55.032,  22.157,   7.051,  58.601,  38.664,
        18.783])

Let's define a test that will draw a list of sample pairs and determine
if they're correlated. We'll use `scipy.stats.pearsonr` which takes two arrays
and returns a correlation coefficient and a p-value representing the
probability the two distributions are correlated.

>>> from scipy.stats import pearsonr
>>> f = lambda x: pearsonr(x[0], x[1])[1]

Now, let's use random sampling to estimate the power of our test on
the first distribution.

>>> samples = [ind, dep]
>>> print("%.3e" % f(samples))
3.646e-08

In `subsample_power`, we can maintain a paired relationship between samples
by setting `draw_mode` to "matched". We can also set our critical value, so
that we estimate power for a critical value of :math:`\alpha = 0.05`, an
estimate for the critical value of 0.01, and a critical value of 0.001.

>>> from skbio.power import subsample_power
>>> pwr_100, counts_100 = subsample_power(test=f,
...                                       samples=samples,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.1,
...                                       num_iter=25)
>>> pwr_010, counts_010 = subsample_power(test=f,
...                                       samples=samples,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.01,
...                                       num_iter=25)
>>> pwr_001, counts_001 = subsample_power(test=f,
...                                       samples=samples,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.001,
...                                       num_iter=25)
>>> counts_100
array([3, 4, 5, 6, 7, 8, 9])
>>> pwr_100.mean(0)
array([ 0.484,  0.844,  0.932,  0.984,  1.   ,  1.   ,  1.   ])
>>> pwr_010.mean(0)
array([ 0.044,  0.224,  0.572,  0.836,  0.928,  0.996,  1.   ])
>>> pwr_001.mean(0)
array([ 0.   ,  0.016,  0.108,  0.332,  0.572,  0.848,  0.956])

Based on this power estimate, as we increase our confidence that we have not
committed a type I error and identified a false positive, the number of samples
we need to be confident that we have not committed a type II error increases.


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._power import (
    subsample_power,
    subsample_paired_power,
    confidence_bound,
    paired_subsamples,
)


__all__ = [
    "subsample_power",
    "subsample_paired_power",
    "confidence_bound",
    "paired_subsamples",
]
