"""
Alpha diversity measures (:mod:`skbio.diversity.alpha`)
=======================================================

.. currentmodule:: skbio.diversity.alpha

This package provides implementations of various alpha diversity measures,
including measures of richness, dominance, and evenness. Some functions
generate confidence intervals (CIs). These functions have the suffix ``_ci``.
A driver function is also provided which supports computing alpha diversity for
many samples at once.

All alpha diversity measures accept a vector of counts within a single sample,
where each count is, for example, the number of observations of a particular
Operational Taxonomic Unit, or OTU. We use the term "OTU" here very loosely, as
these could be counts of any type of feature/observation (e.g., bacterial
species). We'll refer to this vector as the *counts vector* or simply *counts*
throughout the documentation. Some of these metrics incorporate phylogeny.
These metrics take a tree (``skbio.TreeNode``) and a list of OTU ids mapping
the counts to tips in the tree, in addition to the counts vector.

The counts vector must be one-dimensional and contain integers representing the
number of individuals seen (or *counted*) for a particular OTU. Negative values
are not allowed; the counts vector may only contain integers greater than or
equal to zero.

The counts vector is `array_like`: anything that can be converted into a 1-D
numpy array is acceptable input. For example, you can provide a numpy array or
a native Python list and the results should be identical.

If the input to an alpha diversity measure does not meet the above
requirements, the function will raise either a ``ValueError`` or a
``TypeError``, depending on the condition that is violated.

.. note:: There are different ways that samples are represented in the
   ecological literature and in related software. The alpha diversity measures
   provided here *always* assume that the input contains abundance data: each
   count represents the number of individuals seen for a particular OTU in the
   sample. For example, if you have two OTUs, where 3 individuals were observed
   from one of the OTUs and only a single individual was observed from the
   other, you could represent this data in the following forms (among others):

   As a vector of counts. This is the expected type of input for the alpha
   diversity measures in this module. There are 3 individuals from the OTU at
   index 0, and 1 individual from the OTU at index 1:

   >>> counts = [3, 1]

   As a vector of indices. The OTU at index 0 is observed 3 times, while the
   OTU at index 1 is observed 1 time:

   >>> indices = [0, 0, 0, 1]

   As a vector of frequencies. We have 1 OTU that is a singleton and 1 OTU that
   is a tripleton. We do not have any 0-tons or doubletons:

   >>> frequencies = [0, 1, 0, 1]

   Always use the first representation (a counts vector) with this module.

Functions
---------

.. autosummary::
   :toctree: generated/

   alpha_diversity
   ace
   berger_parker_d
   brillouin_d
   chao1
   chao1_ci
   dominance
   doubles
   enspie
   esty_ci
   faith_pd
   fisher_alpha
   gini_index
   goods_coverage
   heip_e
   kempton_taylor_q
   lladser_ci
   lladser_pe
   margalef
   mcintosh_d
   mcintosh_e
   menhinick
   michaelis_menten_fit
   observed_otus
   osd
   pielou_e
   robbins
   shannon
   simpson
   simpson_e
   singles
   strong

Examples
--------

>>> from skbio.diversity.alpha import observed_otus
>>> import numpy as np

Assume we have the following abundance data for a sample, represented as a
counts vector:

>>> counts = [1, 0, 0, 4, 1, 2, 3, 0]

We can count the number of OTUs:

>>> observed_otus(counts)
5

Note that OTUs with counts of zero are ignored.

In the previous example, we provided a Python list as input. We can also
provide other types of input that are `array_like`:

>>> observed_otus((1, 0, 0, 4, 1, 2, 3, 0)) # tuple
5
>>> observed_otus(np.array([1, 0, 0, 4, 1, 2, 3, 0])) # numpy array
5

All of the alpha diversity measures work in this manner.

Other metrics include ``singles``, which tells us how many OTUs are observed
exactly one time (i.e., are *singleton* OTUs), and ``doubles``, which tells us
how many OTUs are observed exactly two times (i.e., are *doubleton* OTUs).
Let's see how many singletons and doubletons there are in the sample:

>>> from skbio.diversity.alpha import singles, doubles
>>> singles(counts)
2
>>> doubles(counts)
1

Phylogenetic diversity metrics additionally incorporate the relative
relatedness of the OTUs in the calculation, and therefore require a tree and
a mapping of counts to OTU (tip) ids in the tree. Here we'll apply Faith's
Phylogenetic Diversity (PD) metric to the sample:

>>> from skbio import TreeNode
>>> from skbio.diversity.alpha import faith_pd
>>> from io import StringIO
>>> tree = TreeNode.read(StringIO(
...                      u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,'
...                      u'(OTU4:0.75,(OTU5:0.5,((OTU6:0.33,OTU7:0.62):0.5,'
...                      u'OTU8:0.5):0.5):0.5):1.25):0.0)root;'))
>>> otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5', 'OTU6', 'OTU7', 'OTU8']
>>> pd = faith_pd(counts, otu_ids, tree)
>>> print(round(pd, 2))
6.95

If you're calculating alpha diversity for more than one sample, you should
use the ``alpha_diversity`` function, which takes a matrix of per-sample
count vectors as input. For some metrics (notably ``faith_pd``), this will
internally use an optimized function so that calling ``alpha_diversity`` on
all of your samples will be much faster than calling the metric on each of your
samples individually.

Create a table containing 7 OTUs and 6 samples:

>>> from skbio.diversity.alpha import alpha_diversity
>>> import numpy as np
>>> data = [[23, 64, 14, 0, 0, 3, 1],
...         [0, 3, 35, 42, 0, 12, 1],
...         [0, 5, 5, 0, 40, 40, 0],
...         [44, 35, 9, 0, 1, 0, 0],
...         [0, 2, 8, 0, 35, 45, 1],
...         [0, 0, 25, 35, 0, 19, 0]]
>>> ids = list('ABCDEF')

Compute observed OTUs for each sample:

>>> alpha_diversity('observed_otus', data, ids)
A    5
B    5
C    4
D    4
E    5
F    3
dtype: int64

Next we'll compute Faith's PD on the same samples. Since this is a phylogenetic
diversity metric, we'll first need to create a tree and an ordered list of
OTU identifiers:

>>> from skbio import TreeNode
>>> from io import StringIO
>>> tree = TreeNode.read(StringIO(
...                      '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,'
...                      '(OTU4:0.75,(OTU5:0.5,(OTU6:0.5,OTU7:0.5):0.5):'
...                      '0.5):1.25):0.0)root;'))
>>> otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5', 'OTU6', 'OTU7']

Because ``faith_pd`` takes ``otu_ids`` and ``tree`` as additional parameters,
those are required to be passed as ``kwargs`` to ``alpha_diversity``.

>>> alpha_diversity('faith_pd', data, ids=ids, otu_ids=otu_ids, tree=tree)
A    6.75
B    7.00
C    6.25
D    5.75
E    6.75
F    5.50
dtype: float64

Note that we passed ``'faith_pd'`` as a string to ``alpha_diversity``. While we
could have passed the function itself (i.e., ``metric=faith_pd``) passing it as
a string results in an optmized verison of the function being used. Wherever
possible, you should pass ``metric`` as a string which will result in an
optimized version being used if available. Passing ``metric`` as a string may
not be possible if you're passing a metric that scikit-bio doesn't know about,
such as a custom one that you've developed.

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util import TestRunner

from ._ace import ace
from ._chao1 import chao1, chao1_ci
from ._base import (
    berger_parker_d, brillouin_d, dominance, doubles, enspie,
    esty_ci, faith_pd, fisher_alpha, goods_coverage, heip_e, kempton_taylor_q,
    margalef, mcintosh_d, mcintosh_e, menhinick, michaelis_menten_fit,
    observed_otus, osd, pielou_e, robbins, shannon, simpson, simpson_e,
    singles, strong, alpha_diversity)
from ._gini import gini_index
from ._lladser import lladser_pe, lladser_ci


__all__ = ['ace', 'chao1', 'chao1_ci', 'berger_parker_d',
           'brillouin_d', 'dominance', 'doubles', 'enspie', 'esty_ci',
           'faith_pd', 'fisher_alpha', 'gini_index', 'goods_coverage',
           'heip_e', 'kempton_taylor_q', 'margalef', 'mcintosh_d',
           'mcintosh_e', 'menhinick', 'michaelis_menten_fit', 'observed_otus',
           'osd', 'pielou_e', 'robbins', 'shannon', 'simpson', 'simpson_e',
           'singles', 'strong', 'lladser_pe', 'lladser_ci', 'alpha_diversity']

test = TestRunner(__file__).test
