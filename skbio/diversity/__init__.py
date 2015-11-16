"""
Diversity calculations (:mod:`skbio.diversity`)
===============================================

.. currentmodule:: skbio.diversity

This package provides functionality for calculating community diversity,
including various alpha- and beta-diversity measures.

Subpackages
-----------

.. autosummary::
   :toctree: generated/

   alpha
   beta

Functions
---------

.. autosummary::
   :toctree: generated/

    alpha_diversity
    beta_diversity

Examples
--------

Create a table containing 7 OTUs and 6 samples:

.. plot::
   :context:

>>> import numpy as np
>>> data = [[23, 64, 14, 0, 0, 3, 1],
...         [0, 3, 35, 42, 0, 12, 1],
...         [0, 5, 5, 0, 40, 40, 0],
...         [44, 35, 9, 0, 1, 0, 0],
...         [0, 2, 8, 0, 35, 45, 1],
...         [0, 0, 25, 35, 0, 19, 0]]
>>> ids = list('ABCDEF')

Compute observed OTUs for each sample:

>>> from skbio.diversity import alpha_diversity
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


The value that you provide for to ``beta_diversity`` for ``metric`` can
be either a string (e.g., "unweighted_unifrac") or a function
(e.g., ``skbio.diversity.beta.unweighted_unifrac``). The metric should
generally be passed as a string, as this often uses an optimized version
of the metric. For example, passing  ``"unweighted_unifrac"`` (a string)
will be hundreds of times faster than passing the function
``skbio.diversity.beta.unweighted_unifrac``. The latter is faster if
computing only one or a few distances, but in these cases the difference in
runtime is negligible, so it's safer to just err on the side of passing
``metric`` as a string.

Compute Bray-Curtis distances between all pairs of samples and return a
``DistanceMatrix`` object.

>>> from skbio.diversity import beta_diversity
>>> bc_dm = beta_diversity("braycurtis", data, ids)
>>> print(bc_dm)
6x6 distance matrix
IDs:
'A', 'B', 'C', 'D', 'E', 'F'
Data:
[[ 0.          0.78787879  0.86666667  0.30927835  0.85714286  0.81521739]
 [ 0.78787879  0.          0.78142077  0.86813187  0.75        0.1627907 ]
 [ 0.86666667  0.78142077  0.          0.87709497  0.09392265  0.71597633]
 [ 0.30927835  0.86813187  0.87709497  0.          0.87777778  0.89285714]
 [ 0.85714286  0.75        0.09392265  0.87777778  0.          0.68235294]
 [ 0.81521739  0.1627907   0.71597633  0.89285714  0.68235294  0.        ]]

Compute weighted UniFrac distances between all pairs of samples and return a
``DistanceMatrix`` object. Because weighted UniFrac is a phylogenetic beta
diversity metric, we'll need to create a ``skbio.TreeNode`` object that
contains all of the tips in the tree, and pass that along with the ids of
the OTUs corresponding to the counts in ``data``.

>>> from skbio import TreeNode
>>> from io import StringIO
>>> tree = TreeNode.read(StringIO(
...                      '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,'
...                      '(OTU4:0.75,(OTU5:0.5,(OTU6:0.5,OTU7:0.5):0.5):0.5'
...                      '):1.25):0.0)root;'))
>>> otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5', 'OTU6', 'OTU7']
>>> wu_dm = beta_diversity("weighted_unifrac", data, ids, tree=tree,
...                      otu_ids=otu_ids)
>>> print(wu_dm)
6x6 distance matrix
IDs:
'A', 'B', 'C', 'D', 'E', 'F'
Data:
[[ 0.          2.77549923  3.82857143  0.42512039  3.8547619   3.10937312]
 [ 2.77549923  0.          2.26433692  2.98435423  2.24270353  0.46774194]
 [ 3.82857143  2.26433692  0.          3.95224719  0.16025641  1.86111111]
 [ 0.42512039  2.98435423  3.95224719  0.          3.98796148  3.30870431]
 [ 3.8547619   2.24270353  0.16025641  3.98796148  0.          1.82967033]
 [ 3.10937312  0.46774194  1.86111111  3.30870431  1.82967033  0.        ]]

Determine if the resulting distance matrices are significantly correlated
by computing the Mantel correlation between them. Then determine if the
p-value is significant based on an alpha of 0.05:

>>> from skbio.stats.distance import mantel
>>> r, p_value, n = mantel(wu_dm, bc_dm)
>>> print(r)
0.922404392093
>>> print(p_value < 0.05)
True

Compute PCoA for both distance matrices, and then find the Procrustes
M-squared value that results from comparing the coordinate matrices.

>>> from skbio.stats.ordination import pcoa
>>> bc_pc = pcoa(bc_dm)
>>> wu_pc = pcoa(wu_dm)
>>> from skbio.stats.spatial import procrustes
>>> print(procrustes(bc_pc.samples.values, wu_pc.samples.values)[2])
0.096574934963

All of this only gets interesting in the context of sample metadata, so
let's define some:

>>> import pandas as pd
>>> sample_md = [
...    ('A', ['gut', 's1']),
...    ('B', ['skin', 's1']),
...    ('C', ['tongue', 's1']),
...    ('D', ['gut', 's2']),
...    ('E', ['tongue', 's2']),
...    ('F', ['skin', 's2'])]
>>> sample_md = pd.DataFrame.from_items(
...     sample_md, columns=['body_site', 'subject'], orient='index')
>>> sample_md
  body_site subject
A       gut      s1
B      skin      s1
C    tongue      s1
D       gut      s2
E    tongue      s2
F      skin      s2

Now let's plot our PCoA results, coloring each sample by the subject it
was taken from:

>>> fig = wu_pc.plot(sample_md, 'subject',
...                  axis_labels=('PC 1', 'PC 2', 'PC 3'),
...                  title='Samples colored by subject', cmap='jet', s=50)

.. plot::
:context:

We don't see any clustering/grouping of samples. If we were to instead color
the samples by the body site they were taken from, we see that the samples
form three separate groups:

>>> import matplotlib.pyplot as plt
>>> plt.close('all') # not necessary for normal use
>>> fig = wu_pc.plot(sample_md, 'body_site',
...                  axis_labels=('PC 1', 'PC 2', 'PC 3'),
...                  title='Samples colored by body site', cmap='jet', s=50)

Ordination techniques, such as PCoA, are useful for exploratory analysis. The
next step is to quantify the strength of the grouping/clustering that we see in
ordination plots. There are many statistical methods available to accomplish
this; many operate on distance matrices. Let's use ANOSIM to quantify the
strength of the clustering we see in the ordination plots above, using our
weighted UniFrac distance matrix and sample metadata.

First test the grouping of samples by subject:

>>> from skbio.stats.distance import anosim
>>> results = anosim(wu_dm, sample_md, column='subject', permutations=999)
>>> results['test statistic']
-0.33333333333333331
>>> results['p-value'] < 0.1
False

The negative value of ANOSIM's R statistic indicates anti-clustering and the
p-value is insignificant at an alpha of 0.1.

Now let's test the grouping of samples by body site:

>>> results = anosim(wu_dm, sample_md, column='body_site', permutations=999)
>>> results['test statistic']
1.0
>>> results['p-value'] < 0.1
True

The R statistic of 1.0 indicates strong separation of samples based on body
site. The p-value is significant at an alpha of 0.1.

References
----------
.. [1] http://matplotlib.org/examples/mplot3d/scatter3d_demo.html



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

from .alpha._base import alpha_diversity
from .beta._base import beta_diversity

__all__ = ["alpha_diversity", "beta_diversity"]

test = TestRunner(__file__).test
