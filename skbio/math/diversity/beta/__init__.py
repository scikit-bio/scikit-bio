"""
Beta diversity measures (:mod:`skbio.math.diversity.beta`)
==========================================================

.. currentmodule:: skbio.math.diversity.beta

This package contains helper functions for working with scipy's pairwise
distance (``pdist``) functions in scikit-bio, and will eventually be expanded
to contain pairwise distance/dissimilarity methods that are not implemented
(or planned to be implemented) in scipy.

The functions in this package currently support applying ``pdist`` functions
to all pairs of samples in a sample by observation count or abundance matrix
and returning an ``skbio.DistanceMatrix`` object. This application is
illustrated below for a few different forms of input.

Functions
---------

.. autosummary::
   :toctree: generated/

    pw_distances
    pw_distances_from_table

Examples
--------
Create a table containing 7 OTUs and 6 samples:

.. plot::
   :context:

   >>> from skbio.math.diversity.beta import pw_distances
   >>> import numpy as np
   >>> data = [[23, 64, 14, 0, 0, 3, 1],
   ...         [0, 3, 35, 42, 0, 12, 1],
   ...         [0, 5, 5, 0, 40, 40, 0],
   ...         [44, 35, 9, 0, 1, 0, 0],
   ...         [0, 2, 8, 0, 35, 45, 1],
   ...         [0, 0, 25, 35, 0, 19, 0]]
   >>> ids = list('ABCDEF')

   Compute Bray-Curtis distances between all pairs of samples and return a
   ``DistanceMatrix`` object:

   >>> bc_dm = pw_distances(data, ids, "braycurtis")
   >>> print(bc_dm)
   6x6 distance matrix
   IDs:
   A, B, C, D, E, F
   Data:
   [[ 0.          0.78787879  0.86666667  0.30927835  0.85714286  0.81521739]
    [ 0.78787879  0.          0.78142077  0.86813187  0.75        0.1627907 ]
    [ 0.86666667  0.78142077  0.          0.87709497  0.09392265  0.71597633]
    [ 0.30927835  0.86813187  0.87709497  0.          0.87777778  0.89285714]
    [ 0.85714286  0.75        0.09392265  0.87777778  0.          0.68235294]
    [ 0.81521739  0.1627907   0.71597633  0.89285714  0.68235294  0.        ]]

   Compute Jaccard distances between all pairs of samples and return a
   ``DistanceMatrix`` object:

   >>> j_dm = pw_distances(data, ids, "jaccard")
   >>> print(j_dm)
   6x6 distance matrix
   IDs:
   A, B, C, D, E, F
   Data:
   [[ 0.          0.83333333  1.          1.          0.83333333  1.        ]
    [ 0.83333333  0.          1.          1.          0.83333333  1.        ]
    [ 1.          1.          0.          1.          1.          1.        ]
    [ 1.          1.          1.          0.          1.          1.        ]
    [ 0.83333333  0.83333333  1.          1.          0.          1.        ]
    [ 1.          1.          1.          1.          1.          0.        ]]

   Determine if the resulting distance matrices are significantly correlated
   by computing the Mantel correlation between them. Then determine if the
   p-value is significant based on an alpha of 0.05:

   >>> from skbio.math.stats.distance import mantel
   >>> r, p_value = mantel(j_dm, bc_dm)
   >>> print(r)
   -0.209362157621
   >>> print(p_value < 0.05)
   False

   Compute PCoA for both distance matrices, and then find the Procrustes
   M-squared value that results from comparing the coordinate matrices:

   >>> from skbio.math.stats.ordination import PCoA
   >>> bc_pc = PCoA(bc_dm).scores()
   >>> j_pc = PCoA(j_dm).scores()
   >>> from skbio.math.stats.spatial import procrustes
   >>> print(procrustes(bc_pc.site, j_pc.site)[2])
   0.466134984787

   All of this only gets interesting in the context of sample metadata, so
   let's define some:

   >>> import pandas as pd
   >>> try:
   ...     # not necessary for normal use
   ...     pd.set_option('show_dimensions', True)
   ... except KeyError:
   ...     pass
   >>> sample_md = {
   ...    'A': {'body_site': 'gut', 'subject': '1'},
   ...    'B': {'body_site': 'skin', 'subject': '1'},
   ...    'C': {'body_site': 'tongue', 'subject': '1'},
   ...    'D': {'body_site': 'gut', 'subject': '2'},
   ...    'E': {'body_site': 'tongue', 'subject': '2'},
   ...    'F': {'body_site': 'skin', 'subject': '2'}}
   >>> sample_md = pd.DataFrame.from_dict(sample_md, orient='index')
   >>> sample_md
     subject body_site
   A       1       gut
   B       1      skin
   C       1    tongue
   D       2       gut
   E       2    tongue
   F       2      skin
   <BLANKLINE>
   [6 rows x 2 columns]

   We'll put a quick 3D plotting function together. This function is adapted
   from the matplotlib gallery [1]_.

   >>> import matplotlib.pyplot as plt
   >>> from mpl_toolkits.mplot3d import Axes3D
   >>> def scatter_3d(ord_results, df, column, color_map, title='', axis1=0,
   ...                axis2=1, axis3=2):
   ...    coord_matrix = ord_results.site.T
   ...    ids = ord_results.site_ids
   ...    colors = [color_map[df[column][id_]] for id_ in ord_results.site_ids]
   ...
   ...    fig = plt.figure()
   ...    ax = fig.add_subplot(111, projection='3d')
   ...
   ...    xs = coord_matrix[axis1]
   ...    ys = coord_matrix[axis2]
   ...    zs = coord_matrix[axis3]
   ...    plot = ax.scatter(xs, ys, zs, c=colors)
   ...
   ...    ax.set_xlabel('PC %d' % (axis1 + 1))
   ...    ax.set_ylabel('PC %d' % (axis2 + 1))
   ...    ax.set_zlabel('PC %d' % (axis3 + 1))
   ...    ax.set_xticklabels([])
   ...    ax.set_yticklabels([])
   ...    ax.set_zticklabels([])
   ...    ax.set_title(title)
   ...    return fig

   Now let's plot our PCoA results, coloring each sample by the subject it
   was taken from:

   >>> fig = scatter_3d(bc_pc, sample_md, 'subject', {'1': 'b', '2': 'r'},
   ...                  'Samples colored by subject')

.. plot::
   :context:

   We don't see any clustering/grouping of samples. If we were to instead color
   the samples by the body site they were taken from, we see that the samples
   form three separate groups:

   >>> plt.close('all') # not necessary for normal use
   >>> fig = scatter_3d(bc_pc, sample_md, 'body_site',
   ...                  {'gut': 'b', 'skin': 'r', 'tongue': 'g'},
   ...                  'Samples colored by body site')

Ordination techniques, such as PCoA, are useful for exploratory analysis. The
next step is to quantify the strength of the grouping/clustering that we see in
ordination plots. There are many statistical methods available to accomplish
this; many operate on distance matrices. Let's use ANOSIM to quantify the
strength of the clustering we see in the ordination plots above, using our
Bray-Curtis distance matrix and sample metadata.

First test the grouping of samples by subject:

>>> from skbio.math.stats.distance import ANOSIM
>>> anosim = ANOSIM(bc_dm, sample_md, column='subject')
>>> results = anosim(999)
>>> results.statistic
-0.4074074074074075
>>> results.p_value < 0.05
False

The negative value of ANOSIM's R statistic indicates anti-clustering and the
p-value is insignificant at an alpha of 0.05.

Now let's test the grouping of samples by body site:

>>> anosim = ANOSIM(bc_dm, sample_md, column='body_site')
>>> results = anosim(999)
>>> results.statistic
1.0
>>> results.p_value < 0.1
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

from .base import pw_distances, pw_distances_from_table

__all__ = ["pw_distances", "pw_distances_from_table"]

from numpy.testing import Tester
test = Tester().test
