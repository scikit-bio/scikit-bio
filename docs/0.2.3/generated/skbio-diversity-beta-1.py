from skbio.diversity.beta import pw_distances
import numpy as np
data = [[23, 64, 14, 0, 0, 3, 1],
        [0, 3, 35, 42, 0, 12, 1],
        [0, 5, 5, 0, 40, 40, 0],
        [44, 35, 9, 0, 1, 0, 0],
        [0, 2, 8, 0, 35, 45, 1],
        [0, 0, 25, 35, 0, 19, 0]]
ids = list('ABCDEF')

# Compute Bray-Curtis distances between all pairs of samples and return a
# ``DistanceMatrix`` object:

bc_dm = pw_distances(data, ids, "braycurtis")
print(bc_dm)
# 6x6 distance matrix
# IDs:
# 'A', 'B', 'C', 'D', 'E', 'F'
# Data:
# [[ 0.          0.78787879  0.86666667  0.30927835  0.85714286  0.81521739]
# [ 0.78787879  0.          0.78142077  0.86813187  0.75        0.1627907 ]
# [ 0.86666667  0.78142077  0.          0.87709497  0.09392265  0.71597633]
# [ 0.30927835  0.86813187  0.87709497  0.          0.87777778  0.89285714]
# [ 0.85714286  0.75        0.09392265  0.87777778  0.          0.68235294]
# [ 0.81521739  0.1627907   0.71597633  0.89285714  0.68235294  0.        ]]

# Compute Jaccard distances between all pairs of samples and return a
# ``DistanceMatrix`` object:

j_dm = pw_distances(data, ids, "jaccard")
print(j_dm)
# 6x6 distance matrix
# IDs:
# 'A', 'B', 'C', 'D', 'E', 'F'
# Data:
# [[ 0.          0.83333333  1.          1.          0.83333333  1.        ]
# [ 0.83333333  0.          1.          1.          0.83333333  1.        ]
# [ 1.          1.          0.          1.          1.          1.        ]
# [ 1.          1.          1.          0.          1.          1.        ]
# [ 0.83333333  0.83333333  1.          1.          0.          1.        ]
# [ 1.          1.          1.          1.          1.          0.        ]]

# Determine if the resulting distance matrices are significantly correlated
# by computing the Mantel correlation between them. Then determine if the
# p-value is significant based on an alpha of 0.05:

from skbio.stats.distance import mantel
r, p_value, n = mantel(j_dm, bc_dm)
print(r)
# -0.209362157621
print(p_value < 0.05)
# False

# Compute PCoA for both distance matrices, and then find the Procrustes
# M-squared value that results from comparing the coordinate matrices.

from skbio.stats.ordination import PCoA
bc_pc = PCoA(bc_dm).scores()
j_pc = PCoA(j_dm).scores()
from skbio.stats.spatial import procrustes
print(procrustes(bc_pc.site, j_pc.site)[2])
# 0.466134984787

# All of this only gets interesting in the context of sample metadata, so
# let's define some:

import pandas as pd
try:
    # not necessary for normal use
    pd.set_option('show_dimensions', True)
except KeyError:
    pass
sample_md = {
   'A': {'body_site': 'gut', 'subject': 's1'},
   'B': {'body_site': 'skin', 'subject': 's1'},
   'C': {'body_site': 'tongue', 'subject': 's1'},
   'D': {'body_site': 'gut', 'subject': 's2'},
   'E': {'body_site': 'tongue', 'subject': 's2'},
   'F': {'body_site': 'skin', 'subject': 's2'}}
sample_md = pd.DataFrame.from_dict(sample_md, orient='index')
sample_md
# subject body_site
# A      s1       gut
# B      s1      skin
# C      s1    tongue
# D      s2       gut
# E      s2    tongue
# F      s2      skin
# <BLANKLINE>
# [6 rows x 2 columns]

# Now let's plot our PCoA results, coloring each sample by the subject it
# was taken from:

fig = bc_pc.plot(sample_md, 'subject',
                 axis_labels=('PC 1', 'PC 2', 'PC 3'),
                 title='Samples colored by subject', cmap='jet', s=50)
