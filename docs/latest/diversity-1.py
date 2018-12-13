data = [[23, 64, 14, 0, 0, 3, 1],
        [0, 3, 35, 42, 0, 12, 1],
        [0, 5, 5, 0, 40, 40, 0],
        [44, 35, 9, 0, 1, 0, 0],
        [0, 2, 8, 0, 35, 45, 1],
        [0, 0, 25, 35, 0, 19, 0]]
ids = list('ABCDEF')

# First, we'll compute observed OTUs, an alpha diversity metric, for each
# sample using the ``alpha_diversity`` driver function:

from skbio.diversity import alpha_diversity
adiv_obs_otus = alpha_diversity('observed_otus', data, ids)
adiv_obs_otus
# A    5
# B    5
# C    4
# D    4
# E    5
# F    3
# dtype: int64

# Next we'll compute Faith's PD on the same samples. Since this is a
# phylogenetic diversity metric, we'll first create a tree and an ordered
# list of OTU identifiers.

from skbio import TreeNode
from io import StringIO
tree = TreeNode.read(StringIO(
                     '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,'
                     '(OTU4:0.75,(OTU5:0.5,(OTU6:0.5,OTU7:0.5):0.5):'
                     '0.5):1.25):0.0)root;'))
otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5', 'OTU6', 'OTU7']
adiv_faith_pd = alpha_diversity('faith_pd', data, ids=ids,
                                otu_ids=otu_ids, tree=tree)
adiv_faith_pd
# A    6.75
# B    7.00
# C    6.25
# D    5.75
# E    6.75
# F    5.50
# dtype: float64

# Now we'll compute Bray-Curtis distances, a beta diversity metric, between
# all pairs of samples. Notice that the ``data`` and ``ids`` parameters
# provided to ``beta_diversity`` are the same as those provided to
# ``alpha_diversity``.

from skbio.diversity import beta_diversity
bc_dm = beta_diversity("braycurtis", data, ids)
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

# Next, we'll compute weighted UniFrac distances between all pairs of samples.
# Because weighted UniFrac is a phylogenetic beta diversity metric, we'll need
# to pass the ``skbio.TreeNode`` and list of OTU ids that we created above.
# Again, these are the same values that were provided to ``alpha_diversity``.

wu_dm = beta_diversity("weighted_unifrac", data, ids, tree=tree,
                       otu_ids=otu_ids)
print(wu_dm)
# 6x6 distance matrix
# IDs:
# 'A', 'B', 'C', 'D', 'E', 'F'
# Data:
# [[ 0.          2.77549923  3.82857143  0.42512039  3.8547619   3.10937312]
# [ 2.77549923  0.          2.26433692  2.98435423  2.24270353  0.46774194]
# [ 3.82857143  2.26433692  0.          3.95224719  0.16025641  1.86111111]
# [ 0.42512039  2.98435423  3.95224719  0.          3.98796148  3.30870431]
# [ 3.8547619   2.24270353  0.16025641  3.98796148  0.          1.82967033]
# [ 3.10937312  0.46774194  1.86111111  3.30870431  1.82967033  0.        ]]

# Next we'll do some work with these beta diversity distance matrices. First,
# we'll determine if the UniFrac and Bray-Curtis distance matrices are
# significantly correlated by computing the Mantel correlation between them.
# Then we'll determine if the p-value is significant based on an alpha of
# 0.05.

from skbio.stats.distance import mantel
r, p_value, n = mantel(wu_dm, bc_dm)
print(r)
# 0.922404392093
alpha = 0.05
print(p_value < alpha)
# True

# Next, we'll perform principal coordinates analysis (PCoA) on our weighted
# UniFrac distance matrix.

from skbio.stats.ordination import pcoa
wu_pc = pcoa(wu_dm)

# PCoA plots are only really interesting in the context of sample metadata, so
# let's define some before we visualize these results.

import pandas as pd
sample_md = pd.DataFrame([
   ['gut', 's1'],
   ['skin', 's1'],
   ['tongue', 's1'],
   ['gut', 's2'],
   ['tongue', 's2'],
   ['skin', 's2']],
   index=['A', 'B', 'C', 'D', 'E', 'F'],
   columns=['body_site', 'subject'])
sample_md
# body_site subject
# A       gut      s1
# B      skin      s1
# C    tongue      s1
# D       gut      s2
# E    tongue      s2
# F      skin      s2

# Now let's plot our PCoA results, coloring each sample by the subject it
# was taken from:

fig = wu_pc.plot(sample_md, 'subject',
                 axis_labels=('PC 1', 'PC 2', 'PC 3'),
                 title='Samples colored by subject', cmap='jet', s=50)
