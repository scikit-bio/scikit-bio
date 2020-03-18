# Ordination techniques, such as PCoA, are useful for exploratory analysis.
# The next step is to quantify the strength of the grouping/clustering that we
# see in ordination plots. There are many statistical methods available to
# accomplish this; many operate on distance matrices. Let's use ANOSIM to
# quantify the strength of the clustering we see in the ordination plots
# above, using our weighted UniFrac distance matrix and sample metadata.

# First test the grouping of samples by subject:

from skbio.stats.distance import anosim
results = anosim(wu_dm, sample_md, column='subject', permutations=999)
results['test statistic']
# -0.33333333333333331
results['p-value'] < 0.1
# False

# The negative value of ANOSIM's R statistic indicates anti-clustering and the
# p-value is insignificant at an alpha of 0.1.

# Now let's test the grouping of samples by body site:

results = anosim(wu_dm, sample_md, column='body_site', permutations=999)
results['test statistic']
# 1.0
results['p-value'] < 0.1
# True

# The R statistic indicates strong separation of samples based on body site.
# The p-value is significant at an alpha of 0.1.

# We can also explore the alpha diversity in the context of sample metadata.
# To do this, let's add the Observed OTU and Faith PD data to our sample
# metadata. This is straight-forward beause ``alpha_diversity`` returns a
# Pandas ``Series`` object, and we're representing our sample metadata in a
# Pandas ``DataFrame`` object.

sample_md['Observed OTUs'] = adiv_obs_otus
sample_md['Faith PD'] = adiv_faith_pd
sample_md
# body_site subject  Observed OTUs  Faith PD
# A       gut      s1              5      6.75
# B      skin      s1              5      7.00
# C    tongue      s1              4      6.25
# D       gut      s2              4      5.75
# E    tongue      s2              5      6.75
# F      skin      s2              3      5.50

# We can investigate these alpha diversity data in the context of our metadata
# categories. For example, we can generate boxplots showing Faith PD by body
# site.

import matplotlib.pyplot as plt
plt.close('all') # not necessary for normal use
fig = sample_md.boxplot(column='Faith PD', by='body_site')
