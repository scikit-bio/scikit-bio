# Define a distance matrix with four samples labelled A-D:

from skbio import DistanceMatrix
dm = DistanceMatrix([[0., 0.21712454, 0.5007512, 0.91769271],
                     [0.21712454, 0., 0.45995501, 0.80332382],
                     [0.5007512, 0.45995501, 0., 0.65463348],
                     [0.91769271, 0.80332382, 0.65463348, 0.]],
                    ['A', 'B', 'C', 'D'])

# Define metadata for each sample in a ``pandas.DataFrame``:

import pandas as pd
metadata = {
    'A': {'body_site': 'skin'},
    'B': {'body_site': 'gut'},
    'C': {'body_site': 'gut'},
    'D': {'body_site': 'skin'}}
df = pd.DataFrame.from_dict(metadata, orient='index')

# Run principal coordinate analysis (PCoA) on the distance matrix:

from skbio.stats.ordination import pcoa
pcoa_results = pcoa(dm)

# Plot the ordination results, where each sample is colored by body
# site (a categorical variable):

fig = pcoa_results.plot(df=df, column='body_site',
                        title='Samples colored by body site',
                        cmap='Set1', s=50)
