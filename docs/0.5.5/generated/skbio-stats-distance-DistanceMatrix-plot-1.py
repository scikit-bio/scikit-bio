# Define a dissimilarity matrix with five objects labeled A-E:

from skbio.stats.distance import DissimilarityMatrix
dm = DissimilarityMatrix([[0, 1, 2, 3, 4], [1, 0, 1, 2, 3],
                          [2, 1, 0, 1, 2], [3, 2, 1, 0, 1],
                          [4, 3, 2, 1, 0]],
                         ['A', 'B', 'C', 'D', 'E'])

# Plot the dissimilarity matrix as a heatmap:

fig = dm.plot(cmap='Reds', title='Example heatmap')
