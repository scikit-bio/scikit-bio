from skbio.draw import grouped_distributions
fig = grouped_distributions('bar',
                            [[[2, 2, 1,], [0, 1, 4]],
                            [[1, 1, 1], [4, 4.5]],
                            [[2.2, 2.4, 2.7, 1.0], [0, 0.2]]],
                            distribution_labels=['Treatment 1',
                                                 'Treatment 2'])
