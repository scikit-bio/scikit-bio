from skbio.draw.distributions import boxplots
fig = boxplots(
    [[2, 2, 1, 3], [0, -1, 0, 0.1, 0.3], [4, 5, 6, 3]],
    x_tick_labels=('Control', 'Treatment 1', 'Treatment 2'),
    box_colors=('green', 'blue', 'red'))
