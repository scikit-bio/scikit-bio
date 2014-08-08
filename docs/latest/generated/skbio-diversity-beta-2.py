# We don't see any clustering/grouping of samples. If we were to instead color
# the samples by the body site they were taken from, we see that the samples
# form three separate groups:

plt.close('all') # not necessary for normal use
fig = scatter_3d(bc_pc, sample_md, 'body_site',
                 {'gut': 'b', 'skin': 'r', 'tongue': 'g'},
                 'Samples colored by body site')
