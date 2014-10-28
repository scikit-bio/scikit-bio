# We don't see any clustering/grouping of samples. If we were to instead color
# the samples by the body site they were taken from, we see that the samples
# form three separate groups:

import matplotlib.pyplot as plt
plt.close('all') # not necessary for normal use
fig = bc_pc.plot(sample_md, 'body_site',
                 axis_labels=('PC 1', 'PC 2', 'PC 3'),
                 title='Samples colored by body site', cmap='jet', s=50)
