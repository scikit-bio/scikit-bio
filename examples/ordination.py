import os
import numpy as np

from skbio.math.stats.ordination import CA, RDA, CCA


path = os.path.dirname(os.path.abspath(__file__))


def get_path(fn):
    return os.path.join(path, os.pardir, 'math', 'stats', 'ordination',
                        'test', 'data', fn)

X = np.loadtxt(get_path('L&L_CA_data'))
ordint = CA(X)
ordint.biplot(1)
ordint.biplot(2)

Y = np.loadtxt(get_path('example2_Y'))
X = np.loadtxt(get_path('example2_X')).reshape(-1, 4, order='F')
ordint = RDA(Y, X)
ordint.biplot()

Y = np.loadtxt(get_path('example3_Y'))
X = np.loadtxt(get_path('example3_X')).reshape(-1, 4, order='F')
ordint = CCA(Y, X)
ordint.biplot()
