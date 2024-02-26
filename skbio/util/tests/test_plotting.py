# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import unittest

try:
    import matplotlib
except ImportError:
    has_matplotlib = False
else:
    has_matplotlib = True

from skbio.util._plotting import PlottableMixin


@unittest.skipUnless(has_matplotlib, "Matplotlib not available.")
class TestPlottableMixin(unittest.TestCase):
    def setUp(self):

        def _plot(self):
            fig, ax = self.plt.subplots()
            ax.plot(1, 1, color='k', marker='o')
            return fig

        PlottableMixin.plot = _plot

    def test_get_mpl_plt(self):
        obj = PlottableMixin()

        # hasn't imported yet
        self.assertFalse(hasattr(obj, 'mpl'))

        # import Matplotlib if available
        obj._get_mpl_plt()
        self.assertEqual(obj.mpl.__name__, 'matplotlib')
        self.assertEqual(obj.plt.__name__, 'matplotlib.pyplot')

        # make Matplotlib unimportable
        delattr(obj, 'mpl')
        backup = sys.modules['matplotlib']
        sys.modules['matplotlib'] = None
        with self.assertRaises(ImportError):
            obj._get_mpl_plt()

        # won't try again if already failed
        sys.modules['matplotlib'] = backup
        self.assertIsNone(obj.mpl)
        with self.assertRaises(ImportError):
            obj._get_mpl_plt()

    def test_figure_data(self):
        obj = PlottableMixin()
        obj._get_mpl_plt()

        # PNG data are bytes
        obs = obj._figure_data('png')
        self.assertIsInstance(obs, bytes)
        self.assertTrue(len(obs) > 0)

        # SVG data are string
        obs = obj._figure_data('svg')
        self.assertIsInstance(obs, str)
        self.assertTrue(len(obs) > 0)

        # plotting backend not available
        delattr(obj, 'mpl')
        backup = sys.modules['matplotlib']
        sys.modules['matplotlib'] = None
        self.assertIsNone(obj._figure_data())
        sys.modules['matplotlib'] = backup

    def test_repr_png(self):
        obj = PlottableMixin()
        obj._get_mpl_plt()
        obs = obj._repr_png_()
        self.assertIsInstance(obs, bytes)
        self.assertTrue(len(obs) > 0)

    def test_repr_svg(self):
        obj = PlottableMixin()
        obj._get_mpl_plt()
        obs = obj._repr_svg_()
        self.assertIsInstance(obs, str)
        self.assertTrue(len(obs) > 0)

    def test_png(self):
        obj = PlottableMixin()
        obj._get_mpl_plt()
        obs = obj.png
        self.assertIsInstance(obs, bytes)
        self.assertTrue(len(obs) > 0)

    def test_svg(self):
        obj = PlottableMixin()
        obj._get_mpl_plt()
        obs = obj.svg
        self.assertIsInstance(obs, str)
        self.assertTrue(len(obs) > 0)


if __name__ == '__main__':
    unittest.main()
