# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import importlib
from io import StringIO, BytesIO

from ._decorator import experimental


@experimental(as_of="0.5.10")
class PlottableMixin():
    """A plottable object.
    """

    @experimental(as_of="0.5.10")
    def _get_mpl_plt(self):
        """Import matplotlib and its plotting interface.
        """
        msg = 'Plotting requires matplotlib installed in the system.'
        if hasattr(self, 'mpl'):
            if self.mpl is None:
                raise ImportError(msg)
            return
        try:
            self.mpl = importlib.import_module('matplotlib')
        except ModuleNotFoundError:
            self.mpl = None
            raise ImportError(msg)
        else:
            self.plt = importlib.import_module('matplotlib.pyplot')

    @experimental(as_of="0.5.10")
    def _figure_data(self, format='png'):
        """Get figure data of a plottable object.

        Parameters
        ----------
        format : str, optional
            Image format supported by the plotting backend. Examples include
            'png' (default), 'svg', 'pdf', and 'eps'.

        Returns
        -------
        str or bytes
            Figure data.
        """
        self._get_mpl_plt()

        # call default plotting method
        fig = self.plot()
        fig.tight_layout()

        # get figure data
        try:
            fig.savefig(f := StringIO(), format=format)
        except TypeError:
            fig.savefig(f := BytesIO(), format=format)

        # close figure to avoid double display in IPython
        self.plt.close(fig)

        return f.getvalue()

    @experimental(as_of="0.5.10")
    def _repr_png_(self):
        """Generate a PNG format figure for display in IPython.
        """
        return self._figure_data('png')

    @experimental(as_of="0.5.10")
    def _repr_svg_(self):
        """Generate an SVG format figure for display in IPython.
        """
        return self._figure_data('svg')

    @property
    @experimental(as_of="0.4.0")
    def png(self):
        """Get figure data in PNG format.

        Returns
        -------
        bytes
            Figure data in PNG format.
        """
        return self._repr_png_()

    @property
    @experimental(as_of="0.4.0")
    def svg(self):
        """Get figure data in SVG format.

        Returns
        -------
        str
            Figure data in SVG format.
        """
        return self._repr_svg_()
