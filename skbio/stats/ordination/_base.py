# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip

import warnings
from functools import partial
from importlib import import_module

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# avoid flake8 unused import error
Axes3D
from IPython.core.pylabtools import print_figure
from IPython.display import Image, SVG

# This will be the responsibility of the ABC in the future.
import_module('skbio.io')


class OrdinationResults(object):
    """Store ordination results, providing serialization and plotting support.

    Stores various components of ordination results. Provides methods for
    serializing/deserializing results, as well as generation of basic
    matplotlib 3-D scatterplots. Will automatically display PNG/SVG
    representations of itself within the IPython Notebook.

    Attributes
    ----------
    eigvals : 1-D numpy array
        The result eigenvalues
    species : 2-D numpy array
        The result coordinates for each species
    site : 2-D numpy array
        The results coordinates for each site
    biplot : 2-D numpy array
        The result biplot coordinates
    site_constraints : 2-D numpy array
        The result coordinates for each site constraint
    proportion_explained : 1-D numpy array
        The proportion explained by each eigenvector
    species_ids : list of str
        The species identifiers
    site_ids : list of str
        The site identifiers
    png
    svg

    """
    default_write_format = 'ordres'

    def __init__(self, eigvals, species=None, site=None, biplot=None,
                 site_constraints=None, proportion_explained=None,
                 species_ids=None, site_ids=None):
        self.eigvals = eigvals
        self.species = species
        self.site = site
        self.biplot = biplot
        self.site_constraints = site_constraints
        self.proportion_explained = proportion_explained
        self.species_ids = species_ids
        self.site_ids = site_ids

    @classmethod
    def from_file(cls, ord_res_f):
        """Load ordination results from text file.

        .. note:: Deprecated in scikit-bio 0.2.0-dev
           ``from_file`` will be removed in scikit-bio 0.3.0. It is replaced by
           ``read``, which is a more general method for deserializing
           ordination results. ``read`` supports multiple file formats,
           automatic file format detection, etc. by taking advantage of
           scikit-bio's I/O registry system. See :mod:`skbio.io` for more
           details.

        Creates an ``OrdinationResults`` instance from a ``ordres`` formatted
        file. See :mod:`skbio.io.ordres` for the format specification.

        Parameters
        ----------
        ord_res_f: filepath or filehandle
            File to read from.

        Returns
        -------
        OrdinationResults
            Instance of type `cls` containing the parsed contents of
            `ord_res_f`.

        Raises
        ------
        OrdResFormatError
            If the format of the file is not valid, or if the shapes of the
            different sections of the file are not consistent.

        See Also
        --------
        read

        """
        warnings.warn(
            "OrdinationResults.from_file is deprecated and will be removed in "
            "scikit-bio 0.3.0. Please update your code to use "
            "OrdinationResults.read.", UserWarning)
        return cls.read(ord_res_f, format='ordres')

    def to_file(self, out_f):
        """Save ordination results to file in text format.

        .. note:: Deprecated in scikit-bio 0.2.0-dev
           ``to_file`` will be removed in scikit-bio 0.3.0. It is replaced by
           ``write``, which is a more general method for serializing ordination
           results. ``write`` supports multiple file formats by taking
           advantage of scikit-bio's I/O registry system. See :mod:`skbio.io`
           for more details.

        Serializes ordination results as an ``ordres`` formatted file. See
        :mod:`skbio.io.ordres` for the format specification.

        Parameters
        ----------
        out_f : filepath or filehandle
            File to write to.

        See Also
        --------
        write

        """
        warnings.warn(
            "OrdinationResults.to_file is deprecated and will be removed in "
            "scikit-bio 0.3.0. Please update your code to use "
            "OrdinationResults.write.", UserWarning)
        self.write(out_f, format='ordres')

    def plot(self, df=None, column=None, axes=(0, 1, 2), axis_labels=None,
             title='', cmap=None, s=20):
        """Create a 3-D scatterplot of ordination results colored by metadata.

        Creates a 3-D scatterplot of the ordination results, where each point
        represents a site. Optionally, these points can be colored by metadata
        (see `df` and `column` below).

        Parameters
        ----------
        df : pandas.DataFrame, optional
            ``DataFrame`` containing site metadata. Must be indexed by site ID,
            and all site IDs in the ordination results must exist in the
            ``DataFrame``. If ``None``, sites (i.e., points) will not be
            colored by metadata.
        column : str, optional
            Column name in `df` to color sites (i.e., points in the plot) by.
            Cannot have missing data (i.e., ``np.nan``). `column` can be
            numeric or categorical. If numeric, all values in the column will
            be cast to ``float`` and mapped to colors using `cmap`. A colorbar
            will be included to serve as a legend. If categorical (i.e., not
            all values in `column` could be cast to ``float``), colors will be
            chosen for each category using evenly-spaced points along `cmap`. A
            legend will be included. If ``None``, sites (i.e., points) will not
            be colored by metadata.
        axes : iterable of int, optional
            Indices of site coordinates to plot on the x-, y-, and z-axes. For
            example, if plotting PCoA results, ``axes=(0, 1, 2)`` will plot
            PC 1 on the x-axis, PC 2 on the y-axis, and PC 3 on the z-axis.
            Must contain exactly three elements.
        axis_labels : iterable of str, optional
            Labels for the x-, y-, and z-axes. If ``None``, labels will be the
            values of `axes` cast as strings.
        title : str, optional
            Plot title.
        cmap : str or matplotlib.colors.Colormap, optional
            Name or instance of matplotlib colormap to use for mapping `column`
            values to colors. If ``None``, defaults to the colormap specified
            in the matplotlib rc file. Qualitative colormaps (e.g., ``Set1``)
            are recommended for categorical data, while sequential colormaps
            (e.g., ``Greys``) are recommended for numeric data. See [1]_ for
            these colormap classifications.
        s : scalar or iterable of scalars, optional
            Size of points. See matplotlib's ``Axes3D.scatter`` documentation
            for more details.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing the scatterplot and legend/colorbar if metadata
            were provided.

        Raises
        ------
        ValueError
            Raised on invalid input, including the following situations:

            - there are not at least three dimensions to plot
            - there are not exactly three values in `axes`, they are not
              unique, or are out of range
            - there are not exactly three values in `axis_labels`
            - either `df` or `column` is provided without the other
            - `column` is not in the ``DataFrame``
            - site IDs in the ordination results are not in `df` or have
              missing data in `column`

        See Also
        --------
        mpl_toolkits.mplot3d.Axes3D.scatter

        Notes
        -----
        This method creates basic plots of ordination results, and is intended
        to provide a quick look at the results in the context of metadata
        (e.g., from within the IPython Notebook). For more customization and to
        generate publication-quality figures, we recommend EMPeror [2]_.

        References
        ----------
        .. [1] http://matplotlib.org/examples/color/colormaps_reference.html
        .. [2] EMPeror: a tool for visualizing high-throughput microbial
           community data. Vazquez-Baeza Y, Pirrung M, Gonzalez A, Knight R.
           Gigascience. 2013 Nov 26;2(1):16. http://biocore.github.io/emperor/

        Examples
        --------
        .. plot::

           Define a distance matrix with four sites labelled A-D:

           >>> from skbio import DistanceMatrix
           >>> dm = DistanceMatrix([[0., 0.21712454, 0.5007512, 0.91769271],
           ...                      [0.21712454, 0., 0.45995501, 0.80332382],
           ...                      [0.5007512, 0.45995501, 0., 0.65463348],
           ...                      [0.91769271, 0.80332382, 0.65463348, 0.]],
           ...                     ['A', 'B', 'C', 'D'])

           Define metadata for each site in a ``pandas.DataFrame``:

           >>> import pandas as pd
           >>> metadata = {
           ...     'A': {'body_site': 'skin'},
           ...     'B': {'body_site': 'gut'},
           ...     'C': {'body_site': 'gut'},
           ...     'D': {'body_site': 'skin'}}
           >>> df = pd.DataFrame.from_dict(metadata, orient='index')

           Run principal coordinate analysis (PCoA) on the distance matrix:

           >>> from skbio.stats.ordination import PCoA
           >>> pcoa_results = PCoA(dm).scores()

           Plot the ordination results, where each site is colored by body site
           (a categorical variable):

           >>> fig = pcoa_results.plot(df=df, column='body_site',
           ...                         title='Sites colored by body site',
           ...                         cmap='Set1', s=50)

        """
        # Note: New features should not be added to this method and should
        # instead be added to EMPeror (http://biocore.github.io/emperor/).
        # Only bug fixes and minor updates should be made to this method.

        coord_matrix = self.site.T
        self._validate_plot_axes(coord_matrix, axes)

        # derived from
        # http://matplotlib.org/examples/mplot3d/scatter3d_demo.html
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        xs = coord_matrix[axes[0]]
        ys = coord_matrix[axes[1]]
        zs = coord_matrix[axes[2]]

        point_colors, category_to_color = self._get_plot_point_colors(
            df, column, self.site_ids, cmap)

        scatter_fn = partial(ax.scatter, xs, ys, zs, s=s)
        if point_colors is None:
            plot = scatter_fn()
        else:
            plot = scatter_fn(c=point_colors, cmap=cmap)

        if axis_labels is None:
            axis_labels = ['%d' % axis for axis in axes]
        elif len(axis_labels) != 3:
            raise ValueError("axis_labels must contain exactly three elements "
                             "(found %d elements)." % len(axis_labels))

        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        ax.set_zlabel(axis_labels[2])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.set_title(title)

        # create legend/colorbar
        if point_colors is not None:
            if category_to_color is None:
                fig.colorbar(plot)
            else:
                self._plot_categorical_legend(ax, category_to_color)

        return fig

    def _validate_plot_axes(self, coord_matrix, axes):
        """Validate `axes` against coordinates matrix."""
        num_dims = coord_matrix.shape[0]
        if num_dims < 3:
            raise ValueError("At least three dimensions are required to plot "
                             "ordination results. There are only %d "
                             "dimension(s)." % num_dims)
        if len(axes) != 3:
            raise ValueError("axes must contain exactly three elements (found "
                             "%d elements)." % len(axes))
        if len(set(axes)) != 3:
            raise ValueError("The values provided for axes must be unique.")

        for idx, axis in enumerate(axes):
            if axis < 0 or axis >= num_dims:
                raise ValueError("axes[%d] must be >= 0 and < %d." %
                                 (idx, num_dims))

    def _get_plot_point_colors(self, df, column, ids, cmap):
        """Return a list of colors for each plot point given a metadata column.

        If `column` is categorical, additionally returns a dictionary mapping
        each category (str) to color (used for legend creation).

        """
        if ((df is None and column is not None) or (df is not None and
                                                    column is None)):
            raise ValueError("Both df and column must be provided, or both "
                             "must be None.")
        elif df is None and column is None:
            point_colors, category_to_color = None, None
        else:
            if column not in df:
                raise ValueError("Column '%s' not in data frame." % column)

            col_vals = df.loc[ids, column]

            if col_vals.isnull().any():
                raise ValueError("One or more IDs in the ordination results "
                                 "are not in the data frame, or there is "
                                 "missing data in the data frame's '%s' "
                                 "column." % column)

            category_to_color = None
            try:
                point_colors = col_vals.astype(float)
            except ValueError:
                # we have categorical data, so choose a color for each
                # category, where colors are evenly spaced across the
                # colormap.
                # derived from http://stackoverflow.com/a/14887119
                categories = col_vals.unique()
                cmap = plt.get_cmap(cmap)
                category_colors = cmap(np.linspace(0, 1, len(categories)))

                category_to_color = dict(zip(categories, category_colors))
                point_colors = col_vals.apply(lambda x: category_to_color[x])

            point_colors = point_colors.tolist()

        return point_colors, category_to_color

    def _plot_categorical_legend(self, ax, color_dict):
        """Add legend to plot using specified mapping of category to color."""
        # derived from http://stackoverflow.com/a/20505720
        proxies = []
        labels = []
        for category in color_dict:
            proxy = mpl.lines.Line2D([0], [0], linestyle='none',
                                     c=color_dict[category], marker='o')
            proxies.append(proxy)
            labels.append(category)

        # place legend outside of the axes (centered)
        # derived from http://matplotlib.org/users/legend_guide.html
        ax.legend(proxies, labels, numpoints=1, loc=6,
                  bbox_to_anchor=(1.05, 0.5), borderaxespad=0.)

    # Here we define the special repr methods that provide the IPython display
    # protocol. Code derived from:
    #     https://github.com/ipython/ipython/blob/2.x/examples/Notebook/
    #         Custom%20Display%20Logic.ipynb
    # See licenses/ipython.txt for more details.

    def _repr_png_(self):
        return self._figure_data('png')

    def _repr_svg_(self):
        return self._figure_data('svg')

    # We expose the above reprs as properties, so that the user can see them
    # directly (since otherwise the client dictates which one it shows by
    # default)
    @property
    def png(self):
        """Display basic 3-D scatterplot in IPython Notebook as PNG."""
        return Image(self._repr_png_(), embed=True)

    @property
    def svg(self):
        """Display basic 3-D scatterplot in IPython Notebook as SVG."""
        return SVG(self._repr_svg_())

    def _figure_data(self, format):
        fig = self.plot()
        data = print_figure(fig, format)
        # We MUST close the figure, otherwise IPython's display machinery
        # will pick it up and send it as output, resulting in a double display
        plt.close(fig)
        return data


class Ordination(object):
    short_method_name = 'Overwrite in subclass!'
    long_method_name = 'Overwrite in subclass!'
