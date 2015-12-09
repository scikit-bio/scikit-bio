# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import with_metaclass
from future.builtins import zip

import abc
import copy
import functools

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
from IPython.core.pylabtools import print_figure
from IPython.core.display import Image, SVG

from skbio.stats._misc import _pprint_strs
from skbio.util._decorator import stable, experimental


class SkbioObject(with_metaclass(abc.ABCMeta, object)):
    """Abstract base class defining core API common to all scikit-bio objects.

    Public scikit-bio classes should subclass this class to ensure a common,
    core API is present. All abstract methods and properties defined here must
    be implemented in subclasses, otherwise they will not be instantiable.

    """
    @abc.abstractmethod
    def __str__(self):
        pass


class MetadataMixin(with_metaclass(abc.ABCMeta, object)):
    @property
    @stable(as_of="0.4.0")
    def metadata(self):
        """``dict`` containing metadata which applies to the entire object.

        Notes
        -----
        This property can be set and deleted. When setting new metadata a
        shallow copy of the dictionary is made.

        Examples
        --------
        .. note:: scikit-bio objects with metadata share a common interface for
           accessing and manipulating their metadata. The following examples
           use scikit-bio's ``Sequence`` class to demonstrate metadata
           behavior. These examples apply to all other scikit-bio objects
           storing metadata.

        Create a sequence with metadata:

        >>> from pprint import pprint
        >>> from skbio import Sequence
        >>> seq = Sequence('ACGT', metadata={'id': 'seq-id',
        ...                                  'description': 'seq description'})

        Retrieve metadata:

        >>> pprint(seq.metadata) # using pprint to display dict in sorted order
        {'description': 'seq description', 'id': 'seq-id'}

        Update metadata:

        >>> seq.metadata['id'] = 'new-id'
        >>> seq.metadata['pubmed'] = 12345
        >>> pprint(seq.metadata)
        {'description': 'seq description', 'id': 'new-id', 'pubmed': 12345}

        Set metadata:

        >>> seq.metadata = {'abc': 123}
        >>> seq.metadata
        {'abc': 123}

        Delete metadata:

        >>> seq.has_metadata()
        True
        >>> del seq.metadata
        >>> seq.metadata
        {}
        >>> seq.has_metadata()
        False

        """
        if self._metadata is None:
            # Not using setter to avoid copy.
            self._metadata = {}
        return self._metadata

    @metadata.setter
    def metadata(self, metadata):
        if not isinstance(metadata, dict):
            raise TypeError("metadata must be a dict")
        # Shallow copy.
        self._metadata = metadata.copy()

    @metadata.deleter
    def metadata(self):
        self._metadata = None

    @abc.abstractmethod
    def __init__(self, metadata=None):
        pass

    def _init_(self, metadata=None):
        if metadata is None:
            self._metadata = None
        else:
            self.metadata = metadata

    @abc.abstractmethod
    def __eq__(self, other):
        pass

    def _eq_(self, other):
        # We're not simply comparing self.metadata to other.metadata in order
        # to avoid creating "empty" metadata representations on the objects if
        # they don't have metadata.
        if self.has_metadata() and other.has_metadata():
            if self.metadata != other.metadata:
                return False
        elif not (self.has_metadata() or other.has_metadata()):
            # Both don't have metadata.
            pass
        else:
            # One has metadata while the other does not.
            return False

        return True

    @abc.abstractmethod
    def __ne__(self, other):
        pass

    def _ne_(self, other):
        return not (self == other)

    @abc.abstractmethod
    def __copy__(self):
        pass

    def _copy_(self):
        if self.has_metadata():
            return self.metadata.copy()
        else:
            return None

    @abc.abstractmethod
    def __deepcopy__(self, memo):
        pass

    def _deepcopy_(self, memo):
        if self.has_metadata():
            return copy.deepcopy(self.metadata, memo)
        else:
            return None

    @stable(as_of="0.4.0")
    def has_metadata(self):
        """Determine if the object has metadata.

        An object has metadata if its ``metadata`` dictionary is not empty
        (i.e., has at least one key-value pair).

        Returns
        -------
        bool
            Indicates whether the object has metadata.

        Examples
        --------
        .. note:: scikit-bio objects with metadata share a common interface for
           accessing and manipulating their metadata. The following examples
           use scikit-bio's ``Sequence`` class to demonstrate metadata
           behavior. These examples apply to all other scikit-bio objects
           storing metadata.

        >>> from skbio import Sequence
        >>> seq = Sequence('ACGT')
        >>> seq.has_metadata()
        False
        >>> seq = Sequence('ACGT', metadata={})
        >>> seq.has_metadata()
        False
        >>> seq = Sequence('ACGT', metadata={'id': 'seq-id'})
        >>> seq.has_metadata()
        True

        """
        return self._metadata is not None and bool(self.metadata)


class PositionalMetadataMixin(with_metaclass(abc.ABCMeta, object)):
    @abc.abstractmethod
    def _positional_metadata_axis_len_(self):
        """Return length of axis that positional metadata applies to.

        Returns
        -------
        int
            Positional metadata axis length.

        """
        pass

    @property
    @stable(as_of="0.4.0")
    def positional_metadata(self):
        """``pd.DataFrame`` containing metadata along an axis.

        Notes
        -----
        This property can be set and deleted. When setting new positional
        metadata a shallow copy is made.

        Examples
        --------
        .. note:: scikit-bio objects with positional metadata share a common
           interface for accessing and manipulating their positional metadata.
           The following examples use scikit-bio's ``DNA`` class to demonstrate
           positional metadata behavior. These examples apply to all other
           scikit-bio objects storing positional metadata.

        Create a DNA sequence with positional metadata:

        >>> from skbio import DNA
        >>> seq = DNA(
        ...     'ACGT',
        ...     positional_metadata={'quality': [3, 3, 20, 11],
        ...                          'exons': [True, True, False, True]})
        >>> seq
        DNA
        -----------------------------
        Positional metadata:
            'exons': <dtype: bool>
            'quality': <dtype: int64>
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            GC-content: 50.00%
        -----------------------------
        0 ACGT

        Retrieve positional metadata:

        >>> seq.positional_metadata
           exons  quality
        0   True        3
        1   True        3
        2  False       20
        3   True       11

        Update positional metadata:

        >>> seq.positional_metadata['gaps'] = seq.gaps()
        >>> seq.positional_metadata
           exons  quality   gaps
        0   True        3  False
        1   True        3  False
        2  False       20  False
        3   True       11  False

        Set positional metadata:

        >>> seq.positional_metadata = {'degenerates': seq.degenerates()}
        >>> seq.positional_metadata
          degenerates
        0       False
        1       False
        2       False
        3       False

        Delete positional metadata:

        >>> seq.has_positional_metadata()
        True
        >>> del seq.positional_metadata
        >>> seq.positional_metadata
        Empty DataFrame
        Columns: []
        Index: [0, 1, 2, 3]
        >>> seq.has_positional_metadata()
        False

        """
        if self._positional_metadata is None:
            # Not using setter to avoid copy.
            self._positional_metadata = pd.DataFrame(
                index=np.arange(self._positional_metadata_axis_len_()))
        return self._positional_metadata

    @positional_metadata.setter
    def positional_metadata(self, positional_metadata):
        try:
            # Pass copy=True to copy underlying data buffer.
            positional_metadata = pd.DataFrame(positional_metadata, copy=True)
        except pd.core.common.PandasError as e:
            raise TypeError(
                "Invalid positional metadata. Must be consumable by "
                "`pd.DataFrame` constructor. Original pandas error message: "
                "\"%s\"" % e)

        num_rows = len(positional_metadata.index)
        axis_len = self._positional_metadata_axis_len_()
        if num_rows != axis_len:
            raise ValueError(
                "Number of positional metadata values (%d) must match the "
                "positional metadata axis length (%d)."
                % (num_rows, axis_len))

        positional_metadata.reset_index(drop=True, inplace=True)
        self._positional_metadata = positional_metadata

    @positional_metadata.deleter
    def positional_metadata(self):
        self._positional_metadata = None

    @abc.abstractmethod
    def __init__(self, positional_metadata=None):
        pass

    def _init_(self, positional_metadata=None):
        if positional_metadata is None:
            self._positional_metadata = None
        else:
            self.positional_metadata = positional_metadata

    @abc.abstractmethod
    def __eq__(self, other):
        pass

    def _eq_(self, other):
        # We're not simply comparing self.positional_metadata to
        # other.positional_metadata in order to avoid creating "empty"
        # positional metadata representations on the objects if they don't have
        # positional metadata.
        if self.has_positional_metadata() and other.has_positional_metadata():
            if not self.positional_metadata.equals(other.positional_metadata):
                return False
        elif not (self.has_positional_metadata() or
                  other.has_positional_metadata()):
            # Both don't have positional metadata.
            pass
        else:
            # One has positional metadata while the other does not.
            return False

        return True

    @abc.abstractmethod
    def __ne__(self, other):
        pass

    def _ne_(self, other):
        return not (self == other)

    @abc.abstractmethod
    def __copy__(self):
        pass

    def _copy_(self):
        if self.has_positional_metadata():
            # deep=True makes a shallow copy of the underlying data buffer.
            return self.positional_metadata.copy(deep=True)
        else:
            return None

    @abc.abstractmethod
    def __deepcopy__(self, memo):
        pass

    def _deepcopy_(self, memo):
        if self.has_positional_metadata():
            return copy.deepcopy(self.positional_metadata, memo)
        else:
            return None

    @stable(as_of="0.4.0")
    def has_positional_metadata(self):
        """Determine if the object has positional metadata.

        An object has positional metadata if its ``positional_metadata``
        ``pd.DataFrame`` has at least one column.

        Returns
        -------
        bool
            Indicates whether the object has positional metadata.

        Examples
        --------
        .. note:: scikit-bio objects with positional metadata share a common
           interface for accessing and manipulating their positional metadata.
           The following examples use scikit-bio's ``DNA`` class to demonstrate
           positional metadata behavior. These examples apply to all other
           scikit-bio objects storing positional metadata.

        >>> import pandas as pd
        >>> from skbio import DNA
        >>> seq = DNA('ACGT')
        >>> seq.has_positional_metadata()
        False
        >>> seq = DNA('ACGT', positional_metadata=pd.DataFrame(index=range(4)))
        >>> seq.has_positional_metadata()
        False
        >>> seq = DNA('ACGT', positional_metadata={'quality': range(4)})
        >>> seq.has_positional_metadata()
        True

        """
        return (self._positional_metadata is not None and
                len(self.positional_metadata.columns) > 0)


class OrdinationResults(SkbioObject):
    """Store ordination results, providing serialization and plotting support.

    Stores various components of ordination results. Provides methods for
    serializing/deserializing results, as well as generation of basic
    matplotlib 3-D scatterplots. Will automatically display PNG/SVG
    representations of itself within the IPython Notebook.

    Attributes
    ----------
    short_method_name : str
        Abbreviated ordination method name.
    long_method_name : str
        Ordination method name.
    eigvals : pd.Series
        The resulting eigenvalues.  The index corresponds to the ordination
        axis labels
    samples : pd.DataFrame
        The position of the samples in the ordination space, row-indexed by the
        sample id.
    features : pd.DataFrame
        The position of the features in the ordination space, row-indexed by
        the feature id.
    biplot_scores : pd.DataFrame
        Correlation coefficients of the samples with respect to the features.
    sample_constraints : pd.DataFrame
        Site constraints (linear combinations of constraining variables):
        coordinates of the sites in the space of the explanatory variables X.
        These are the fitted site scores
    proportion_explained : pd.Series
        Proportion explained by each of the dimensions in the ordination space.
        The index corresponds to the ordination axis labels
    png
    svg

    See Also
    --------
    ca
    cca
    pcoa
    rda
    """
    default_write_format = 'ordination'

    @experimental(as_of="0.4.0")
    def __init__(self, short_method_name, long_method_name, eigvals,
                 samples, features=None, biplot_scores=None,
                 sample_constraints=None, proportion_explained=None):

        self.short_method_name = short_method_name
        self.long_method_name = long_method_name

        self.eigvals = eigvals
        self.samples = samples
        self.features = features
        self.biplot_scores = biplot_scores
        self.sample_constraints = sample_constraints
        self.proportion_explained = proportion_explained

    @experimental(as_of="0.4.0")
    def __str__(self):
        """Return a string representation of the ordination results.

        String representation lists ordination results attributes and indicates
        whether or not they are present. If an attribute is present, its
        dimensions are listed. A truncated list of features and sample IDs are
        included (if they are present).

        Returns
        -------
        str
            String representation of the ordination results.

        .. shownumpydoc

        """
        lines = ['Ordination results:']
        method = '%s (%s)' % (self.long_method_name, self.short_method_name)
        lines.append(self._format_attribute(method, 'Method', str))

        attrs = [(self.eigvals, 'Eigvals'),
                 (self.proportion_explained, 'Proportion explained'),
                 (self.features, 'Features'),
                 (self.samples, 'Samples'),
                 (self.biplot_scores, 'Biplot Scores'),
                 (self.sample_constraints, 'Sample constraints')]
        for attr, attr_label in attrs:
            def formatter(e):
                return 'x'.join(['%d' % s for s in e.shape])

            lines.append(self._format_attribute(attr, attr_label, formatter))

        lines.append(self._format_attribute(
            self.features, 'Feature IDs',
            lambda e: _pprint_strs(e.index.tolist())))
        lines.append(self._format_attribute(
            self.samples, 'Sample IDs',
            lambda e: _pprint_strs(e.index.tolist())))

        return '\n'.join(lines)

    @experimental(as_of="0.4.0")
    def plot(self, df=None, column=None, axes=(0, 1, 2), axis_labels=None,
             title='', cmap=None, s=20):
        """Create a 3-D scatterplot of ordination results colored by metadata.

        Creates a 3-D scatterplot of the ordination results, where each point
        represents a sample. Optionally, these points can be colored by
        metadata (see `df` and `column` below).

        Parameters
        ----------
        df : pd.DataFrame, optional
            ``DataFrame`` containing sample metadata. Must be indexed by sample
            ID, and all sample IDs in the ordination results must exist in the
            ``DataFrame``. If ``None``, samples (i.e., points) will not be
            colored by metadata.
        column : str, optional
            Column name in `df` to color samples (i.e., points in the plot) by.
            Cannot have missing data (i.e., ``np.nan``). `column` can be
            numeric or categorical. If numeric, all values in the column will
            be cast to ``float`` and mapped to colors using `cmap`. A colorbar
            will be included to serve as a legend. If categorical (i.e., not
            all values in `column` could be cast to ``float``), colors will be
            chosen for each category using evenly-spaced points along `cmap`. A
            legend will be included. If ``None``, samples (i.e., points) will
            not be colored by metadata.
        axes : iterable of int, optional
            Indices of sample coordinates to plot on the x-, y-, and z-axes.
            For example, if plotting PCoA results, ``axes=(0, 1, 2)`` will plot
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
            - sample IDs in the ordination results are not in `df` or have
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

           Define a distance matrix with four samples labelled A-D:

           >>> from skbio import DistanceMatrix
           >>> dm = DistanceMatrix([[0., 0.21712454, 0.5007512, 0.91769271],
           ...                      [0.21712454, 0., 0.45995501, 0.80332382],
           ...                      [0.5007512, 0.45995501, 0., 0.65463348],
           ...                      [0.91769271, 0.80332382, 0.65463348, 0.]],
           ...                     ['A', 'B', 'C', 'D'])

           Define metadata for each sample in a ``pandas.DataFrame``:

           >>> import pandas as pd
           >>> metadata = {
           ...     'A': {'body_site': 'skin'},
           ...     'B': {'body_site': 'gut'},
           ...     'C': {'body_site': 'gut'},
           ...     'D': {'body_site': 'skin'}}
           >>> df = pd.DataFrame.from_dict(metadata, orient='index')

           Run principal coordinate analysis (PCoA) on the distance matrix:

           >>> from skbio.stats.ordination import pcoa
           >>> pcoa_results = pcoa(dm)

           Plot the ordination results, where each sample is colored by body
           site (a categorical variable):

           >>> fig = pcoa_results.plot(df=df, column='body_site',
           ...                         title='Samples colored by body site',
           ...                         cmap='Set1', s=50)

        """
        # Note: New features should not be added to this method and should
        # instead be added to EMPeror (http://biocore.github.io/emperor/).
        # Only bug fixes and minor updates should be made to this method.

        coord_matrix = self.samples.values.T
        self._validate_plot_axes(coord_matrix, axes)

        # derived from
        # http://matplotlib.org/examples/mplot3d/scatter3d_demo.html
        fig = plt.figure()
        # create the axes, leaving room for a legend as described here:
        # http://stackoverflow.com/a/9651897/3424666
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75], projection='3d')

        xs = coord_matrix[axes[0]]
        ys = coord_matrix[axes[1]]
        zs = coord_matrix[axes[2]]

        point_colors, category_to_color = self._get_plot_point_colors(
            df, column, self.samples.index, cmap)

        scatter_fn = functools.partial(ax.scatter, xs, ys, zs, s=s)
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
            raise ValueError("`axes` must contain exactly three elements "
                             "(found %d elements)." % len(axes))
        if len(set(axes)) != 3:
            raise ValueError("The values provided for `axes` must be unique.")

        for idx, axis in enumerate(axes):
            if axis < 0 or axis >= num_dims:
                raise ValueError("`axes[%d]` must be >= 0 and < %d." %
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
    @experimental(as_of="0.4.0")
    def png(self):
        """Display basic 3-D scatterplot in IPython Notebook as PNG."""
        return Image(self._repr_png_(), embed=True)

    @property
    @experimental(as_of="0.4.0")
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

    def _format_attribute(self, attr, attr_label, formatter):
        if attr is None:
            formatted_attr = 'N/A'
        else:
            formatted_attr = formatter(attr)
        return '\t%s: %s' % (attr_label, formatted_attr)


class ElasticLines(object):
    """Store blocks of content separated by dashed lines.

    Each dashed line (separator) is as long as the longest content
    (non-separator) line.

    """

    def __init__(self):
        self._lines = []
        self._separator_idxs = []
        self._max_line_len = -1

    def add_line(self, line):
        line_len = len(line)
        if line_len > self._max_line_len:
            self._max_line_len = line_len
        self._lines.append(line)

    def add_lines(self, lines):
        for line in lines:
            self.add_line(line)

    def add_separator(self):
        self._lines.append(None)
        self._separator_idxs.append(len(self._lines) - 1)

    def to_str(self):
        separator = '-' * self._max_line_len
        for idx in self._separator_idxs:
            self._lines[idx] = separator
        return '\n'.join(self._lines)
