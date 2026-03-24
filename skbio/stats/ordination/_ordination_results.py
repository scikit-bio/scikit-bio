# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import functools
from warnings import warn

import numpy as np
import pandas as pd

from skbio._base import SkbioObject
from skbio.stats._misc import _pprint_strs
from skbio.util._plotting import PlottableMixin
from skbio.io.descriptors import Read, Write
from skbio.table._tabular import _extract_row_ids


class OrdinationResults(SkbioObject, PlottableMixin):
    """Store ordination results, providing serialization and plotting support.

    Stores various components of ordination results. Provides methods for
    serializing/deserializing results, as well as generation of basic
    matplotlib 3-D scatterplots using the
    :meth:`plot` method.

    Attributes
    ----------
    short_method_name : str
        Abbreviated ordination method name.
    long_method_name : str
        Ordination method name.
    eigvals : table_like
        The resulting eigenvalues. The index corresponds to the ordination
        axis labels. See :ref:`table_output` for details.
    samples : table_like
        The position of the samples in the ordination space, row-indexed by the
        sample id. See :ref:`table_output` for details.
    features : table_like
        The position of the features in the ordination space, row-indexed by
        the feature id. See :ref:`table_output` for details.
    biplot_scores : table_like
        Correlation coefficients of the samples with respect to the features.
        See :ref:`table_output` for details.
    sample_constraints : table_like
        Site constraints (linear combinations of constraining variables):
        coordinates of the sites in the space of the explanatory variables X.
        These are the fitted site scores. See :ref:`table_output` for details.
    proportion_explained : table_like
        Proportion explained by each of the dimensions in the ordination space.
        The index corresponds to the ordination axis labels. See
        :ref:`table_output` for details.
    sample_ids, feature_ids, constraint_ids, output_format : optional
        Standard table parameters. See :ref:`table_params` for details.

    See Also
    --------
    ca
    cca
    pcoa
    rda

    """

    default_write_format = "ordination"
    """Default write format for this object: ``ordination``."""

    read = Read()
    write = Write()

    def __init__(
        self,
        short_method_name,
        long_method_name,
        eigvals,
        samples,
        sample_ids=None,
        features=None,
        feature_ids=None,
        biplot_scores=None,
        sample_constraints=None,
        constraint_ids=None,
        proportion_explained=None,
    ):
        self.short_method_name = short_method_name
        self.long_method_name = long_method_name
        self.eigvals = eigvals

        self.samples = samples
        if sample_ids is None:
            no_samp_ids = True
            self.sample_ids = _extract_row_ids(samples)
        else:
            no_samp_ids = False
            self.sample_ids = sample_ids

        self.features = features
        if feature_ids is None and features is not None:
            self.feature_ids = _extract_row_ids(features, warn_ids=no_samp_ids)
        else:
            self.feature_ids = feature_ids

        self.biplot_scores = biplot_scores
        self.sample_constraints = sample_constraints
        self.constraint_ids = constraint_ids
        self.proportion_explained = proportion_explained

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

        """
        lines = ["Ordination results:"]
        method = "%s (%s)" % (self.long_method_name, self.short_method_name)
        lines.append(self._format_attribute(method, "Method", str))

        attrs = [
            (self.eigvals, "Eigvals"),
            (self.proportion_explained, "Proportion explained"),
            (self.features, "Features"),
            (self.samples, "Samples"),
            (self.biplot_scores, "Biplot Scores"),
            (self.sample_constraints, "Sample constraints"),
        ]
        for attr, attr_label in attrs:

            def formatter(e):
                return "x".join(["%d" % s for s in e.shape])

            lines.append(self._format_attribute(attr, attr_label, formatter))

        lines.append(
            self._add_id_line(
                attr_name="Feature IDs", ids=self.feature_ids, data=self.features
            )
        )
        lines.append(
            self._add_id_line(
                attr_name="Sample IDs", ids=self.sample_ids, data=self.samples
            )
        )

        return "\n".join(lines)

    def _add_id_line(self, attr_name, ids, data):
        """Helper to append ids to str."""
        if ids is not None:
            return "\t%s: %s" % (attr_name, _pprint_strs(ids))
        elif data is not None:
            return self._format_attribute(
                data, attr_name, _pprint_strs(_extract_row_ids)
            )
        else:
            return "\t%s: N/A" % attr_name

    def _format_attribute(self, attr, attr_label, formatter):
        if attr is None:
            formatted_attr = "N/A"
        else:
            formatted_attr = formatter(attr)
        return "\t%s: %s" % (attr_label, formatted_attr)

    def plot(
        self,
        df=None,
        column=None,
        axes=None,
        axis_labels=None,
        title="",
        cmap=None,
        s=20,
        centroids=False,
        confidence_ellipses=False,
    ):
        """Create a scatterplot of ordination results colored by metadata.

        Creates a scatterplot of the ordination results, where each point
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
            Indices of sample coordinates to plot. Must contain exactly two or three
            elements, for 2D or 3D plots, respectively. For example, ``axes=(0, 1, 2)``
            (default if there are three or more dimensions) will create a 3D plot with
            PC1 on the x-axis, PC2 on the y-axis, and PC3 on the z-axis.
            ``axes=(0, 1)`` (default if there are only two dimensions) will create a
            2D plot with PC1 on the x-axis and PC2 on the y-axis.
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
        centroids : bool, optional
            If True, plot the centroids of each category in `column`.
        confidence_ellipses : bool, optional
            If True, plot confidence ellipses for each category in `column`.
            Ellipses are calculated using to fit an interval of 2 standard deviations
            using the covariance of the points. Only supported for 2D plots.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing the scatterplot and legend/colorbar if metadata
            were provided.

        Raises
        ------
        ValueError
            Raised on invalid input, including the following situations:

            - there are not at least two dimensions to plot
            - there are not exactly two or three values in `axes`, they are not unique,
              or are out of range
            - there are not exactly two or three values in `axis_labels`
            - either `df` or `column` is provided without the other
            - `column` is not in the ``DataFrame``
            - sample IDs in the ordination results are not in `df` or have missing data
              in `column`
            - confidence ellipses are requested for 3D plots


        Notes
        -----
        This method creates basic plots of ordination results, and is intended
        to provide a quick look at the results in the context of metadata
        (e.g., from within the Jupyter Lab). For more customization and to
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

        >>> fig = pcoa_results.plot(
        ...     df=df, column='body_site',
        ...     title='Samples colored by body site',
        ...     cmap='Set1', s=50
        ... )  # doctest: +SKIP

        """
        # Note: New features should not be added to this method and should
        # instead be added to EMPeror (http://biocore.github.io/emperor/).
        # Only bug fixes and minor updates should be made to this method.

        self._get_mpl_plt()

        # This handles any input, numpy/pandas/polars
        coord_matrix = np.atleast_2d(self.samples).T

        # Determine default axes based on available dimensions
        if axes is None:
            num_dims = coord_matrix.shape[0]
            if num_dims >= 3:
                axes = [0, 1, 2]
            elif num_dims == 2:
                axes = [0, 1]
            else:
                raise ValueError(
                    "At least two dimensions are required to plot "
                    "ordination results. There is only %d dimension." % num_dims
                )

        self._validate_plot_axes(coord_matrix, axes)

        point_colors, category_to_color = self._get_plot_point_colors(
            df, column, self.sample_ids, cmap
        )

        # Validate metadata requirements for centroids and ellipses
        if category_to_color is None and centroids is True:
            raise ValueError("Metadata must be provided to plot centroids.")

        if category_to_color is None and confidence_ellipses is True:
            raise ValueError("Metadata must be provided to plot confidence ellipses.")

        # Check if confidence ellipses are requested for 3D plot
        if len(axes) == 3 and confidence_ellipses is True:
            raise ValueError("Confidence ellipses can only be plotted in 2D.")

        # Create figure and axes
        is_3d = len(axes) == 3
        if is_3d:
            fig = self.plt.figure()
            ax = fig.add_subplot(projection="3d")
            xs, ys, zs = (
                coord_matrix[axes[0]],
                coord_matrix[axes[1]],
                coord_matrix[axes[2]],
            )
        else:
            fig, ax = self.plt.subplots()
            xs, ys = coord_matrix[axes[0]], coord_matrix[axes[1]]
            zs = None

        # Create scatter plot
        if zs is None:
            scatter_fn = functools.partial(ax.scatter, xs, ys, s=s)
        else:
            scatter_fn = functools.partial(ax.scatter, xs, ys, zs, s=s)

        if point_colors is None:
            plot = scatter_fn()
        else:
            plot = scatter_fn(c=point_colors)

        # Add centroids if requested
        if centroids and category_to_color:
            self._plot_centroids(ax, df, column, category_to_color, axes, is_3d)

        # Add confidence ellipses if requested (2D only)
        if confidence_ellipses and not is_3d and category_to_color:
            self._plot_confidence_ellipses(ax, df, column, category_to_color, axes)

        # Set axis labels
        if axis_labels is None:
            axis_labels = ["%d" % axis for axis in axes]
        elif len(axis_labels) != len(axes):
            raise ValueError(
                f"axis_labels ({len(axis_labels)} elements) "
                f"must contain the same number of elements as "
                f"axes ({len(axes)} elements)."
            )

        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        if is_3d:
            ax.set_zlabel(axis_labels[2])
            ax.set_zticklabels([])

        ax.set_title(title)

        # Create legend/colorbar
        if point_colors is not None:
            if category_to_color is None:
                fig.colorbar(plot)
            else:
                self._plot_categorical_legend(
                    ax,
                    category_to_color,
                    centroids=centroids,
                    confidence_ellipses=confidence_ellipses,
                )

        return fig

    def _plot_centroids(self, ax, df, column, category_to_color, axes, is_3d):
        """Plot centroids for each category.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes object to plot on.
        df : pd.DataFrame
            DataFrame containing sample metadata.
        column : str
            Column name in df for grouping.
        category_to_color : dict
            Mapping of category labels to colors.
        axes : list of int
            Indices of axes to plot.
        is_3d : bool
            Whether the plot is 3D.

        """
        current_centroids = self.samples.groupby(df[column]).mean()

        for label, color in category_to_color.items():
            if label not in current_centroids.index:
                continue

            centroid = current_centroids.loc[label]

            if is_3d:
                ax.scatter(
                    centroid.iloc[axes[0]],
                    centroid.iloc[axes[1]],
                    centroid.iloc[axes[2]],
                    color=color,
                    marker="x",
                    s=30,
                    label=f"'{label}' centroid",
                )
            else:
                ax.scatter(
                    centroid.iloc[axes[0]],
                    centroid.iloc[axes[1]],
                    color=color,
                    marker="x",
                    s=30,
                    label=f"'{label}' centroid",
                )

    def _plot_confidence_ellipses(self, ax, df, column, category_to_color, axes):
        """Plot confidence ellipses for each category in 2D.

        Ellipses are calculated to fit an interval of 2 standard deviations
        using the covariance of the points.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes object to plot on.
        df : pd.DataFrame
            DataFrame containing sample metadata.
        column : str
            Column name in df for grouping.
        category_to_color : dict
            Mapping of category labels to colors.
        axes : list of int
            Indices of axes to plot (must be length 2).

        Notes
        -----
        Derived from:
        https://matplotlib.org/stable/gallery/statistics/confidence_ellipse.html

        """
        for label, color in category_to_color.items():
            group = self.samples[df[column] == label]

            if len(group) < 3:
                continue  # can't draw ellipse with less than 3 points

            x_vals = group.iloc[:, axes[0]]
            y_vals = group.iloc[:, axes[1]]

            # Covariance matrix
            cov = np.cov(x_vals, y_vals)

            # Pearson correlation coefficient
            pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])

            # Ellipse radii
            ell_radius_x = np.sqrt(1 + pearson)
            ell_radius_y = np.sqrt(1 - pearson)

            # Means
            mean_x = x_vals.mean()
            mean_y = y_vals.mean()

            # 2 standard deviations for confidence interval
            scale_x = np.sqrt(cov[0, 0]) * 2
            scale_y = np.sqrt(cov[1, 1]) * 2

            # Rotation angle
            angle = 0.5 * np.degrees(np.arctan2(2 * cov[0, 1], cov[0, 0] - cov[1, 1]))

            # Create and transform ellipse
            ellipse = self.mpl.patches.Ellipse(
                (0, 0),
                width=ell_radius_x * 2,
                height=ell_radius_y * 2,
                facecolor="none",
                edgecolor=color,
                lw=2,
                label=f"'{label}' ellipse",
            )

            transf = (
                self.mpl.transforms.Affine2D()
                .rotate_deg(angle)
                .scale(scale_x, scale_y)
                .translate(mean_x, mean_y)
            )
            ellipse.set_transform(transf + ax.transData)
            ax.add_patch(ellipse)

    def _validate_plot_axes(self, coord_matrix, axes):
        """Validate `axes` against coordinates matrix."""
        num_dims = coord_matrix.shape[0]
        if num_dims < 2:
            raise ValueError(
                "At least two dimensions are required to plot "
                "ordination results. There are only %d "
                "dimension(s)." % num_dims
            )
        if len(axes) not in [2, 3]:
            raise ValueError(
                "`axes` must contain exactly two or three elements "
                "(found %d elements)." % len(axes)
            )
        if len(set(axes)) != len(axes):
            raise ValueError("The values provided for `axes` must be unique.")

        for idx, axis in enumerate(axes):
            if axis < 0 or axis >= num_dims:
                raise ValueError("`axes[%d]` must be >= 0 and < %d." % (idx, num_dims))

    def _get_plot_point_colors(self, df, column, ids, cmap):
        """Return a list of colors for each plot point given a metadata column.

        If `column` is categorical, additionally returns a dictionary mapping
        each category (str) to color (used for legend creation).

        """
        if (df is None and column is not None) or (df is not None and column is None):
            raise ValueError(
                "Both df and column must be provided, or both must be None."
            )
        elif df is None and column is None:
            point_colors, category_to_color = None, None
        else:
            if column not in df:
                raise ValueError("Column '%s' not in data frame." % column)

            col_vals = df.reindex(ids, axis=0).loc[:, column]

            if col_vals.isnull().any():
                raise ValueError(
                    "One or more IDs in the ordination results "
                    "are not in the data frame, or there is "
                    "missing data in the data frame's '%s' "
                    "column." % column
                )

            category_to_color = None
            try:
                point_colors = col_vals.astype(float)
            except ValueError:
                # we have categorical data, so choose a color for each
                # category, where colors are evenly spaced across the
                # colormap.
                # derived from http://stackoverflow.com/a/14887119
                categories = col_vals.unique()
                cmap = self.plt.get_cmap(cmap)
                category_colors = cmap(np.linspace(0, 1, len(categories)))

                category_to_color = dict(zip(categories, category_colors))
                point_colors = col_vals.apply(lambda x: category_to_color[x])

            point_colors = point_colors.tolist()

        return point_colors, category_to_color

    def _plot_categorical_legend(
        self, ax, color_dict, centroids=False, confidence_ellipses=False
    ):
        """Add legend to plot using specified mapping of category to color."""
        # derived from http://stackoverflow.com/a/20505720
        proxies = []
        labels = []
        for category in color_dict:
            proxy = self.mpl.lines.Line2D(
                [0], [0], linestyle="none", c=color_dict[category], marker="o"
            )
            proxies.append(proxy)
            labels.append(category)

            # Add centroid markers if enabled
            if centroids:
                proxy = self.mpl.lines.Line2D(
                    [0],
                    [0],
                    linestyle="none",
                    c=color_dict[category],
                    marker="x",
                    markersize=np.sqrt(30),
                )
                proxies.append(proxy)
                labels.append(f"'{category}' centroid")
            # Add confidence ellipse lines if enabled
            if confidence_ellipses:
                proxy = self.mpl.lines.Line2D([0], [0], c=color_dict[category], lw=2)
                proxies.append(proxy)
                labels.append(f"'{category}' ellipse")

        # place legend outside of the axes (centered)
        # derived from http://matplotlib.org/users/legend_guide.html
        ax.legend(
            proxies,
            labels,
            numpoints=1,
            loc=6,
            bbox_to_anchor=(1.05, 0.5),
            borderaxespad=0.0,
        )

    def rename(self, mapper, matrix="samples", strict=True):
        r"""Rename sample or feature IDs in the data matrix.

        Parameters
        ----------
        mapper : dict or callable
            A dictionary or function that maps current IDs to new IDs.
        matrix : str, optional
            Specifies which matrix contains the IDs to be renamed. Either
            "samples" (default) or "features".
        strict : bool, optional
           If ``True`` (default), every ID in the matrix must be included in
           ``mapper``. If ``False``, only the specified IDs will be renamed.

        Raises
        ------
        ValueError
            If ``mapper`` does not contain all of the same IDs in the matrix
            whereas in strict mode.
        ValueError
            If renaming features but self does not contain features.
        ValueError
            If ``matrix`` is neither "samples" nor "features".

        """
        if matrix not in ("samples", "features"):
            raise ValueError('Matrix must be either "samples" or "features".')
        df = getattr(self, matrix)
        if matrix == "features" and df is None:
            raise ValueError(
                "`features` were not provided on the construction of this object."
            )
        if strict and isinstance(mapper, dict) and not set(df.index).issubset(mapper):
            raise ValueError(
                "The IDs in mapper do not include all IDs in the %s matrix." % matrix
            )
        df.rename(index=mapper, inplace=True)
