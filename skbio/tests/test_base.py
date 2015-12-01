# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
import six
from six import binary_type, text_type

import unittest
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.testing as npt
import pandas as pd
from IPython.core.display import Image, SVG
from nose.tools import assert_is_instance, assert_true

from skbio import OrdinationResults
from skbio._base import (SkbioObject, MetadataMixin, PositionalMetadataMixin,
                         ElasticLines)
from skbio.util._decorator import overrides
from skbio.util._testing import (ReallyEqualMixin, MetadataMixinTests,
                                 PositionalMetadataMixinTests)


class TestSkbioObject(unittest.TestCase):
    def test_no_instantiation(self):
        class Foo(SkbioObject):
            pass

        with self.assertRaises(TypeError):
            Foo()


class TestMetadataMixin(unittest.TestCase, ReallyEqualMixin,
                        MetadataMixinTests):
    def setUp(self):
        class ExampleMetadataMixin(MetadataMixin):
            def __init__(self, metadata=None):
                MetadataMixin._init_(self, metadata=metadata)

            def __eq__(self, other):
                return MetadataMixin._eq_(self, other)

            def __ne__(self, other):
                return MetadataMixin._ne_(self, other)

            def __copy__(self):
                copy = self.__class__(metadata=None)
                copy._metadata = MetadataMixin._copy_(self)
                return copy

            def __deepcopy__(self, memo):
                copy = self.__class__(metadata=None)
                copy._metadata = MetadataMixin._deepcopy_(self, memo)
                return copy

        self._metadata_constructor_ = ExampleMetadataMixin


class TestPositionalMetadataMixin(unittest.TestCase, ReallyEqualMixin,
                                  PositionalMetadataMixinTests):
    def setUp(self):
        class ExamplePositionalMetadataMixin(PositionalMetadataMixin):
            @overrides(PositionalMetadataMixin)
            def _positional_metadata_axis_len_(self):
                return self._axis_len

            def __init__(self, axis_len, positional_metadata=None):
                self._axis_len = axis_len

                PositionalMetadataMixin._init_(
                    self, positional_metadata=positional_metadata)

            def __eq__(self, other):
                return PositionalMetadataMixin._eq_(self, other)

            def __ne__(self, other):
                return PositionalMetadataMixin._ne_(self, other)

            def __copy__(self):
                copy = self.__class__(self._axis_len, positional_metadata=None)
                copy._positional_metadata = \
                    PositionalMetadataMixin._copy_(self)
                return copy

            def __deepcopy__(self, memo):
                copy = self.__class__(self._axis_len, positional_metadata=None)
                copy._positional_metadata = \
                    PositionalMetadataMixin._deepcopy_(self, memo)
                return copy

        self._positional_metadata_constructor_ = ExamplePositionalMetadataMixin


class TestOrdinationResults(unittest.TestCase):
    def setUp(self):
        # Define in-memory CA results to serialize and deserialize.
        eigvals = pd.Series([0.0961330159181, 0.0409418140138], ['CA1', 'CA2'])
        features = np.array([[0.408869425742, 0.0695518116298],
                             [-0.1153860437, -0.299767683538],
                             [-0.309967102571, 0.187391917117]])
        samples = np.array([[-0.848956053187, 0.882764759014],
                            [-0.220458650578, -1.34482000302],
                            [1.66697179591, 0.470324389808]])
        features_ids = ['Species1', 'Species2', 'Species3']
        sample_ids = ['Site1', 'Site2', 'Site3']

        samples_df = pd.DataFrame(samples, index=sample_ids,
                                  columns=['CA1', 'CA2'])
        features_df = pd.DataFrame(features, index=features_ids,
                                   columns=['CA1', 'CA2'])

        self.ordination_results = OrdinationResults(
            'CA', 'Correspondance Analysis', eigvals=eigvals,
            samples=samples_df, features=features_df)

        # DataFrame for testing plot method. Has a categorical column with a
        # mix of numbers and strings. Has a numeric column with a mix of ints,
        # floats, and strings that can be converted to floats. Has a numeric
        # column with missing data (np.nan).
        self.df = pd.DataFrame([['foo', '42', 10],
                                [22, 0, 8],
                                [22, -4.2, np.nan],
                                ['foo', '42.19', 11]],
                               index=['A', 'B', 'C', 'D'],
                               columns=['categorical', 'numeric', 'nancolumn'])

        # Minimal ordination results for easier testing of plotting method.
        # Paired with df above.
        eigvals = np.array([0.50, 0.25, 0.25])
        samples = np.array([[0.1, 0.2, 0.3],
                            [0.2, 0.3, 0.4],
                            [0.3, 0.4, 0.5],
                            [0.4, 0.5, 0.6]])
        samples_df = pd.DataFrame(samples, ['A', 'B', 'C', 'D'],
                                  ['PC1', 'PC2', 'PC3'])

        self.min_ord_results = OrdinationResults(
            'PCoA', 'Principal Coordinate Analysis', eigvals, samples_df)

    def test_str(self):
        exp = ("Ordination results:\n"
               "\tMethod: Correspondance Analysis (CA)\n"
               "\tEigvals: 2\n"
               "\tProportion explained: N/A\n"
               "\tFeatures: 3x2\n"
               "\tSamples: 3x2\n"
               "\tBiplot Scores: N/A\n"
               "\tSample constraints: N/A\n"
               "\tFeature IDs: 'Species1', 'Species2', 'Species3'\n"
               "\tSample IDs: 'Site1', 'Site2', 'Site3'")
        obs = str(self.ordination_results)
        self.assertEqual(obs, exp)

        # all optional attributes missing
        exp = ("Ordination results:\n"
               "\tMethod: Principal Coordinate Analysis (PCoA)\n"
               "\tEigvals: 1\n"
               "\tProportion explained: N/A\n"
               "\tFeatures: N/A\n"
               "\tSamples: 2x1\n"
               "\tBiplot Scores: N/A\n"
               "\tSample constraints: N/A\n"
               "\tFeature IDs: N/A\n"
               "\tSample IDs: 0, 1")
        samples_df = pd.DataFrame(np.array([[1], [2]]))
        obs = str(OrdinationResults('PCoA', 'Principal Coordinate Analysis',
                                    pd.Series(np.array([4.2])), samples_df))
        self.assertEqual(obs.split('\n'), exp.split('\n'))

    def check_basic_figure_sanity(self, fig, exp_num_subplots, exp_title,
                                  exp_legend_exists, exp_xlabel, exp_ylabel,
                                  exp_zlabel):
        # check type
        assert_is_instance(fig, mpl.figure.Figure)

        # check number of subplots
        axes = fig.get_axes()
        npt.assert_equal(len(axes), exp_num_subplots)

        # check title
        ax = axes[0]
        npt.assert_equal(ax.get_title(), exp_title)

        # shouldn't have tick labels
        for tick_label in (ax.get_xticklabels() + ax.get_yticklabels() +
                           ax.get_zticklabels()):
            npt.assert_equal(tick_label.get_text(), '')

        # check if legend is present
        legend = ax.get_legend()
        if exp_legend_exists:
            assert_true(legend is not None)
        else:
            assert_true(legend is None)

        # check axis labels
        npt.assert_equal(ax.get_xlabel(), exp_xlabel)
        npt.assert_equal(ax.get_ylabel(), exp_ylabel)
        npt.assert_equal(ax.get_zlabel(), exp_zlabel)

    def test_plot_no_metadata(self):
        fig = self.min_ord_results.plot()
        self.check_basic_figure_sanity(fig, 1, '', False, '0', '1', '2')

    def test_plot_with_numeric_metadata_and_plot_options(self):
        fig = self.min_ord_results.plot(
            self.df, 'numeric', axes=(1, 0, 2),
            axis_labels=['PC 2', 'PC 1', 'PC 3'], title='a title', cmap='Reds')
        self.check_basic_figure_sanity(
            fig, 2, 'a title', False, 'PC 2', 'PC 1', 'PC 3')

    def test_plot_with_categorical_metadata_and_plot_options(self):
        fig = self.min_ord_results.plot(
            self.df, 'categorical', axes=[2, 0, 1], title='a title',
            cmap='Accent')
        self.check_basic_figure_sanity(fig, 1, 'a title', True, '2', '0', '1')

    def test_plot_with_invalid_axis_labels(self):
        with six.assertRaisesRegex(self, ValueError, 'axis_labels.*4'):
            self.min_ord_results.plot(axes=[2, 0, 1],
                                      axis_labels=('a', 'b', 'c', 'd'))

    def test_validate_plot_axes_valid_input(self):
        # shouldn't raise an error on valid input. nothing is returned, so
        # nothing to check here
        samples = self.min_ord_results.samples.values.T
        self.min_ord_results._validate_plot_axes(samples, (1, 2, 0))

    def test_validate_plot_axes_invalid_input(self):
        # not enough dimensions
        with six.assertRaisesRegex(self, ValueError, '2 dimension\(s\)'):
            self.min_ord_results._validate_plot_axes(
                np.asarray([[0.1, 0.2, 0.3], [0.2, 0.3, 0.4]]), (0, 1, 2))

        coord_matrix = self.min_ord_results.samples.values.T

        # wrong number of axes
        with six.assertRaisesRegex(self, ValueError, 'exactly three.*found 0'):
            self.min_ord_results._validate_plot_axes(coord_matrix, [])
        with six.assertRaisesRegex(self, ValueError, 'exactly three.*found 4'):
            self.min_ord_results._validate_plot_axes(coord_matrix,
                                                     (0, 1, 2, 3))

        # duplicate axes
        with six.assertRaisesRegex(self, ValueError, 'must be unique'):
            self.min_ord_results._validate_plot_axes(coord_matrix, (0, 1, 0))

        # out of range axes
        with six.assertRaisesRegex(self, ValueError, 'axes\[1\].*3'):
            self.min_ord_results._validate_plot_axes(coord_matrix, (0, -1, 2))
        with six.assertRaisesRegex(self, ValueError, 'axes\[2\].*3'):
            self.min_ord_results._validate_plot_axes(coord_matrix, (0, 2, 3))

    def test_get_plot_point_colors_invalid_input(self):
        # column provided without df
        with npt.assert_raises(ValueError):
            self.min_ord_results._get_plot_point_colors(None, 'numeric',
                                                        ['B', 'C'], 'jet')

        # df provided without column
        with npt.assert_raises(ValueError):
            self.min_ord_results._get_plot_point_colors(self.df, None,
                                                        ['B', 'C'], 'jet')

        # column not in df
        with six.assertRaisesRegex(self, ValueError, 'missingcol'):
            self.min_ord_results._get_plot_point_colors(self.df, 'missingcol',
                                                        ['B', 'C'], 'jet')

        # id not in df
        with six.assertRaisesRegex(self, ValueError, 'numeric'):
            self.min_ord_results._get_plot_point_colors(
                self.df, 'numeric', ['B', 'C', 'missingid', 'A'], 'jet')

        # missing data in df
        with six.assertRaisesRegex(self, ValueError, 'nancolumn'):
            self.min_ord_results._get_plot_point_colors(self.df, 'nancolumn',
                                                        ['B', 'C', 'A'], 'jet')

    def test_get_plot_point_colors_no_df_or_column(self):
        obs = self.min_ord_results._get_plot_point_colors(None, None,
                                                          ['B', 'C'], 'jet')
        npt.assert_equal(obs, (None, None))

    def test_get_plot_point_colors_numeric_column(self):
        # subset of the ids in df
        exp = [0.0, -4.2, 42.0]
        obs = self.min_ord_results._get_plot_point_colors(
            self.df, 'numeric', ['B', 'C', 'A'], 'jet')
        npt.assert_almost_equal(obs[0], exp)
        assert_true(obs[1] is None)

        # all ids in df
        exp = [0.0, 42.0, 42.19, -4.2]
        obs = self.min_ord_results._get_plot_point_colors(
            self.df, 'numeric', ['B', 'A', 'D', 'C'], 'jet')
        npt.assert_almost_equal(obs[0], exp)
        assert_true(obs[1] is None)

    def test_get_plot_point_colors_categorical_column(self):
        # subset of the ids in df
        exp_colors = [[0., 0., 0.5, 1.], [0., 0., 0.5, 1.], [0.5, 0., 0., 1.]]
        exp_color_dict = {
            'foo': [0.5, 0., 0., 1.],
            22: [0., 0., 0.5, 1.]
        }
        obs = self.min_ord_results._get_plot_point_colors(
            self.df, 'categorical', ['B', 'C', 'A'], 'jet')
        npt.assert_almost_equal(obs[0], exp_colors)
        npt.assert_equal(obs[1], exp_color_dict)

        # all ids in df
        exp_colors = [[0., 0., 0.5, 1.], [0.5, 0., 0., 1.], [0.5, 0., 0., 1.],
                      [0., 0., 0.5, 1.]]
        obs = self.min_ord_results._get_plot_point_colors(
            self.df, 'categorical', ['B', 'A', 'D', 'C'], 'jet')
        npt.assert_almost_equal(obs[0], exp_colors)
        # should get same color dict as before
        npt.assert_equal(obs[1], exp_color_dict)

    def test_plot_categorical_legend(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # we shouldn't have a legend yet
        assert_true(ax.get_legend() is None)

        self.min_ord_results._plot_categorical_legend(
            ax, {'foo': 'red', 'bar': 'green'})

        # make sure we have a legend now
        legend = ax.get_legend()
        assert_true(legend is not None)

        # do some light sanity checking to make sure our input labels and
        # colors are present. we're not using nose.tools.assert_items_equal
        # because it isn't available in Python 3.
        labels = [t.get_text() for t in legend.get_texts()]
        npt.assert_equal(sorted(labels), ['bar', 'foo'])

        colors = [l.get_color() for l in legend.get_lines()]
        npt.assert_equal(sorted(colors), ['green', 'red'])

    def test_repr_png(self):
        obs = self.min_ord_results._repr_png_()
        assert_is_instance(obs, binary_type)
        assert_true(len(obs) > 0)

    def test_repr_svg(self):
        obs = self.min_ord_results._repr_svg_()
        # print_figure(format='svg') can return text or bytes depending on the
        # version of IPython
        assert_true(isinstance(obs, text_type) or isinstance(obs, binary_type))
        assert_true(len(obs) > 0)

    def test_png(self):
        assert_is_instance(self.min_ord_results.png, Image)

    def test_svg(self):
        assert_is_instance(self.min_ord_results.svg, SVG)


class TestElasticLines(unittest.TestCase):
    def setUp(self):
        self.el = ElasticLines()

    def test_empty(self):
        self.assertEqual(self.el.to_str(), '')

    def test_add_line(self):
        self.el.add_line('foo')
        self.assertEqual(self.el.to_str(), 'foo')

    def test_add_lines(self):
        self.el = ElasticLines()
        self.el.add_lines(['alice', 'bob', 'carol'])
        self.assertEqual(self.el.to_str(), 'alice\nbob\ncarol')

    def test_add_separator(self):
        self.el.add_separator()
        self.assertEqual(self.el.to_str(), '')

        self.el.add_line('foo')
        self.assertEqual(self.el.to_str(), '---\nfoo')

        self.el.add_separator()
        self.el.add_lines(['bar', 'bazzzz'])
        self.el.add_separator()

        self.assertEqual(self.el.to_str(),
                         '------\nfoo\n------\nbar\nbazzzz\n------')


if __name__ == '__main__':
    unittest.main()
