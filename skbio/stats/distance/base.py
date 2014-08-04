# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils.six import StringIO

import csv

import numpy as np
import pandas as pd

from skbio.core.distance import DistanceMatrix
from skbio.stats import p_value_to_str


class CategoricalStats(object):
    """Base class for categorical statistical methods.

    Categorical statistical methods generally test for significant differences
    between discrete groups of objects, as determined by a categorical variable
    (grouping vector).

    See Also
    --------
    ANOSIM
    PERMANOVA

    """

    short_method_name = ''
    long_method_name = ''
    test_statistic_name = ''

    def __init__(self, distance_matrix, grouping, column=None):
        if not isinstance(distance_matrix, DistanceMatrix):
            raise TypeError("Input must be a DistanceMatrix.")

        if isinstance(grouping, pd.DataFrame):
            if column is None:
                raise ValueError("Must provide a column name if supplying a "
                                 "data frame.")
            else:
                grouping = self._df_to_vector(distance_matrix, grouping,
                                              column)
        elif column is not None:
            raise ValueError("Must provide a data frame if supplying a column "
                             "name.")

        if len(grouping) != distance_matrix.shape[0]:
            raise ValueError("Grouping vector size must match the number of "
                             "IDs in the distance matrix.")

        # Find the group labels and convert grouping to an integer vector
        # (factor).
        groups, grouping = np.unique(grouping, return_inverse=True)

        if len(groups) == len(grouping):
            raise ValueError("All values in the grouping vector are unique. "
                             "This method cannot operate on a grouping vector "
                             "with only unique values (e.g., there are no "
                             "'within' distances because each group of "
                             "objects contains only a single object).")
        if len(groups) == 1:
            raise ValueError("All values in the grouping vector are the same. "
                             "This method cannot operate on a grouping vector "
                             "with only a single group of objects (e.g., "
                             "there are no 'between' distances because there "
                             "is only a single group).")

        self._dm = distance_matrix
        self._grouping = grouping
        self._groups = groups
        self._tri_idxs = np.triu_indices(self._dm.shape[0], k=1)

    def _df_to_vector(self, distance_matrix, df, column):
        """Return a grouping vector from a data frame column.

        Parameters
        ----------
        distance_marix : DistanceMatrix
            Distance matrix whose IDs will be mapped to group labels.
        df : pandas.DataFrame
            Data frame (indexed by distance matrix ID).
        column : str
            Column name in `df` containing group labels.

        Returns
        -------
        list
            Grouping vector (vector of labels) based on the IDs in
            `distance_matrix`. Each ID's label is looked up in the data frame
            under the column specified by `column`.

        Raises
        ------
        ValueError
            If `column` is not in the data frame, or a distance matrix ID is
            not in the data frame.

        """
        if column not in df:
            raise ValueError("Column '%s' not in data frame." % column)

        grouping = df.loc[distance_matrix.ids, column]
        if grouping.isnull().any():
            raise ValueError("One or more IDs in the distance matrix are not "
                             "in the data frame.")
        return grouping.tolist()

    def __call__(self, permutations=999):
        """Execute the statistical method.

        Parameters
        ----------
        permutations : int, optional
            Number of permutations to use when calculating statistical
            significance. Must be >= 0. If 0, the resulting p-value will be
            ``None``.

        Returns
        -------
        CategoricalStatsResults
            Results of the method, including test statistic and p-value.

        .. shownumpydoc

        """
        if permutations < 0:
            raise ValueError("Number of permutations must be greater than or "
                             "equal to zero.")

        stat = self._run(self._grouping)

        p_value = None
        if permutations > 0:
            perm_stats = np.empty(permutations, dtype=np.float64)

            for i in range(permutations):
                perm_grouping = np.random.permutation(self._grouping)
                perm_stats[i] = self._run(perm_grouping)

            p_value = ((perm_stats >= stat).sum() + 1) / (permutations + 1)

        return CategoricalStatsResults(self.short_method_name,
                                       self.long_method_name,
                                       self.test_statistic_name,
                                       self._dm.shape[0], self._groups, stat,
                                       p_value, permutations)

    def _run(self, grouping):
        raise NotImplementedError("Subclasses must implement _run().")


class CategoricalStatsResults(object):
    """Statistical method results container.

    Stores the results of running a `CategoricalStats` method a single time,
    and provides a way to format the results.

    Attributes
    ----------
    short_method_name
    long_method_name
    test_statistic_name
    sample_size
    groups
    statistic
    p_value
    permutations

    Notes
    -----
    Users will generally not directly instantiate objects of this class. The
    various categorical statistical methods will return an object of this type
    when they are run.

    """

    def __init__(self, short_method_name, long_method_name,
                 test_statistic_name, sample_size, groups, statistic, p_value,
                 permutations):
        self.short_method_name = short_method_name
        self.long_method_name = long_method_name
        self.test_statistic_name = test_statistic_name
        self.sample_size = sample_size
        self.groups = groups
        self.statistic = statistic
        self.p_value = p_value
        self.permutations = permutations

    def __str__(self):
        """Return pretty-print (fixed width) string."""
        rows = (self._format_header(), self._format_data())

        max_widths = []
        for col_idx in range(len(rows[0])):
            max_widths.append(max(map(lambda e: len(e[col_idx]), rows)))

        results = []
        for row in rows:
            padded_row = []
            for col_idx, val in enumerate(row):
                padded_row.append(val.rjust(max_widths[col_idx]))
            results.append('  '.join(padded_row))

        return '\n'.join(results) + '\n'

    def _repr_html_(self):
        """Return a string containing an HTML table of results.

        This method will be called within the IPython Notebook instead of
        __repr__ to display results.

        """
        header = self._format_header()
        data = self._format_data()
        return pd.DataFrame([data[1:]], columns=header[1:],
                            index=[data[0]])._repr_html_()

    def summary(self, delimiter='\t'):
        """Return a formatted summary of results as a string.

        The string is formatted as delimited text.

        Parameters
        ----------
        delimiter : str, optional
            String to delimit fields by in formatted output. Default is tab
            (TSV).

        Returns
        -------
        str
            Delimited-text summary of results.

        """
        summary = StringIO()
        csv_writer = csv.writer(summary, delimiter=delimiter,
                                lineterminator='\n')
        csv_writer.writerow(self._format_header())
        csv_writer.writerow(self._format_data())
        return summary.getvalue()

    def _format_header(self):
        return ('Method name', 'Sample size', 'Number of groups',
                self.test_statistic_name, 'p-value', 'Number of permutations')

    def _format_data(self):
        p_value_str = p_value_to_str(self.p_value, self.permutations)

        return (self.short_method_name, '%d' % self.sample_size,
                '%d' % len(self.groups), str(self.statistic), p_value_str,
                '%d' % self.permutations)
