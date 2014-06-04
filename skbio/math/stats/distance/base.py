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
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

from skbio.core.distance import DistanceMatrix


class CategoricalStats(object):
    """Base class for categorical statistical methods.

    Categorical statistical methods generally test for significant differences
    between discrete groups of objects, as determined by a categorical variable
    (grouping vector).

    See Also
    --------
    ANOSIM, PERMANOVA

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
        p_value_str = self._format_p_value(self.p_value, self.permutations)

        return (self.short_method_name, '%d' % self.sample_size,
                '%d' % len(self.groups), str(self.statistic), p_value_str,
                '%d' % self.permutations)

    def _format_p_value(self, p_value, permutations):
        """Format p-value as a string with the correct number of decimals.

        Number of decimals is determined by the number of permutations.
        """
        if p_value is None:
            result = 'N/A'
        elif permutations < 10:
            # This can be the last step of a long process, so we don't want to
            # fail.
            result = ('Too few permutations to compute p-value (permutations '
                      '= %d)' % permutations)
        else:
            decimal_places = int(np.log10(permutations + 1))
            result = ('%1.' + '%df' % decimal_places) % p_value

        return result


# This function is defined here instead of in its own bioenv.py module since
# running 'nosetests --with-doctest' fails if a function is defined in a module
# with the same name. When the following issue is fixed in nose, consider
# moving this to its own module:
#
#     https://github.com/nose-devs/nose/issues/92
#
# Either way, this shouldn't affect users as the public API/import stays the
# same regardless of where we put the function.
def bioenv(distance_matrix, data_frame, columns=None):
    """Find subset of variables maximally correlated with distances.

    Finds subsets of variables whose Euclidean distances (after scaling the
    variables; see Notes section below for details) are maximally
    rank-correlated with the distance matrix. For example, the distance matrix
    might contain distances between communities, and the variables might be
    numeric environmental variables (e.g., pH). Correlation between the
    community distance matrix and Euclidean environmental distance matrix is
    computed using Spearman's rank correlation coefficient (:math:`\\rho`).

    Subsets of environmental variables range in size from 1 to the total number
    of variables (inclusive). For example, if there are 3 variables, the "best"
    variable subsets will be computed for subset sizes 1, 2, and 3.

    The "best" subset is chosen by computing the correlation between the
    community distance matrix and all possible Euclidean environmental distance
    matrices at the given subset size. The combination of environmental
    variables with maximum correlation is chosen as the "best" subset.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
    data_frame : pandas.DataFrame
        Contains columns of variables (e.g., numeric environmental variables
        such as pH) associated with the objects in `distance_matrix`. Must be
        indexed by the IDs in `distance_matrix` (i.e., the row labels must be
        distance matrix IDs), but the order of IDs between `distance_matrix`
        and `data_frame` need not be the same. All IDs in the distance matrix
        must be present in `data_frame`. Extra IDs in `data_frame` are allowed
        (they are ignored in the calculations).
    columns : iterable of strs, optional
        Column names in `data_frame` to include as variables in the
        calculations. If not provided, defaults to all columns in `data_frame`.
        The values in each column must be numeric or convertible to a numeric
        type.

    Returns
    -------
    pandas.DataFrame
        Data frame containing the "best" subset of variables at each subset
        size, as well as the correlation coefficient of each.

    Raises
    ------
    TypeError
        If invalid input types are provided, or if one or more specified
        columns in `data_frame` are not numeric.
    ValueError
        If column name(s) or `distance_matrix` IDs cannot be found in
        `data_frame`, if there is missing data (``NaN``) in the environmental
        variables, or if the environmental variables cannot be scaled (e.g.,
        due to zero variance).

    See Also
    --------
    scipy.stats.spearmanr

    Notes
    -----
    See [1]_ for the original method reference (originally called BIO-ENV).
    The general algorithm and interface are similar to ``vegan::bioenv``,
    available in R's vegan package [2]_. This method can also be found in
    PRIMER-E [3]_ (originally called BIO-ENV, but is now called BEST).

    .. warning:: This method can take a *long* time to run if a large number of
       variables are specified, as all possible subsets are evaluated at each
       subset size.

    The variables are scaled before computing the Euclidean distance: each
    column is centered and then scaled by its standard deviation.

    References
    ----------
    .. [1] Clarke, K. R & Ainsworth, M. 1993. "A method of linking multivariate
       community structure to environmental variables". Marine Ecology Progress
       Series, 92, 205-219.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    .. [3] http://www.primer-e.com/primer.htm

    Examples
    --------
    Import the functionality we'll use in the following examples. The call to
    ``pd.set_option`` ensures consistent data frame formatting across
    different versions of pandas. This call is not necessary for normal
    use; it is only included here so that the doctests will pass.

    >>> import pandas as pd
    >>> from skbio.core.distance import DistanceMatrix
    >>> from skbio.math.stats.distance import bioenv
    >>> pd.set_option('show_dimensions', True) # not necessary for normal use

    Load a 4x4 community distance matrix:

    >>> dm = DistanceMatrix([[0.0, 0.5, 0.25, 0.75],
    ...                      [0.5, 0.0, 0.1, 0.42],
    ...                      [0.25, 0.1, 0.0, 0.33],
    ...                      [0.75, 0.42, 0.33, 0.0]],
    ...                     ['A', 'B', 'C', 'D'])

    Load a ``pandas.DataFrame`` with two environmental variables, pH and
    elevation:

    >>> df = pd.DataFrame([[7.0, 400],
    ...                    [8.0, 530],
    ...                    [7.5, 450],
    ...                    [8.5, 810]],
    ...                   index=['A','B','C','D'],
    ...                   columns=['pH', 'Elevation'])

    Note that the data frame is indexed with the same IDs (``'A'``, ``'B'``,
    ``'C'``, and ``'D'``) that are in the distance matrix. This is necessary in
    order to link the environmental variables (metadata) to each of the objects
    in the distance matrix. In this example, the IDs appear in the same order
    in both the distance matrix and data frame, but this is not necessary.

    Find the best subsets of environmental variables that are correlated with
    community distances:

    >>> bioenv(dm, df) # doctest: +NORMALIZE_WHITESPACE
                   size  correlation
    vars
    pH                1     0.771517
    pH, Elevation     2     0.714286
    <BLANKLINE>
    [2 rows x 2 columns]

    We see that in this simple example, pH alone is maximally rank-correlated
    with the community distances (:math:`\\rho=0.771517`).

    """
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("Must provide a DistanceMatrix as input.")
    if not isinstance(data_frame, pd.DataFrame):
        raise TypeError("Must provide a pandas.DataFrame as input.")

    if columns is None:
        columns = data_frame.columns.values.tolist()

    if len(set(columns)) != len(columns):
        raise ValueError("Duplicate column names are not supported.")

    if len(columns) < 1:
        raise ValueError("Must provide at least one column.")

    for column in columns:
        if column not in data_frame:
            raise ValueError("Column '%s' not in data frame." % column)

    # Subset and order the vars data frame to match the IDs in the distance
    # matrix, only keeping the specified columns.
    vars_df = data_frame.loc[distance_matrix.ids, columns]

    if vars_df.isnull().any().any():
        raise ValueError("One or more IDs in the distance matrix are not "
                         "in the data frame, or there is missing data in the "
                         "data frame.")

    try:
        vars_df = vars_df.astype(float)
    except ValueError:
        raise TypeError("All specified columns in the data frame must be "
                        "numeric.")

    # Scale the vars and extract the underlying numpy array from the data
    # frame. We mainly do this for performance as we'll be taking subsets of
    # columns within a tight loop and using a numpy array ends up being ~2x
    # faster.
    vars_array = _scale(vars_df).values
    dm_flat = distance_matrix.condensed_form()

    num_vars = len(columns)
    var_idxs = np.arange(num_vars)

    # For each subset size, store the best combination of variables:
    #     (string identifying best vars, subset size, rho)
    max_rhos = np.empty(num_vars, dtype=[('vars', object),
                                         ('size', int),
                                         ('correlation', float)])
    for subset_size in range(1, num_vars + 1):
        max_rho = None
        for subset_idxs in combinations(var_idxs, subset_size):
            # Compute Euclidean distances using the current subset of
            # variables. pdist returns the distances in condensed form.
            vars_dm_flat = pdist(vars_array[:, subset_idxs],
                                 metric='euclidean')
            rho = spearmanr(dm_flat, vars_dm_flat)[0]

            # If there are ties for the best rho at a given subset size, choose
            # the first one in order to match vegan::bioenv's behavior.
            if max_rho is None or rho > max_rho[0]:
                max_rho = (rho, subset_idxs)

        vars_label = ', '.join([columns[i] for i in max_rho[1]])
        max_rhos[subset_size - 1] = (vars_label, subset_size, max_rho[0])

    return pd.DataFrame.from_records(max_rhos, index='vars')


def _scale(df):
    """Center and scale each column in a data frame.

    Each column is centered (by subtracting the mean) and then scaled by its
    standard deviation.

    """
    # Modified from http://stackoverflow.com/a/18005745
    df = df.copy()
    df -= df.mean()
    df /= df.std()

    if df.isnull().any().any():
        raise ValueError("Column(s) in the data frame could not be scaled, "
                         "likely because the column(s) had no variance.")
    return df
