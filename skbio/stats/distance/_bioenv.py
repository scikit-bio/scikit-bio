# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

from skbio.stats.distance import DistanceMatrix


def bioenv(distance_matrix, data_frame, columns=None):
    r"""Find subset of variables maximally correlated with distances.

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
    Import the functionality we'll use in the following examples:

    >>> import pandas as pd
    >>> from skbio import DistanceMatrix
    >>> from skbio.stats.distance import bioenv

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

    We see that in this simple example, pH alone is maximally rank-correlated
    with the community distances (:math:`\rho=0.771517`).

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
    vars_df = data_frame.reindex(distance_matrix.ids, axis=0).loc[:, columns]

    if vars_df.isnull().any().any():
        raise ValueError(
            "One or more IDs in the distance matrix are not "
            "in the data frame, or there is missing data in the "
            "data frame."
        )

    try:
        vars_df = vars_df.astype(float)
    except ValueError:
        raise TypeError("All specified columns in the data frame must be " "numeric.")

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
    max_rhos = np.empty(
        num_vars, dtype=[("vars", object), ("size", int), ("correlation", float)]
    )
    for subset_size in range(1, num_vars + 1):
        max_rho = None
        for subset_idxs in combinations(var_idxs, subset_size):
            # Compute Euclidean distances using the current subset of
            # variables. pdist returns the distances in condensed form.
            vars_dm_flat = pdist(vars_array[:, subset_idxs], metric="euclidean")
            rho = spearmanr(dm_flat, vars_dm_flat)[0]

            # If there are ties for the best rho at a given subset size, choose
            # the first one in order to match vegan::bioenv's behavior.
            if max_rho is None or rho > max_rho[0]:
                max_rho = (rho, subset_idxs)

        vars_label = ", ".join([columns[i] for i in max_rho[1]])
        max_rhos[subset_size - 1] = (vars_label, subset_size, max_rho[0])

    return pd.DataFrame.from_records(max_rhos, index="vars")


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
        raise ValueError(
            "Column(s) in the data frame could not be scaled, "
            "likely because the column(s) had no variance."
        )
    return df
