# ----------------------------------------------------------------------------
# Copyright (c) 2025--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import ctypes
import numpy as np
from numbers import Integral
from ..binaries._util import (
    get_api_version as _skbb_get_api_version,
    get_dll as _get_skbb_dll,
    py_to_bin_random_seed,
)


def permanova_available(
    distance_matrix,
    grouping,
    permutations,
    seed=None,
):
    """Is binaries permanova available?

    Check if the scikit-bio-binaries shared library provides
    the permanova functionality.

    Parameters
    ----------
    distance_matrix : np.ndarray or DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
    grouping : 1-D np.ndarray
        Vector indicating the assignment of objects to groups.
        These integers denote which group an object belongs to.
        It must be the same length and in the same order
        as the objects in `distance_matrix`.
    permutations : int
        Number of permutations to use when assessing statistical
        significance. Must be greater than zero.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance.

    Returns
    -------
    boolean
        If False, the permanova function will raise an exception.

    Note
    ----
    Internally it uses caching, to minimize overhead.

    """
    # v1 is the minimum version that support permanova
    # includes check for library existence
    if _skbb_get_api_version() >= 1:  # pragma: no cover
        # Do some basic sanity checks
        # check it is a positive number
        if not isinstance(permutations, Integral):
            return False
        elif permutations < 1:
            return False
        if not isinstance(grouping, np.ndarray):
            return False
        return True
    else:  # pragma: no cover
        return False


def permanova(
    distance_matrix,
    grouping,
    permutations,
    seed=None,
):
    r"""Test for significant differences between groups using PERMANOVA.

    Permutational Multivariate Analysis of Variance (PERMANOVA) is a
    non-parametric method that tests whether two or more groups of objects
    (e.g., samples) are significantly different based on a categorical factor.
    It is conceptually similar to ANOVA except that it operates on a distance
    matrix, which allows for multivariate analysis. PERMANOVA computes a
    pseudo-F statistic.

    Statistical significance is assessed via a permutation test. The assignment
    of objects to groups (`grouping`) is randomly permuted a number of times
    (controlled via `permutations`). A pseudo-F statistic is computed for each
    permutation and the p-value is the proportion of permuted pseudo-F
    statistics that are equal to or greater than the original (unpermuted)
    pseudo-F statistic.

    Parameters
    ----------
    distance_matrix : np.ndarray or DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
    grouping : 1-D np.ndarray
        Vector indicating the assignment of objects to groups.
        These integers denote which group an object belongs to.
        It must be the same length and in the same order
        as the objects in `distance_matrix`.
    permutations : int
        Number of permutations to use when assessing statistical
        significance. Must be greater than zero.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance.

    Returns
    -------
    float pair
        ``test statistic`` and ``p-value``.

    Notes
    -----
    See [1]_ for the original method reference, as well as ``vegan::adonis``,
    available in R's vegan package [2]_.

    Finally, the result is undefined if distance_matrix is a np.ndarray
    that does not represent a valid distance matrix.

    References
    ----------
    .. [1] Anderson, Marti J. "A new method for non-parametric multivariate
       analysis of variance." Austral Ecology 26.1 (2001): 32-46.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    """
    # minimum version that support pcoa
    if _skbb_get_api_version() >= 1:  # pragma: no cover
        if permutations < 1:
            raise ValueError("permutations must be a positive number")
        int_seed = py_to_bin_random_seed(seed)
        if isinstance(distance_matrix, np.ndarray):
            # already a raw matrix, just use
            distance_matrix_data = distance_matrix
        else:
            # we are assuming it is a DistanceMatrix object
            # get internal representations
            distance_matrix_data = distance_matrix.data
        distance_matrix_shape0 = distance_matrix_data.shape[0]
        if (
            len(distance_matrix_data.shape) != 2
            or distance_matrix_data.shape[1] != distance_matrix_shape0
        ):
            raise TypeError("distance_matrix not square")
        if not isinstance(grouping, np.ndarray):
            raise TypeError("grouping not ndarray")
        if len(grouping.shape) != 1:
            raise TypeError("grouping not 1D ndarray")
        distance_matrix_shape0 = distance_matrix_data.shape[0]
        if grouping.shape[0] != distance_matrix_shape0:
            raise TypeError("grouping not the same size as distance_matrix")
        # skbb expects uint32 array, so convert if needed (cheap)
        if grouping.dtype == np.dtype("uint32"):
            grouping_data = grouping
        else:
            grouping_data = grouping.astype(np.uint32)
        # ready to call the C functions
        dll = _get_skbb_dll()
        i_mdim = ctypes.c_uint(distance_matrix_shape0)
        i_mat = distance_matrix_data.ctypes.data_as(ctypes.c_void_p)
        i_grp = grouping_data.ctypes.data_as(ctypes.c_void_p)
        i_n_perm = ctypes.c_uint(permutations)
        i_seed = ctypes.c_int(int_seed)
        if distance_matrix_data.dtype == np.dtype("float64"):
            o_fstat = ctypes.c_double()
            o_pvalue = ctypes.c_double()
            dll.skbb_permanova_fp64(
                i_mdim,
                i_mat,
                i_grp,
                i_n_perm,
                i_seed,
                ctypes.byref(o_fstat),
                ctypes.byref(o_pvalue),
            )
        elif distance_matrix_data.dtype == np.dtype("float32"):
            o_fstat = ctypes.c_float()
            o_pvalue = ctypes.c_float()
            dll.skbb_permanova_fp32(
                i_mdim,
                i_mat,
                i_grp,
                i_n_perm,
                i_seed,
                ctypes.byref(o_fstat),
                ctypes.byref(o_pvalue),
            )
        else:
            raise TypeError("distance_matrix type must be either float32 or float64")
        # if we got here, everything went well
        return o_fstat.value, o_pvalue.value
    else:  # pragma: no cover
        raise ImportError("skbb_permanova not available")
