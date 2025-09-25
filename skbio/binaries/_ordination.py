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


def pcoa_fsvd_available(
    distance_matrix,
    dimensions,
    inplace=False,
    seed=None,
):
    """Check if the scikit-bio-binaries shared library provides the pcoa functionality.

    Parameters
    ----------
    distance_matrix : np.ndarray or DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
    dimensions : int
        Dimensions to reduce the distance matrix to. This number determines how many
        eigenvectors and eigenvalues will be returned.
    inplace : bool
        If True, the input distance matrix will be centered in-place to reduce memory
        consumption, at the cost of losing the original distances. Default is False.
    seed : int or np.random.Generator, optional
        A user-provided random seed or random generator instance.

    Returns
    -------
    boolean
        If False, the pcoa_fsvd function will raise an exception.

    Note
    ----
    Internally it uses caching, to minimize overhead.

    """
    # v1 is the minimum version that supports pcoa_fsvd
    # includes check for library existence
    if _skbb_get_api_version() >= 1:  # pragma: no cover
        # Do some basic sanity checks
        # check it is a positive number
        if not isinstance(dimensions, Integral):
            return False
        elif dimensions < 1:
            return False
        return True
    else:  # pragma: no cover
        return False


def pcoa_fsvd(
    distance_matrix,
    dimensions,
    inplace=False,
    seed=None,
):
    r"""Perform Principal Coordinate Analysis (PCoA).

    PCoA is an ordination method similar to Principal Components Analysis (PCA), with
    the difference that it operates on distance matrices, calculated using meaningful
    and typically non-Euclidian methods.

    Parameters
    ----------
    distance_matrix : np.ndarray or DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
    dimensions : int
        Dimensions to reduce the distance matrix to. This number determines how many
        eigenvectors and eigenvalues will be returned.
    inplace : bool
        If True, the input distance matrix will be centered in-place to reduce memory
        consumption, at the cost of losing the original distances. Default is False.
    seed : int or np.random.Generator, optional
        A user-provided random seed or random generator instance.

    Returns
    -------
    tuple of 3 np.ndarray
        eigenvalues, transformed samples and the proportion explained by each of them.


    Notes
    -----
    Principal Coordinate Analysis (PCoA) was first described in [1]_.

    It used the ``fsvd`` method, which performs fast
    singular value decomposition (FSVD) [2]_, an efficient heuristic method that
    allows a custom number of dimensions to be specified to reduce calculation at the
    cost of losing accuracy. The degree of accuracy lost is dependent on dataset.

    Eigenvalues represent the magnitude of individual principal coordinates, and
    they are usually positive. However, negative eigenvalues can occur when the
    distances were calculated using a non-Euclidean metric that does not satisfy
    triangle inequality. If the negative eigenvalues are small in magnitude compared
    to the largest positive eigenvalue, it is usually safe to ignore them. However,
    large negative eigenvalues may indicate result inaccuracy, in which case a warning
    message will be displayed.

    PCoA on Euclidean distances is equivalent to Principal Component Analysis (PCA).
    However, in ecology, the Euclidean distance preserved by PCA is often not a good
    choice because it deals poorly with double zeros. For example, species have
    unimodal distributions along environmental gradients. If a species is absent from
    two sites simultaneously, it can't be known if an environmental variable is too
    high in one of them and too low in the other, or too low in both, etc. On the other
    hand, if a species is present in two sites, that means that the sites are similar.

    Note that the returned eigenvectors are not normalized to unit length.

    Finally, the result is undefined if distance_matrix is a np.ndarray
    that does not represent a valid distance matrix.

    References
    ----------
    .. [1] Gower, J. C. (1966). Some distance properties of latent root and vector
       methods used in multivariate analysis. Biometrika, 53(3-4), 325-338.

    .. [2] Halko, N., Martinsson, P. G., Shkolnisky, Y., & Tygert, M. (2011). An
       algorithm for the principal component analysis of large data sets. SIAM
       Journal on Scientific computing, 33(5), 2580-2594.

    """
    # minimum version that support pcoa
    if _skbb_get_api_version() >= 1:  # pragma: no cover
        if not isinstance(dimensions, Integral):
            raise ValueError("dimensions must be an integer value")
        if dimensions < 1:
            raise ValueError("dimensions must be a positive number")
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
        if dimensions > distance_matrix_shape0:
            raise ValueError("dimensions cannot be larger than matrix size")
        # create output buffers
        eigenvalues = np.ndarray(shape=(dimensions,), dtype=distance_matrix_data.dtype)
        proportion_explained = np.ndarray(
            shape=(dimensions,), dtype=distance_matrix_data.dtype
        )
        samples = np.ndarray(
            shape=(distance_matrix_shape0, dimensions),
            dtype=distance_matrix_data.dtype,
            order="C",
        )
        # ready to call the C functions
        dll = _get_skbb_dll()
        i_mdim = ctypes.c_uint(distance_matrix_shape0)
        i_mat = distance_matrix_data.ctypes.data_as(ctypes.c_void_p)
        i_n_eigh = ctypes.c_uint(dimensions)
        i_seed = ctypes.c_int(int_seed)
        o_ev = eigenvalues.ctypes.data_as(ctypes.c_void_p)
        o_sp = samples.ctypes.data_as(ctypes.c_void_p)
        o_pe = proportion_explained.ctypes.data_as(ctypes.c_void_p)
        if distance_matrix_data.dtype == np.dtype("float64"):
            if inplace:
                dll.skbb_pcoa_fsvd_inplace_fp64(
                    i_mdim, i_mat, i_n_eigh, i_seed, o_ev, o_sp, o_pe
                )
            else:
                dll.skbb_pcoa_fsvd_fp64(
                    i_mdim, i_mat, i_n_eigh, i_seed, o_ev, o_sp, o_pe
                )
        elif distance_matrix_data.dtype == np.dtype("float32"):
            if inplace:
                dll.skbb_pcoa_fsvd_inplace_fp32(
                    i_mdim, i_mat, i_n_eigh, i_seed, o_ev, o_sp, o_pe
                )
            else:
                dll.skbb_pcoa_fsvd_fp32(
                    i_mdim, i_mat, i_n_eigh, i_seed, o_ev, o_sp, o_pe
                )
        else:
            raise TypeError("distance_matrix type must be either float32 or float64")
        # if we got here, everything went well
        return eigenvalues, samples, proportion_explained
    else:  # pragma: no cover
        raise ImportError("skbb_pcoa_fsvd not available")
