# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from numpy import dot, hstack
from numpy.linalg import qr, svd
from numpy.random import standard_normal
from scipy.linalg import eigh
from warnings import warn

from skbio.stats.distance import DistanceMatrix
from skbio.util._decorator import experimental
from ._ordination_results import OrdinationResults
from ._utils import center_distance_matrix, scale


@experimental(as_of="0.4.0")
def pcoa(distance_matrix, method="eigh", number_of_dimensions=0,
         inplace=False):
    r"""Perform Principal Coordinate Analysis.

    Principal Coordinate Analysis (PCoA) is a method similar
    to Principal Components Analysis (PCA) with the difference that PCoA
    operates on distance matrices, typically with non-euclidian and thus
    ecologically meaningful distances like UniFrac in microbiome research.

    In ecology, the euclidean distance preserved by Principal
    Component Analysis (PCA) is often not a good choice because it
    deals poorly with double zeros (Species have unimodal
    distributions along environmental gradients, so if a species is
    absent from two sites at the same site, it can't be known if an
    environmental variable is too high in one of them and too low in
    the other, or too low in both, etc. On the other hand, if an
    species is present in two sites, that means that the sites are
    similar.).

    Note that the returned eigenvectors are not normalized to unit length.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        A distance matrix.
    method : str, optional
        Eigendecomposition method to use in performing PCoA.
        By default, uses SciPy's `eigh`, which computes exact
        eigenvectors and eigenvalues for all dimensions. The alternate
        method, `fsvd`, uses faster heuristic eigendecomposition but loses
        accuracy. The magnitude of accuracy lost is dependent on dataset.
    number_of_dimensions : int, optional
        Dimensions to reduce the distance matrix to. This number determines
        how many eigenvectors and eigenvalues will be returned.
        By default, equal to the number of dimensions of the distance matrix,
        as default eigendecomposition using SciPy's `eigh` method computes
        all eigenvectors and eigenvalues. If using fast heuristic
        eigendecomposition through `fsvd`, a desired number of dimensions
        should be specified. Note that the default eigendecomposition
        method `eigh` does not natively support a specifying number of
        dimensions to reduce a matrix to, so if this parameter is specified,
        all eigenvectors and eigenvalues will be simply be computed with no
        speed gain, and only the number specified by `number_of_dimensions`
        will be returned. Specifying a value of `0`, the default, will
        set `number_of_dimensions` equal to the number of dimensions of the
        specified `distance_matrix`.
    inplace : bool, optional
        If true, centers a distance matrix in-place in a manner that reduces
        memory consumption.

    Returns
    -------
    OrdinationResults
        Object that stores the PCoA results, including eigenvalues, the
        proportion explained by each of them, and transformed sample
        coordinates.

    See Also
    --------
    OrdinationResults

    Notes
    -----
    .. note:: If the distance is not euclidean (for example if it is a
        semimetric and the triangle inequality doesn't hold),
        negative eigenvalues can appear. There are different ways
        to deal with that problem (see Legendre & Legendre 1998, \S
        9.2.3), but none are currently implemented here.
        However, a warning is raised whenever negative eigenvalues
        appear, allowing the user to decide if they can be safely
        ignored.
    """
    distance_matrix = DistanceMatrix(distance_matrix)

    # Center distance matrix, a requirement for PCoA here
    matrix_data = center_distance_matrix(distance_matrix.data, inplace=inplace)

    # If no dimension specified, by default will compute all eigenvectors
    # and eigenvalues
    if number_of_dimensions == 0:
        if method == "fsvd" and matrix_data.shape[0] > 10:
            warn("FSVD: since no value for number_of_dimensions is specified, "
                 "PCoA for all dimensions will be computed, which may "
                 "result in long computation time if the original "
                 "distance matrix is large.", RuntimeWarning)

        # distance_matrix is guaranteed to be square
        number_of_dimensions = matrix_data.shape[0]
    elif number_of_dimensions < 0:
        raise ValueError('Invalid operation: cannot reduce distance matrix '
                         'to negative dimensions using PCoA. Did you intend '
                         'to specify the default value "0", which sets '
                         'the number_of_dimensions equal to the '
                         'dimensionality of the given distance matrix?')

    # Perform eigendecomposition
    if method == "eigh":
        # eigh does not natively support specifying number_of_dimensions, i.e.
        # there are no speed gains unlike in FSVD. Later, we slice off unwanted
        # dimensions to conform the result of eigh to the specified
        # number_of_dimensions.

        eigvals, eigvecs = eigh(matrix_data)
        long_method_name = "Principal Coordinate Analysis"
    elif method == "fsvd":
        eigvals, eigvecs = _fsvd(matrix_data, number_of_dimensions)
        long_method_name = "Approximate Principal Coordinate Analysis " \
                           "using FSVD"
    else:
        raise ValueError(
            "PCoA eigendecomposition method {} not supported.".format(method))

    # cogent makes eigenvalues positive by taking the
    # abs value, but that doesn't seem to be an approach accepted
    # by L&L to deal with negative eigenvalues. We raise a warning
    # in that case. First, we make values close to 0 equal to 0.
    negative_close_to_zero = np.isclose(eigvals, 0)
    eigvals[negative_close_to_zero] = 0
    if np.any(eigvals < 0):
        warn(
            "The result contains negative eigenvalues."
            " Please compare their magnitude with the magnitude of some"
            " of the largest positive eigenvalues. If the negative ones"
            " are smaller, it's probably safe to ignore them, but if they"
            " are large in magnitude, the results won't be useful. See the"
            " Notes section for more details. The smallest eigenvalue is"
            " {0} and the largest is {1}.".format(eigvals.min(),
                                                  eigvals.max()),
            RuntimeWarning
        )

    # eigvals might not be ordered, so we first sort them, then analogously
    # sort the eigenvectors by the ordering of the eigenvalues too
    idxs_descending = eigvals.argsort()[::-1]
    eigvals = eigvals[idxs_descending]
    eigvecs = eigvecs[:, idxs_descending]

    # If we return only the coordinates that make sense (i.e., that have a
    # corresponding positive eigenvalue), then Jackknifed Beta Diversity
    # won't work as it expects all the OrdinationResults to have the same
    # number of coordinates. In order to solve this issue, we return the
    # coordinates that have a negative eigenvalue as 0
    num_positive = (eigvals >= 0).sum()
    eigvecs[:, num_positive:] = np.zeros(eigvecs[:, num_positive:].shape)
    eigvals[num_positive:] = np.zeros(eigvals[num_positive:].shape)

    if method == "fsvd":
        # Since the dimension parameter, hereafter referred to as 'd',
        # restricts the number of eigenvalues and eigenvectors that FSVD
        # computes, we need to use an alternative method to compute the sum
        # of all eigenvalues, used to compute the array of proportions
        # explained. Otherwise, the proportions calculated will only be
        # relative to d number of dimensions computed; whereas we want
        # it to be relative to the entire dimensionality of the
        # centered distance matrix.

        # An alternative method of calculating th sum of eigenvalues is by
        # computing the trace of the centered distance matrix.
        # See proof outlined here: https://goo.gl/VAYiXx
        sum_eigenvalues = np.trace(matrix_data)
    else:
        # Calculate proportions the usual way
        sum_eigenvalues = np.sum(eigvals)

    proportion_explained = eigvals / sum_eigenvalues

    # In case eigh is used, eigh computes all eigenvectors and -values.
    # So if number_of_dimensions was specified, we manually need to ensure
    # only the requested number of dimensions
    # (number of eigenvectors and eigenvalues, respectively) are returned.
    eigvecs = eigvecs[:, :number_of_dimensions]
    eigvals = eigvals[:number_of_dimensions]
    proportion_explained = proportion_explained[:number_of_dimensions]

    # Scale eigenvalues to have length = sqrt(eigenvalue). This
    # works because np.linalg.eigh returns normalized
    # eigenvectors. Each row contains the coordinates of the
    # objects in the space of principal coordinates. Note that at
    # least one eigenvalue is zero because only n-1 axes are
    # needed to represent n points in a euclidean space.
    coordinates = eigvecs * np.sqrt(eigvals)

    axis_labels = ["PC%d" % i for i in range(1, number_of_dimensions + 1)]
    return OrdinationResults(
        short_method_name="PCoA",
        long_method_name=long_method_name,
        eigvals=pd.Series(eigvals, index=axis_labels),
        samples=pd.DataFrame(coordinates, index=distance_matrix.ids,
                             columns=axis_labels),
        proportion_explained=pd.Series(proportion_explained,
                                       index=axis_labels))


def _fsvd(centered_distance_matrix, number_of_dimensions=10):
    """
    Performs singular value decomposition, or more specifically in
    this case eigendecomposition, using fast heuristic algorithm
    nicknamed "FSVD" (FastSVD), adapted and optimized from the algorithm
    described by Halko et al (2011).

    Parameters
    ----------
    centered_distance_matrix : np.array
       Numpy matrix representing the distance matrix for which the
       eigenvectors and eigenvalues shall be computed
    number_of_dimensions : int
       Number of dimensions to keep. Must be lower than or equal to the
       rank of the given distance_matrix.

    Returns
    -------
    np.array
       Array of eigenvectors, each with number_of_dimensions length.
    np.array
       Array of eigenvalues, a total number of number_of_dimensions.

    Notes
    -----
    The algorithm is based on 'An Algorithm for the Principal
    Component analysis of Large Data Sets'
    by N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
    Original Paper: https://arxiv.org/abs/1007.5510

    Ported from MATLAB implementation described here:
    https://stats.stackexchange.com/a/11934/211065
    """

    m, n = centered_distance_matrix.shape

    # Number of levels of the Krylov method to use.
    # For most applications, num_levels=1 or num_levels=2 is sufficient.
    num_levels = 1

    # Changes the power of the spectral norm, thus minimizing the error).
    use_power_method = False

    # Note: a (conjugate) transpose is removed for performance, since we
    # only expect square matrices.
    if m != n:
        raise ValueError('FSVD expects square distance matrix')

    if number_of_dimensions > m or number_of_dimensions > n:
        raise ValueError('FSVD: number_of_dimensions cannot be larger than'
                         ' the dimensionality of the given distance matrix.')

    if number_of_dimensions < 0:
        raise ValueError('Invalid operation: cannot reduce distance matrix '
                         'to negative dimensions using PCoA. Did you intend '
                         'to specify the default value "0", which sets '
                         'the number_of_dimensions equal to the '
                         'dimensionality of the given distance matrix?')

    k = number_of_dimensions + 2

    # Form a real nxl matrix G whose entries are independent, identically
    # distributed Gaussian random variables of zero mean and unit variance
    G = standard_normal(size=(n, k))

    if use_power_method:
        # use only the given exponent
        H = dot(centered_distance_matrix, G)

        for x in range(2, num_levels + 2):
            # enhance decay of singular values
            # note: distance_matrix is no longer transposed, saves work
            # since we're expecting symmetric, square matrices anyway
            # (Daniel McDonald's changes)
            H = dot(centered_distance_matrix, dot(centered_distance_matrix, H))

    else:
        # compute the m x l matrices H^{(0)}, ..., H^{(i)}
        # Note that this is done implicitly in each iteration below.
        H = dot(centered_distance_matrix, G)
        # to enhance performance
        H = hstack(
            (H,
             dot(centered_distance_matrix, dot(centered_distance_matrix, H))))
        for x in range(3, num_levels + 2):
            tmp = dot(centered_distance_matrix,
                      dot(centered_distance_matrix, H))

            H = hstack(
                (H, dot(centered_distance_matrix,
                        dot(centered_distance_matrix, tmp))))

    # Using the pivoted QR-decomposition, form a real m * ((i+1)l) matrix Q
    # whose columns are orthonormal, s.t. there exists a real
    # ((i+1)l) * ((i+1)l) matrix R for which H = QR
    Q, R = qr(H)

    # Compute the n * ((i+1)l) product matrix T = A^T Q
    T = dot(centered_distance_matrix, Q)  # step 3

    # Form an SVD of T
    Vt, St, W = svd(T, full_matrices=False)
    W = W.transpose()

    # Compute the m * ((i+1)l) product matrix
    Ut = dot(Q, W)

    U_fsvd = Ut[:, :number_of_dimensions]

    S = St[:number_of_dimensions]

    # drop imaginary component, if we got one
    # Note:
    #   In cogent, after computing eigenvalues/vectors, the imaginary part
    #   is dropped, if any. We know for a fact that the eigenvalues are
    #   real, so that's not necessary, but eigenvectors can in principle
    #   be complex (see for example
    #   http://math.stackexchange.com/a/47807/109129 for details)
    eigenvalues = S.real
    eigenvectors = U_fsvd.real

    return eigenvalues, eigenvectors


@experimental(as_of="0.5.3")
def pcoa_biplot(ordination, y):
    """Compute the projection of descriptors into a PCoA matrix

    This implementation is as described in Chapter 9 of Legendre & Legendre,
    Numerical Ecology 3rd edition.

    Parameters
    ----------
    ordination: OrdinationResults
        The computed principal coordinates analysis of dimensions (n, c) where
        the matrix ``y`` will be projected onto.
    y: DataFrame
        Samples by features table of dimensions (n, m). These can be
        environmental features or abundance counts. This table should be
        normalized in cases of dimensionally heterogenous physical variables.

    Returns
    -------
    OrdinationResults
        The modified input object that includes projected features onto the
        ordination space in the ``features`` attribute.
    """

    # acknowledge that most saved ordinations lack a name, however if they have
    # a name, it should be PCoA
    if (ordination.short_method_name != '' and
            ordination.short_method_name != 'PCoA'):
        raise ValueError('This biplot computation can only be performed in a '
                         'PCoA matrix.')

    if set(y.index) != set(ordination.samples.index):
        raise ValueError('The eigenvectors and the descriptors must describe '
                         'the same samples.')

    eigvals = ordination.eigvals
    coordinates = ordination.samples
    N = coordinates.shape[0]

    # align the descriptors and eigenvectors in a sample-wise fashion
    y = y.reindex(coordinates.index)

    # S_pc from equation 9.44
    # Represents the covariance matrix between the features matrix and the
    # column-centered eigenvectors of the pcoa.
    spc = (1 / (N - 1)) * y.values.T.dot(scale(coordinates, ddof=1))

    # U_proj from equation 9.55, is the matrix of descriptors to be projected.
    #
    # Only get the power of non-zero values, otherwise this will raise a
    # divide by zero warning. There shouldn't be negative eigenvalues(?)
    Uproj = np.sqrt(N - 1) * spc.dot(np.diag(np.power(eigvals, -0.5,
                                                      where=eigvals > 0)))

    ordination.features = pd.DataFrame(data=Uproj,
                                       index=y.columns.copy(),
                                       columns=coordinates.columns.copy())
    ordination.features.fillna(0.0, inplace=True)

    return ordination
