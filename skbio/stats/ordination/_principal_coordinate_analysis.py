# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from numbers import Integral
from warnings import warn

import numpy as np
import pandas as pd
from numpy import dot, hstack
from numpy.linalg import qr, svd
from scipy.linalg import eigh

from skbio.util import get_rng
from skbio.stats.distance import DistanceMatrix
from ._ordination_results import OrdinationResults
from ._utils import center_distance_matrix, scale


def pcoa(
    distance_matrix,
    method="eigh",
    number_of_dimensions=0,
    inplace=False,
    seed=None,
    warn_neg_eigval=0.01,
):
    r"""Perform Principal Coordinate Analysis (PCoA).

    PCoA is an ordination method similar to Principal Components Analysis (PCA), with
    the difference that it operates on distance matrices, calculated using meaningful
    and typically non-Euclidian methods.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        The input distance matrix.
    method : str, optional
        Matrix decomposition method to use. Default is "eigh" (eigendecomposition),
        which computes exact eigenvectors and eigenvalues for all dimensions. The
        alternate is "fsvd" (fast singular value decomposition), a heuristic that can
        compute only a given number of dimensions.
    number_of_dimensions : int or float, optional
        Dimensions to reduce the distance matrix to. This number determines how many
        eigenvectors and eigenvalues will be returned. If an integer is provided, the
        exact number of dimensions will be retained. If a float between 0 and 1, it
        represents the fractional cumulative variance to be retained. Default is 0,
        which will retain the same number of dimensions as the distance matrix.
    inplace : bool, optional
        If True, the input distance matrix will be centered in-place to reduce memory
        consumption, at the cost of losing the original distances. Default is False.
    seed : int or np.random.Generator, optional
        A user-provided random seed or random generator instance for method "fsvd".
        See :func:`details <skbio.util.get_rng>`.

        .. versionadded:: 0.6.3

    warn_neg_eigval : bool or float, optional
        Raise a warning if any negative eigenvalue is obtained and its magnitude
        exceeds the specified fraction threshold compared to the largest positive
        eigenvalue, which suggests potential inaccuracy in the PCoA result. Default is
        0.01. Set True to warn regardless of the magnitude. Set False to disable
        warning completely.

        .. versionadded:: 0.6.3

    Returns
    -------
    OrdinationResults
        Object that stores the PCoA results, including eigenvalues, the proportion
        explained by each of them, and transformed sample coordinates.

    See Also
    --------
    OrdinationResults

    Notes
    -----
    Principal Coordinate Analysis (PCoA) was first described in [1]_.

    This function uses a choice of two methods for matrix decomposition: The default
    method, ``eigh``, performs eigendecomposition, an exact method that computes all
    eigenvectors and eigenvalues. The alternative method, ``fsvd``, performs fast
    singular value decomposition (FSVD) [2]_, an efficient heuristic method that
    allows a custom number of dimensions to be specified to reduce calculation at the
    cost of losing accuracy. The degree of accuracy lost is dependent on dataset.

    Note that the default method ``eigh`` does not natively support a given number of
    dimensions to reduce a matrix to. Therefore, if this parameter is specified, all
    eigenvectors and eigenvalues will be simply be computed with no speed gain, and
    only the specified number of dimensions will be returned.

    Eigenvalues represent the magnitude of individual principal coordinates, and
    they are usually positive. However, negative eigenvalues can occur when the
    distances were calculated using a non-Euclidean metric that does not satisfy
    triangle inequality. If the negative eigenvalues are small in magnitude compared
    to the largest positive eigenvalue, it is usually safe to ignore them. However,
    large negative eigenvalues may indicate result inaccuracy, in which case a warning
    message will be displayed. The paramter ``warn_neg_eigval`` controls the threshold
    for the warning.

    PCoA on Euclidean distances is equivalent to Principal Component Analysis (PCA).
    However, in ecology, the Euclidean distance preserved by PCA is often not a good
    choice because it deals poorly with double zeros. For example, species have
    unimodal distributions along environmental gradients. If a species is absent from
    two sites simultaneously, it can't be known if an environmental variable is too
    high in one of them and too low in the other, or too low in both, etc. On the other
    hand, if a species is present in two sites, that means that the sites are similar.

    Note that the returned eigenvectors are not normalized to unit length.

    References
    ----------
    .. [1] Gower, J. C. (1966). Some distance properties of latent root and vector
       methods used in multivariate analysis. Biometrika, 53(3-4), 325-338.

    .. [2] Halko, N., Martinsson, P. G., Shkolnisky, Y., & Tygert, M. (2011). An
       algorithm for the principal component analysis of large data sets. SIAM
       Journal on Scientific computing, 33(5), 2580-2594.

    """
    distance_matrix = DistanceMatrix(distance_matrix)

    # Center distance matrix, a requirement for PCoA here
    matrix_data = center_distance_matrix(distance_matrix.data, inplace=inplace)

    # If no dimension specified, by default will compute all eigenvectors
    # and eigenvalues
    if number_of_dimensions == 0:
        if method == "fsvd" and matrix_data.shape[0] > 10:
            warn(
                "FSVD: since no value for number_of_dimensions is specified, "
                "PCoA for all dimensions will be computed, which may "
                "result in long computation time if the original "
                "distance matrix is large.",
                RuntimeWarning,
            )

        # distance_matrix is guaranteed to be square
        number_of_dimensions = matrix_data.shape[0]
    elif number_of_dimensions < 0:
        raise ValueError(
            "Invalid operation: cannot reduce distance matrix "
            "to negative dimensions using PCoA. Did you intend "
            'to specify the default value "0", which sets '
            "the number_of_dimensions equal to the "
            "dimensionality of the given distance matrix?"
        )
    elif not isinstance(number_of_dimensions, Integral) and number_of_dimensions > 1:
        raise ValueError(
            "Invalid operation: A floating-point number greater than 1 cannot be "
            "supplied as the number of dimensions."
        )

    if warn_neg_eigval and not 0 <= warn_neg_eigval <= 1:
        raise ValueError(
            "warn_neg_eigval must be Boolean or a floating-point number between 0 "
            "and 1."
        )

    # Perform eigendecomposition
    if method == "eigh":
        # eigh does not natively support specifying number_of_dimensions, i.e.
        # there are no speed gains unlike in FSVD. Later, we slice off unwanted
        # dimensions to conform the result of eigh to the specified
        # number_of_dimensions.

        eigvals, eigvecs = eigh(matrix_data)
        long_method_name = "Principal Coordinate Analysis"
    elif method == "fsvd":
        # new parameter for num_dimensions = number of dimensions (accounting for
        # non-int values)
        num_dimensions = number_of_dimensions
        if 0 < number_of_dimensions < 1:
            warn(
                "FSVD: since value for number_of_dimensions is specified as float, "
                "PCoA for all dimensions will be computed, which may "
                "result in long computation time if the original "
                "distance matrix is large. "
                "Consider specifying an integer value to optimize performance.",
                RuntimeWarning,
            )
            num_dimensions = matrix_data.shape[0]
        eigvals, eigvecs = _fsvd(matrix_data, num_dimensions, seed=seed)
        long_method_name = "Approximate Principal Coordinate Analysis using FSVD"
    else:
        raise ValueError(
            "PCoA eigendecomposition method {} not supported.".format(method)
        )
    # Ensure number_of_dimensions does not exceed available dimensions
    # number_of_dimensions = min(number_of_dimensions, eigvals.shape[0])

    # cogent makes eigenvalues positive by taking the
    # abs value, but that doesn't seem to be an approach accepted
    # by L&L to deal with negative eigenvalues. We raise a warning
    # in that case. First, we make values close to 0 equal to 0.
    negative_close_to_zero = np.isclose(eigvals, 0)
    eigvals[negative_close_to_zero] = 0

    # eigvals might not be ordered, so we first sort them, then analogously
    # sort the eigenvectors by the ordering of the eigenvalues too
    idxs_descending = eigvals.argsort()[::-1]
    eigvals = eigvals[idxs_descending]
    eigvecs = eigvecs[:, idxs_descending]

    # large negative eigenvalues suggest result inaccuracy
    # see: https://github.com/scikit-bio/scikit-bio/issues/1410
    if warn_neg_eigval and eigvals[-1] < 0:
        if warn_neg_eigval is True or -eigvals[-1] > eigvals[0] * warn_neg_eigval:
            warn(
                "The result contains negative eigenvalues that are large in magnitude,"
                " which may suggest result inaccuracy. See Notes for details. The"
                " negative-most eigenvalue is {0} whereas the largest positive one is"
                " {1}.".format(eigvals[-1], eigvals[0]),
                RuntimeWarning,
            )

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
    if 0 < number_of_dimensions < 1:
        cumulative_variance = np.cumsum(proportion_explained)
        num_dimensions = (
            np.searchsorted(cumulative_variance, number_of_dimensions, side="left") + 1
        )
        # gives the number of dimensions needed to reach specified variance
        # updates number of dimensions to reach the requirement of variance.
        number_of_dimensions = num_dimensions

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
        samples=pd.DataFrame(
            coordinates, index=distance_matrix.ids, columns=axis_labels
        ),
        proportion_explained=pd.Series(proportion_explained, index=axis_labels),
    )


def _fsvd(centered_distance_matrix, number_of_dimensions=10, seed=None):
    """Perform singular value decomposition.

    More specifically in this case eigendecomposition, using fast heuristic algorithm
    nicknamed "FSVD" (FastSVD), adapted and optimized from the algorithm described
    by Halko et al (2011).

    Parameters
    ----------
    centered_distance_matrix : np.array
       Numpy matrix representing the distance matrix for which the
       eigenvectors and eigenvalues shall be computed
    number_of_dimensions : int
       Number of dimensions to keep. Must be lower than or equal to the
       rank of the given distance_matrix.
    seed : int or np.random.Generator, optional
        A user-provided random seed or random generator instance.

    Returns
    -------
    np.array
       Array of eigenvectors, each with number_of_dimensions length.
    np.array
       Array of eigenvalues, a total number of number_of_dimensions.

    Notes
    -----
    The algorithm is based on [1]_.

    Ported from MATLAB implementation described in [2]_.

    References
    ----------
    .. [1] Halko, N., Martinsson, P. G., Shkolnisky, Y., & Tygert, M. (2011). An
       algorithm for the principal component analysis of large data sets. SIAM
       Journal on Scientific computing, 33(5), 2580-2594.

    .. [2] https://stats.stackexchange.com/a/11934/211065

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
        raise ValueError("FSVD expects square distance matrix")

    if number_of_dimensions > m or number_of_dimensions > n:
        raise ValueError(
            "FSVD: number_of_dimensions cannot be larger than"
            " the dimensionality of the given distance matrix."
        )

    if number_of_dimensions < 0:
        raise ValueError(
            "Invalid operation: cannot reduce distance matrix "
            "to negative dimensions using PCoA. Did you intend "
            'to specify the default value "0", which sets '
            "the number_of_dimensions equal to the "
            "dimensionality of the given distance matrix?"
        )

    k = number_of_dimensions + 2

    # Form a real nxl matrix G whose entries are independent, identically
    # distributed Gaussian random variables of zero mean and unit variance
    rng = get_rng(seed)
    G = rng.standard_normal(size=(n, k))

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
        H = hstack((H, dot(centered_distance_matrix, dot(centered_distance_matrix, H))))
        for x in range(3, num_levels + 2):
            tmp = dot(centered_distance_matrix, dot(centered_distance_matrix, H))

            H = hstack(
                (H, dot(centered_distance_matrix, dot(centered_distance_matrix, tmp)))
            )

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


def pcoa_biplot(ordination, y):
    """Compute the projection of descriptors into a PCoA matrix.

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

    Notes
    -----
    This implementation is as described in Chapter 9 of [1]_.

    References
    ----------
    .. [1] Legendre P. and Legendre L. 1998. Numerical Ecology. Elsevier, Amsterdam.

    """
    # acknowledge that most saved ordinations lack a name, however if they have
    # a name, it should be PCoA
    if ordination.short_method_name != "" and ordination.short_method_name != "PCoA":
        raise ValueError(
            "This biplot computation can only be performed in a " "PCoA matrix."
        )

    if set(y.index) != set(ordination.samples.index):
        raise ValueError(
            "The eigenvectors and the descriptors must describe " "the same samples."
        )

    eigvals = ordination.eigvals.values
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
    Uproj = np.sqrt(N - 1) * spc.dot(
        np.diag(np.power(eigvals, -0.5, where=eigvals > 0))
    )

    ordination.features = pd.DataFrame(
        data=Uproj, index=y.columns.copy(), columns=coordinates.columns.copy()
    )
    ordination.features.fillna(0.0, inplace=True)

    return ordination
