# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from sys import modules
from warnings import warn
from typing import Any, TYPE_CHECKING

import numpy as np

from skbio.util._decorator import aliased, register_aliases, params_aliased
from skbio.util._array import ingest_array

if TYPE_CHECKING:  # pragma: no cover
    from types import ModuleType
    from skbio.util._typing import ArrayLike, StdArray


def _check_composition(
    xp: ModuleType,
    mat: StdArray,
    axis: int = -1,
    nozero: bool = False,
    maxdim: int | None = None,
    allnum: bool = True,
):
    r"""Check if the input matrix contain valid compositions.

    Parameters
    ----------
    xp : namespace
        The array API compatible namespace corresponding ``mat``.
    mat : array of shape (..., n_components, ...)
        A matrix of proportions.
    axis : int, optional
        Axis that represents each composition. Default is the last axis (-1).
    nozero : bool, optional
        If True, matrix cannot have zero values.
    maxdim : int, optional
        Maximum number of dimensions allowed. Default is None.
    allnum : bool, optionsl
        If True, matrix cannot have NaN values.

    Raises
    ------
    TypeError
        If the matrix is not numeric.
    ValueError
        If the matrix contains nan or infinite values.
    ValueError
        If any values in the matrix are negative.
    ValueError
        If there are compositions that have all zeros.
    ValueError
        If the matrix has more than maximum number of dimensions.

    """
    if not xp.isdtype(mat.dtype, "numeric"):
        raise TypeError("Input matrix must have a numeric data type.")
    if allnum:
        # Don't allow infinite or NaN values.
        if not xp.all(xp.isfinite(mat)):
            raise ValueError("Input matrix cannot have infinite or NaN values.")
    else:
        # Allow NaN, but not infinite.
        if xp.any(xp.isinf(mat)):
            raise ValueError("Input matrix cannot have infinite values.")
    if nozero:
        if xp.any(mat <= 0):
            raise ValueError("Input matrix cannot have negative or zero components.")
    else:
        if xp.any(mat < 0):
            raise ValueError("Input matrix cannot have negative components.")
        if xp.any(~xp.any(mat, axis=axis)):
            raise ValueError("Input matrix cannot have compositions with all zeros.")
    if maxdim is not None and mat.ndim > maxdim:
        raise ValueError(f"Input matrix can only have {maxdim} dimensions or less.")


def closure(mat: ArrayLike, axis: int = -1, validate: bool = True) -> StdArray:
    r"""Perform closure to ensure that all components of each composition sum to 1.

    .. versionchanged:: 0.7.0
        The function now works on any dimension in arrays of any number of dimensions.

    Parameters
    ----------
    mat : array_like of shape (..., n_components, ...)
        A matrix of compositions.
    axis : int, optional
        Axis along which closure will be performed. That is, each vector along this
        axis is considered as a composition. Default is the last axis (-1).

        .. versionadded:: 0.7.0

    validate : bool, default True
        Check if the compositions are legitimate.

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray of shape (..., n_components, ...)
        The matrix where all components of each composition sum to 1.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import closure
    >>> X = np.array([[2, 2, 6], [4, 4, 2]])
    >>> closure(X)
    array([[ 0.2,  0.2,  0.6],
           [ 0.4,  0.4,  0.2]])

    """
    xp, mat = ingest_array(mat)
    if validate:
        _check_composition(xp, mat, axis)
    return _closure(xp, mat, axis)


def _closure(xp: ModuleType, mat: StdArray, axis: int = -1) -> StdArray:
    """Perform closure."""
    row_sums = _nansum(xp, mat, axis=axis, keepdims=True)
    return mat / row_sums


def _nansum(
    xp: ModuleType, arr: StdArray, axis: int, keepdims: bool = False
) -> StdArray:
    """Sum array elements along axis, ignoring NaN values."""
    # Create mask of non-NaN values
    nan_mask = xp.isnan(arr)
    # Replace NaN with 0 for summation
    arr_no_nan = xp.where(nan_mask, xp.asarray(0.0, dtype=arr.dtype), arr)
    return xp.sum(arr_no_nan, axis=axis, keepdims=keepdims)


@aliased("multiplicative_replacement", "0.6.0", True)
def multi_replace(mat, delta=None):
    r"""Replace all zeros with small non-zero values.

    It uses the multiplicative replacement strategy [1]_, replacing zeros with
    a small positive :math:`\delta` and ensuring that the compositions still
    add up to 1.

    Parameters
    ----------
    mat : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    delta : float, optional
        A small number to be used to replace zeros. If not specified, the
        default value is :math:`\delta = \frac{1}{N^2}` where :math:`N` is the
        number of components.

    Returns
    -------
    ndarray of shape (n_compositions, n_components)
        The matrix where all of the values are non-zero and each composition
        (row) adds up to 1.

    Raises
    ------
    ValueError
        If negative proportions are created due to a large ``delta``.

    Notes
    -----
    This method will result in negative proportions if a large delta is chosen.

    References
    ----------
    .. [1] J. A. Martin-Fernandez. "Dealing With Zeros and Missing Values in
           Compositional Data Sets Using Nonparametric Imputation"

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import multi_replace
    >>> X = np.array([[.2, .4, .4, 0], [0, .5, .5, 0]])
    >>> multi_replace(X)
    array([[ 0.1875,  0.375 ,  0.375 ,  0.0625],
           [ 0.0625,  0.4375,  0.4375,  0.0625]])

    """
    mat = closure(mat)
    z_mat = mat == 0

    num_feats = mat.shape[-1]
    tot = z_mat.sum(axis=-1, keepdims=True)

    if delta is None:
        delta = (1.0 / num_feats) ** 2

    zcnts = 1 - tot * delta
    if (zcnts < 0).any():
        raise ValueError(
            "Multiplicative replacement created negative proportions. Consider "
            "using a smaller `delta`."
        )
    mat = np.where(z_mat, delta, zcnts * mat)
    return mat.squeeze()


def _closure_two(x, y, validate):
    xp, x, y = ingest_array(x, y)
    if validate:
        _check_composition(xp, x)
        _check_composition(xp, y)
    return xp, _closure(xp, x), _closure(xp, y)


def perturb(x: ArrayLike, y: ArrayLike, validate: bool = True) -> StdArray:
    r"""Perform the perturbation operation.

    This operation is defined as:

    .. math::
        x \oplus y = C[x_1 y_1, \ldots, x_D y_D]

    :math:`C[x]` is the closure operation defined as:

    .. math::
        C[x] = \left[\frac{x_1}{\sum_{i=1}^{D} x_i},\ldots,
                     \frac{x_D}{\sum_{i=1}^{D} x_i} \right]

    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.

    Parameters
    ----------
    x : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    y : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    validate : bool, default True
        Check if the compositions are legitimate.

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray of shape (n_compositions, n_components)
       A matrix of proportions where all of the values are non-zero and each
       composition (row) adds up to 1.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import perturb

    Consider a very simple environment with only three species. The species in
    the environment are evenly distributed and their proportions are equal:

    >>> before = np.array([1/3, 1/3, 1/3])

    Suppose that an antibiotic kills off half of the population for the first
    two species, but doesn't harm the third species. Then the perturbation
    vector would be as follows:

    >>> after = np.array([1/2, 1/2, 1])

    And the resulting perturbation would be:

    >>> perturb(before, after)
    array([ 0.25,  0.25,  0.5 ])

    """
    xp, cx, cy = _closure_two(x, y, validate)
    return _closure(xp, cx * cy)


def perturb_inv(x: ArrayLike, y: ArrayLike, validate: bool = True) -> StdArray:
    r"""Perform the inverse perturbation operation.

    This operation is defined as:

    .. math::
        x \ominus y = C[x_1 y_1^{-1}, \ldots, x_D y_D^{-1}]

    :math:`C[x]` is the closure operation defined as:

    .. math::
        C[x] = \left[\frac{x_1}{\sum_{i=1}^{D} x_i},\ldots,
                     \frac{x_D}{\sum_{i=1}^{D} x_i} \right]

    for some :math:`D` dimensional real vector :math:`x` and :math:`D` is the
    number of components for every composition.

    Parameters
    ----------
    x : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    y : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    validate : bool, default True
        Check if the compositions are legitimate.

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray of shape (n_compositions, n_components)
        A matrix of proportions where all of the values are non-zero and each
        composition (row) adds up to 1.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import perturb_inv
    >>> x = np.array([.1, .3, .4, .2])
    >>> y = np.array([1/6, 1/6, 1/3, 1/3])
    >>> perturb_inv(x, y)
    array([ 0.14285714,  0.42857143,  0.28571429,  0.14285714])

    """
    xp, cx, cy = _closure_two(x, y, validate)
    return _closure(xp, cx / cy)


def power(x: ArrayLike, a: float, validate: bool = True) -> StdArray:
    r"""Perform the power operation.

    This operation is defined as follows:

    .. math::
        `x \odot a = C[x_1^a, \ldots, x_D^a]

    :math:`C[x]` is the closure operation defined as:

    .. math::
        C[x] = \left[\frac{x_1}{\sum_{i=1}^{D} x_i},\ldots,
                     \frac{x_D}{\sum_{i=1}^{D} x_i} \right]

    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.

    Parameters
    ----------
    x : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    a : float
        A scalar exponent.
    validate : bool, default True
        Check if the compositions are legitimate.

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray of shape (n_compositions, n_components)
       The matrix where all of the values are non-zero and each composition
       (row) adds up to 1.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import power
    >>> x = np.array([.1, .3, .4, .2])
    >>> power(x, .1)
    array([ 0.23059566,  0.25737316,  0.26488486,  0.24714631])

    """
    xp, x = ingest_array(x)
    if validate:
        _check_composition(xp, x)
    cx = _closure(xp, x)
    return _closure(xp, cx**a).squeeze()


def inner(x: ArrayLike, y: ArrayLike, validate: bool = True) -> StdArray:
    r"""Calculate the Aitchson inner product.

    This inner product is defined as follows:

    .. math::
        \langle x, y \rangle_a =
        \frac{1}{2D} \sum\limits_{i=1}^{D} \sum\limits_{j=1}^{D}
        \ln\left(\frac{x_i}{x_j}\right) \ln\left(\frac{y_i}{y_j}\right)

    Parameters
    ----------
    x : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    y : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    validate : bool, default True
        Check if the compositions are legitimate.

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray or scalar of shape (n_compositions, n_compositions)
        Inner product result.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import inner
    >>> x = np.array([.1, .3, .4, .2])
    >>> y = np.array([.2, .4, .2, .2])
    >>> inner(x, y)  # doctest: +ELLIPSIS
    0.2107852473...

    """
    xp, cx, cy = _closure_two(x, y, validate)
    clrx, clry = _clr(xp, cx, axis=-1), _clr(xp, cy, axis=-1)
    return xp.matmul(clrx, clry.T)


def clr(mat: ArrayLike, axis: int = -1, validate: bool = True) -> StdArray:
    r"""Perform centre log ratio (CLR) transformation.

    This function transforms compositions from Aitchison geometry to the real
    space. The :math:`clr` transform is both an isometry and an isomorphism
    defined on the following spaces:

    .. math::
        clr: S^D \rightarrow U

    where :math:`U=
    \{x :\sum\limits_{i=1}^D x = 0 \; \forall x \in \mathbb{R}^D\}`

    It is defined for a composition :math:`x` as follows:

    .. math::
        clr(x) = \ln\left[\frac{x_1}{g_m(x)}, \ldots, \frac{x_D}{g_m(x)}\right]

    where :math:`g_m(x) = (\prod\limits_{i=1}^{D} x_i)^{1/D}` is the geometric
    mean of :math:`x`.

    .. versionchanged:: 0.7.0
        The function now works on any dimension in arrays of any number of dimensions.

    Parameters
    ----------
    mat : array_like of shape (..., n_components, ...)
        A matrix of positive proportions.
    axis : int, optional
        Axis along which CLR transformation will be performed. Each vector on this axis
        is considered as a composition. Default is the last axis (-1).

        .. versionadded:: 0.7.0

    validate : bool, default True
        Check if the matrix consists of strictly positive values.

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray of shape (..., n_components, ...)
        CLR-transformed matrix.

    See Also
    --------
    clr_inv

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import clr
    >>> x = np.array([.1, .3, .4, .2])
    >>> clr(x)
    array([-0.79451346,  0.30409883,  0.5917809 , -0.10136628])

    """
    xp, mat = ingest_array(mat)
    if validate:
        _check_composition(xp, mat, nozero=True)
    return _clr(xp, mat, axis)


def _clr(xp: ModuleType, mat: StdArray, axis: int) -> StdArray:
    """Perform CLR transform."""
    return (lmat := xp.log(mat)) - xp.mean(lmat, axis=axis, keepdims=True)


def clr_inv(mat: ArrayLike, axis: int = -1, validate: bool = True) -> StdArray:
    r"""Perform inverse centre log ratio (CLR) transformation.

    This function transforms compositions from the real space to Aitchison
    geometry. The :math:`clr^{-1}` transform is both an isometry, and an
    isomorphism defined on the following spaces:

    .. math::
        clr^{-1}: U \rightarrow S^D

    where :math:`U=
    \{x :\sum\limits_{i=1}^D x = 0 \; \forall x \in \mathbb{R}^D\}`

    This transformation is defined as follows:

    .. math::
        clr^{-1}(x) = C[\exp( x_1, \ldots, x_D)]

    .. versionchanged:: 0.7.0
        The function now works on any dimension in arrays of any number of dimensions.

    Parameters
    ----------
    mat : array_like of shape (..., n_components, ...)
        A matrix of CLR-transformed data.
    axis : int, optional
        Axis along which inverse CLR transformation will be performed. Each vector on
        this axis is considered as a CLR-transformed composition. Default is the last
        axis (-1).

        .. versionadded:: 0.7.0

    validate: bool, optional
        Check if the matrix has been centered at 0. Violation will result in a warning
        rather than an error, for backward compatibility. Defaults to True.

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray of shape (..., n_components, ...)
        Inverse CLR-transformed matrix.

    See Also
    --------
    clr

    Notes
    -----
    The output of ``clr_inv`` is guaranteed to have each composition sum to 1. But this
    property isn't required for the input for ``clr``. Therefore, ``clr_inv`` does not
    completely invert ``clr``. Instead, ``clr_inv(clr(mat))`` and ``closure(mat)`` are
    equal.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import clr_inv
    >>> x = np.array([.1, .3, .4, .2])
    >>> clr_inv(x)
    array([ 0.21383822,  0.26118259,  0.28865141,  0.23632778])

    """
    xp, mat = ingest_array(mat)

    # `1e-8` is taken from `np.allclose`. It's not guaranteed that `xp` has `allclose`,
    # therefore it is manually written here.
    if validate and xp.any(xp.abs(xp.sum(mat, axis=axis)) > 1e-08):
        warn(
            "The input matrix is not in the CLR range, which requires the sum of "
            "values per composition equals to 0.",
            UserWarning,
        )

    return _clr_inv(xp, mat, axis)


def _clr_inv(xp: ModuleType, mat: StdArray, axis: int) -> StdArray:
    """Perform inverse CLR transform."""
    # for numerical stability, shift the values < 1
    diff = xp.exp(mat - xp.max(mat, axis=axis, keepdims=True))
    return _closure(xp, diff, axis)


def rclr(mat: ArrayLike, axis: int = -1, validate: bool = True) -> StdArray:
    r"""Perform robust centre log ratio (rclr) transformation.

    The robust CLR transformation is similar to the standard CLR transformation,
    but it only operates on observed (non-zero) values [1]_. This makes it suitable
    for sparse compositional data.

    For each composition, the transformation computes:

    .. math::

        rclr(x_i) = \ln(x_i) - \frac{1}{|S|} \sum_{j \in S} \ln(x_j)

    where :math:`S` is the set of indices with non-zero values, and :math:`|S|`
    is the number of non-zero values.

    Parameters
    ----------
    mat : array_like of shape (..., n_components, ...)
        A matrix of non-negative values. Zeros are allowed and will become
        NaN in the output. NaN values in the input are preserved (representing
        missing entries).
    axis : int, optional
        Axis along which rclr transformation will be performed. Each vector
        on this axis is considered as a composition. Default is the last
        axis (-1).
    validate : bool, default True
        Check if the matrix consists of non-negative, finite values.
        NaN values are allowed as missing entries.

    Returns
    -------
    ndarray of shape (..., n_components, ...)
        rclr-transformed matrix. Zero values in the input become NaN.

    See Also
    --------
    clr

    Notes
    -----
    The rclr transformation has several advantages for sparse compositional
    data:

    1. It does not require pseudocount addition, which can bias results
    2. It preserves the zero/non-zero structure of the data
    3. It allows for matrix completion methods to be applied

    The geometric mean is computed only over non-zero values in each
    composition, making it "robust" to the presence of zeros.

    References
    ----------
    .. [1] Martino, C., Morton, J. T., Marotz, C. A., Thompson, L. R., Tripathi, A.,
       Knight, R., & Zengler, K. (2019). A novel sparse compositional technique reveals
       microbial perturbations. MSystems, 4(1), 10-1128.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import rclr
    >>> x = np.array([[1, 2, 0, 4],
    ...               [0, 3, 3, 0],
    ...               [2, 2, 2, 2]])
    >>> result = rclr(x)
    >>> np.round(result, 3)
    array([[-0.693,  0.   ,    nan,  0.693],
           [   nan,  0.   ,  0.   ,    nan],
           [ 0.   ,  0.   ,  0.   ,  0.   ]])

    """
    xp, mat = ingest_array(mat)
    if validate:
        _check_composition(xp, mat, allnum=False)
    return _rclr(xp, mat, axis)


def _rclr(xp: ModuleType, mat: StdArray, axis: int) -> StdArray:
    """Perform rclr transform."""
    float_dtype = xp.float64
    mat_float = xp.asarray(mat, dtype=float_dtype)

    # Track which values were observed in the ORIGINAL input
    observed_mask = (mat_float > 0) & ~xp.isnan(mat_float)

    # Normalize to closure (using _nansum internally)
    closed = _closure(xp, mat_float, axis)
    closed_safe = xp.where(
        closed > 0, closed, xp.asarray(float("nan"), dtype=float_dtype)
    )

    # Take log (will give -inf for zeros, NaN for NaN)
    log_closed = xp.log(closed_safe)

    # Count observed values from ORIGINAL mask
    n_observed = xp.sum(observed_mask.astype(float_dtype), axis=axis, keepdims=True)

    # Sum logs for observed values only
    log_masked = xp.where(observed_mask, log_closed, xp.asarray(0.0, dtype=float_dtype))
    log_sum = xp.sum(log_masked, axis=axis, keepdims=True)

    # Geometric mean
    n_observed_safe = xp.where(
        n_observed > 0, n_observed, xp.asarray(1.0, dtype=float_dtype)
    )
    geo_mean_log = log_sum / n_observed_safe

    # Center by geometric mean
    result = log_closed - geo_mean_log

    # Replace non-observed with NaN
    result = xp.where(
        observed_mask, result, xp.asarray(float("nan"), dtype=float_dtype)
    )

    return result


@params_aliased([("validate", "check", "0.7.0", True)])
def ilr(
    mat: ArrayLike,
    basis: ArrayLike | None = None,
    axis: int = -1,
    validate: bool = True,
) -> StdArray:
    r"""Perform isometric log ratio (ILR) transformation.

    This function transforms compositions from Aitchison simplex to the real
    space. The :math:`ilr` transform is both an isometry, and an isomorphism
    defined on the following spaces:

    .. math::
        ilr: S^D \rightarrow \mathbb{R}^{D-1}

    The ilr transformation is defined as follows:

    .. math::
        ilr(x) =
        [\langle x, e_1 \rangle_a, \ldots, \langle x, e_{D-1} \rangle_a]

    where :math:`[e_1,\ldots,e_{D-1}]` is an orthonormal basis in the simplex.

    If an orthornormal basis isn't specified, the J. J. Egozcue orthonormal basis
    derived from Gram-Schmidt orthogonalization [1]_ will be used by default.

    .. versionchanged:: 0.7.0
        The function now works on any dimension in arrays of any number of dimensions.

    Parameters
    ----------
    mat : array_like of shape (..., n_components, ...)
        A matrix of positive proportions.
    basis : ndarray or sparse matrix, optional
        Orthonormal basis for Aitchison simplex. Defaults to J. J. Egozcue
        orthonormal basis.
    axis : int, optional
        Axis along which ILR transformation will be performed. That is, each vector
        along this axis is considered as a composition. Default is the last axis (-1).

        .. versionadded:: 0.7.0

    validate : bool, default True
        Check if i) the matrix is compositional, ii) the basis is orthonormal,
        2-dimensional, and the dimensions are matched.

    Returns
    -------
    ndarray of shape (..., n_components - 1,...)
        ILR-transformed matrix.

    See Also
    --------
    ilr_inv

    Notes
    -----
    If the ``basis`` parameter is specified, it is expected to be a basis in
    the Aitchison simplex. If there are :math:`D - 1` elements specified in
    ``mat``, then the dimensions of the basis needs be :math:`(D-1) \times D`,
    where rows represent basis vectors, and the columns represent proportions.

    References
    ----------
    .. [1] Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barcelo-Vidal,
       C. (2003). Isometric logratio transformations for compositional data analysis.
       Mathematical geology, 35(3), 279-300.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import ilr
    >>> x = np.array([.1, .3, .4, .2])
    >>> ilr(x)
    array([-0.7768362 , -0.68339802,  0.11704769])

    """
    xp, mat = ingest_array(mat)
    if validate:
        _check_composition(xp, mat, nozero=True)
    N = mat.shape[axis]
    if basis is None:
        # NOTE: acc.device(mat) would be nicer
        basis = xp.asarray(
            _gram_schmidt_basis(N), device=mat.device, dtype=xp.float64
        )  # dimension (N-1) x N
    else:
        xp_, basis = ingest_array(basis)
        if validate:
            # the following maybe redundant
            if basis.ndim != 2:
                raise ValueError(
                    f"Basis needs to be a 2-D matrix, not a {basis.ndim}-D matrix."
                )
            _check_basis(xp_, basis, orthonormal=True, subspace_dim=N - 1)
            basis = xp.asarray(basis, device=mat.device, dtype=xp.float64)
    axis %= mat.ndim
    return _ilr(xp, mat, basis, axis)


def _swap_axis(ndim, axis):
    """Create a list of axis indices with one axis swapped with the last axis."""
    res = list(range(ndim))
    res[axis] = ndim - 1
    res[ndim - 1] = axis
    return res


def _ilr(xp: ModuleType, mat: StdArray, basis: StdArray, axis: int) -> StdArray:
    """Perform ILR transform."""
    mat = _clr(xp, mat, axis)
    # tensordot return's shape consists of the non-contracted axes (dimensions) of
    # the first array x1, followed by the non-contracted axes (dimensions) of the
    # second array x2
    prod = xp.tensordot(mat, basis, axes=([axis], [1]))
    return xp.permute_dims(prod, axes=_swap_axis(mat.ndim, axis))


@params_aliased([("validate", "check", "0.7.0", True)])
def ilr_inv(
    mat: ArrayLike,
    basis: ArrayLike | None = None,
    axis: int = -1,
    validate: bool = True,
) -> StdArray:
    r"""Perform inverse isometric log ratio (ILR) transformation.

    This function transforms compositions from the real space to Aitchison
    geometry. The :math:`ilr^{-1}` transform is both an isometry, and an
    isomorphism defined on the following spaces:

    .. math::
        ilr^{-1}: \mathbb{R}^{D-1} \rightarrow S^D

    The inverse ilr transformation is defined as follows:

    .. math::
        ilr^{-1}(x) = \bigoplus\limits_{i=1}^{D-1} x \odot e_i

    where :math:`[e_1,\ldots, e_{D-1}]` is an orthonormal basis in the simplex.

    If an orthonormal basis isn't specified, the J. J. Egozcue orthonormal basis
    derived from Gram-Schmidt orthogonalization [1]_ will be used by default.

    .. versionchanged:: 0.7.0
        The function now works on any dimension in arrays of any number of dimensions.

    Parameters
    ----------
    mat : array_like of shape (..., n_components - 1, ...)
        A matrix of ILR-transformed data.
    basis : ndarray or sparse matrix, optional
        Orthonormal basis for Aitchison simplex. Defaults to J. J. Egozcue
        orthonormal basis.
    axis : int, optional
        Axis along which ILR transformation will be performed. That is, each vector
        along this axis is considered as a ILR transformed composition data.
        Default is the last axis (-1).

        .. versionadded:: 0.7.0

    validate : bool, default True
        Check to see if basis is orthonormal and dimension matches.

    Returns
    -------
    ndarray of shape (..., n_components, ...)
        Inverse ILR-transformed matrix.

    See Also
    --------
    ilr

    Notes
    -----
    If the ``basis`` parameter is specified, it is expected to be a basis in
    the Aitchison simplex. If there are :math:`D - 1` elements specified in
    ``mat``, then the dimensions of the basis needs be :math:`(D-1) \times D`,
    where rows represent basis vectors, and the columns represent proportions.

    References
    ----------
    .. [1] Egozcue, J. J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barcelo-Vidal,
       C. (2003). Isometric logratio transformations for compositional data analysis.
       Mathematical geology, 35(3), 279-300.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import ilr
    >>> x = np.array([.1, .3, .6,])
    >>> ilr_inv(x)
    array([ 0.34180297,  0.29672718,  0.22054469,  0.14092516])

    """
    xp, mat = ingest_array(mat)
    N = mat.shape[axis] + 1

    if basis is None:
        basis = xp.asarray(
            _gram_schmidt_basis(N), device=mat.device, dtype=xp.float64
        )  # dimension (N-1) x N
    elif validate:
        xp_, basis = ingest_array(basis)
        # the following maybe redundant as the orthonrmal implicitly check 2-d
        if basis.ndim != 2:
            raise ValueError(
                f"Basis needs to be a 2-D matrix, not a {basis.ndim}-D matrix."
            )
        _check_basis(xp_, basis, orthonormal=True, subspace_dim=N - 1)
        basis = xp.asarray(basis, device=mat.device, dtype=xp.float64)
    axis %= mat.ndim
    return _ilr_inv(xp, mat, basis, axis)


def _ilr_inv(xp: ModuleType, mat: StdArray, basis: StdArray, axis: int) -> StdArray:
    """Perform ILR transform."""
    prod = xp.tensordot(mat, basis, axes=([axis], [0]))
    perm = xp.permute_dims(prod, axes=_swap_axis(mat.ndim, axis))
    return _clr_inv(xp, perm, axis)


@params_aliased([("ref_idx", "denominator_idx", "0.7.0", False)])
def alr(
    mat: ArrayLike, ref_idx: int = 0, axis: int = -1, validate: bool = True
) -> StdArray:
    r"""Perform additive log ratio (ALR) transformation.

    This function transforms compositions from a D-part Aitchison simplex to
    a non-isometric real space of D-1 dimensions. The argument
    ``ref_idx`` defines the index of the column used as the reference (a common
    denominator). The :math:`alr` transformed data are amenable to multivariate
    analysis as long as statistics don't involve distances.

    .. math::
        alr: S^D \rightarrow \mathbb{R}^{D-1}

    The alr transformation is defined as follows

    .. math::
        alr(x) = \left[ \ln \frac{x_1}{x_D}, \ldots,
        \ln \frac{x_{D-1}}{x_D} \right]

    where :math:`D` is the index of the part used as the reference.

    .. versionchanged:: 0.7.0
        The function now works on any dimension in arrays of any number of dimensions.

    Parameters
    ----------
    mat : array_like of shape (..., n_components, ...)
        A matrix of proportions.
    ref_idx : int, optional
        Index on the target axis which should be used as the reference composition
        (denominator). Default is 0 (the first position).
    axis : int, optional
        Axis along which ALR transformation will be performed. Each vector along this
        axis is considered as a composition. Default is the last axis (-1).

        .. versionadded:: 0.7.0

    validate: bool, default True
        Check whether the input is positive, whether the mat is 2D.

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray of shape (..., n_components - 1, ...)
        ALR-transformed data projected in a non-isometric real space of
        :math:`D - 1` dimensions for a *D*-parts composition.

    See Also
    --------
    alr_inv

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import alr
    >>> x = np.array([.1, .3, .4, .2])
    >>> alr(x)
    array([ 1.09861229,  1.38629436,  0.69314718])

    """
    xp, mat = ingest_array(mat)
    if validate:
        _check_composition(xp, mat, nozero=True)

    # validate and normalize axis and index
    N = mat.shape[axis]
    if N < 2:
        raise ValueError(f"Dimension {axis} of the input matrix is singleton.")
    axis %= mat.ndim
    if ref_idx < -N or ref_idx >= N:
        raise IndexError(f"Invalid index {ref_idx} on dimension {axis}.")
    ref_idx %= N
    return _alr(xp, mat, ref_idx, axis)


def _alr(xp: ModuleType, mat: StdArray, ref_idx: int, axis: int) -> StdArray:
    # Given that: log(numerator / denominator) = log(numerator) - log(denominator)
    # The following code will perform logarithm on the entire matrix, then subtract
    # denominator from numerator. This is also for numerical stability.
    lmat = xp.log(mat)

    # The following code can be replaced with a single NumPy function call:
    #     numerator_matrix = xp.delete(lmat, ref_idx, axis=axis)
    # However, `delete` is not in the Python array API standard. For compatibility with
    # libraries that don't have `delete`, an arbitrary dimension slicing method is
    # is provided below.
    before: Any = [slice(None)] * mat.ndim
    before[axis] = slice(None, ref_idx)
    before = tuple(before)
    after: Any = [slice(None)] * mat.ndim
    after[axis] = slice(ref_idx + 1, None)
    after = tuple(after)
    numerator_matrix = xp.concat((lmat[before], lmat[after]), axis=axis)

    # The following code can be replaced with a single NumPy function call:
    #     denominator_vector = xp.take(lmat, xp.asarray([ref_idx]), axis=axis)
    # `take` is in the Python array API standard. The following code is to keep the
    # style consistent with the code above.
    column: Any = [slice(None)] * mat.ndim
    column[axis] = slice(ref_idx, ref_idx + 1)
    column = tuple(column)
    denominator_vector = lmat[column]

    return numerator_matrix - denominator_vector


@params_aliased([("ref_idx", "denominator_idx", "0.7.0", False)])
def alr_inv(mat: ArrayLike, ref_idx: int = 0, axis: int = -1) -> StdArray:
    r"""Perform inverse additive log ratio (ALR) transform.

    This function transforms compositions from the non-isometric real space of
    ALRs to Aitchison geometry.

    .. math::
        alr^{-1}: \mathbb{R}^{D-1} \rightarrow S^D

    The inverse ALR transformation is defined as follows:

    .. math::
         alr^{-1}(x) = C[exp([y_1, y_2, ..., y_{D-1}, 0])]

    where :math:`C[x]` is the closure operation defined as

    .. math::
        C[x] = \left[\frac{x_1}{\sum_{i=1}^{D} x_i},\ldots,
                     \frac{x_D}{\sum_{i=1}^{D} x_i} \right]

    for some :math:`D` dimensional real vector :math:`x` and
    :math:`D` is the number of components for every composition.

    .. versionchanged:: 0.7.0
        The function now works on any dimension in arrays of any number of dimensions.

    Parameters
    ----------
    mat : array_like of shape (..., n_components - 1, ...)
        A matrix of ALR-transformed data.
    ref_idx : int, optional
        Index on the target axis where the reference composition (denominator) will be
        inserted. Default is 0 (the first position).
    axis : int, optional
        Axis along which inverse ALR transformation will be performed. Each vector on
        this axis is considered as a CLR-transformed composition. Default is the last
        axis (-1).

        .. versionadded:: 0.7.0

    Returns
    -------
    ndarray of shape (..., n_components, ...)
        Inverse ALR-transformed matrix or vector where rows sum to 1.

    See Also
    --------
    alr

    Notes
    -----
    The output of ``alr_inv`` is guaranteed to have each composition sum to 1. But this
    property isn't required for the input for ``alr``. Therefore, ``alr_inv`` does not
    completely invert ``alr``. Instead, ``alr_inv(clr(mat))`` and ``closure(mat)`` are
    equal.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import alr, alr_inv
    >>> x = np.array([.1, .3, .4, .2])
    >>> alr_inv(alr(x))
    array([ 0.1,  0.3,  0.4,  0.2])

    """
    xp, mat = ingest_array(mat)

    # validate and normalize axis and index
    N = mat.shape[axis] + 1
    if N < 2:
        raise ValueError(f"Dimension {axis} of the input matrix has zero length.")
    axis %= mat.ndim
    if ref_idx < -N or ref_idx >= N:
        raise IndexError(f"Invalid index {ref_idx} on dimension {axis}.")
    ref_idx %= N
    return _alr_inv(xp, mat, ref_idx, axis)


def _alr_inv(xp: ModuleType, mat: StdArray, ref_idx: int, axis: int) -> StdArray:
    # The following code can be replaced with a single NumPy function call.
    #     comp = xp.insert(emat, ref_idx, 1.0, axis=axis)
    # However, `insert` is not in the Python array API standard. For compatibility with
    # libraries that don't have `insert`, an arbitrary dimension slicing method is
    # is provided below.
    before: Any = [slice(None)] * mat.ndim
    before[axis] = slice(None, ref_idx)
    before = tuple(before)
    after: Any = [slice(None)] * mat.ndim
    after[axis] = slice(ref_idx, None)
    after = tuple(after)
    shape: Any = list(mat.shape)
    shape[axis] = 1
    shape = tuple(shape)
    zeros = xp.zeros(shape, dtype=mat.dtype, device=mat.device)
    comp = xp.concat((mat[before], zeros, mat[after]), axis=axis)
    comp = xp.exp(comp - xp.max(comp, axis=axis, keepdims=True))

    return _closure(xp, comp, axis)


def centralize(mat: ArrayLike) -> StdArray:
    r"""Center data around its geometric average.

    Parameters
    ----------
    mat : array_like of shape (n_compositions, n_components)
        A matrix of proportions.

    Returns
    -------
    ndarray of shape (n_compositions, n_components)
        Centered composition matrix.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import centralize
    >>> X = np.array([[.1, .3, .4, .2], [.2, .2, .2, .4]])
    >>> centralize(X)
    array([[ 0.17445763,  0.30216948,  0.34891526,  0.17445763],
           [ 0.32495488,  0.18761279,  0.16247744,  0.32495488]])

    """
    from scipy.stats import gmean

    mat = closure(mat)
    cen = gmean(mat, axis=0)
    return perturb_inv(mat, cen)


def _vlr(x, y, ddof):
    r"""Calculate variance log ratio.

    Parameters
    ----------
    x : array_like of shape (n_components,)
        A vector of proportions.
    y : array_like of shape (n_components,)
        A vector of proportions.
    ddof : int
        Degrees of freedom.

    Returns
    -------
    float
        Variance log ratio value.

    """
    # Log transformation
    x = np.log(x)
    y = np.log(y)

    # Variance log ratio
    return np.var(x - y, ddof=ddof)


def _robust_vlr(x, y, ddof):
    r"""Calculate variance log ratio while masking zeros.

    Parameters
    ----------
    x : array_like of shape (n_components,)
        A vector of proportions.
    y : array_like of shape (n_components,)
        A vector of proportions.
    ddof : int
        Degrees of freedom.

    Returns
    -------
    float
        Variance log ratio value.

    """
    # Mask zeros
    x = np.ma.masked_array(x, mask=x == 0)
    y = np.ma.masked_array(y, mask=y == 0)

    # Log transformation
    x = np.ma.log(x)
    y = np.ma.log(y)

    # Variance log ratio
    return np.ma.var(x - y, ddof=ddof)


def vlr(x, y, ddof=1, robust=False):
    r"""Calculate variance log ratio.

    Parameters
    ----------
    x : array_like of shape (n_components,)
        A vector of proportions.
    y : array_like of shape (n_components,)
        A vector of proportions.
    ddof : int
        Degrees of freedom.
    robust : bool
        Whether to mask zeros at the cost of performance.

    Returns
    -------
    float
        Variance log ratio value.

    Notes
    -----
    Variance log ratio was described in [1]_ and [2]_.

    References
    ----------
    .. [1] V. Lovell D, Pawlowsky-Glahn V, Egozcue JJ, Marguerat S,
           Bähler J (2015) Proportionality: A Valid Alternative to
           Correlation for Relative Data. PLoS Comput Biol 11(3): e1004075.
           https://doi.org/10.1371/journal.pcbi.1004075
    .. [2] Erb, I., Notredame, C.
           How should we measure proportionality on relative gene
           expression data?. Theory Biosci. 135, 21-36 (2016).
           https://doi.org/10.1007/s12064-015-0220-8

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import vlr
    >>> x = np.exp([1, 2, 3])
    >>> y = np.exp([2, 3, 4])
    >>> vlr(x, y)  # no zeros
    0.0

    """
    # Convert array_like to numpy array
    x = closure(x)
    y = closure(y)

    # Set up input and parameters
    kwargs = {
        "x": x,
        "y": y,
        "ddof": ddof,
    }

    # Run backend function
    if robust:
        return _robust_vlr(**kwargs)
    else:
        return _vlr(**kwargs)


def _pairwise_vlr(mat, ddof):
    r"""Perform pairwise variance log ratio transformation.

    Parameters
    ----------
    mat : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    ddof : int
        Degrees of freedom.

    Returns
    -------
    ndarray of shape (n_compositions, n_compositions)
        Distance matrix of variance log ratio values.

    """
    # Log Transform
    X_log = np.log(mat)

    # Variance Log Ratio
    covariance = np.cov(X_log.T, ddof=ddof)
    diagonal = np.diagonal(covariance)
    vlr_data = -2 * covariance + diagonal[:, np.newaxis] + diagonal
    return vlr_data


def pairwise_vlr(mat, ids=None, ddof=1, robust=False, validate=True):
    r"""Perform pairwise variance log ratio transformation.

    Parameters
    ----------
    mat : array_like of shape (n_compositions, n_components)
        A matrix of proportions.
    ids : array_like of str of shape (n_components,)
        Component names.
    ddof : int
        Degrees of freedom.
    robust : bool
        Whether to mask zeros at the cost of performance.
    validate : bool
        Whether to validate the distance matrix after construction.

    Returns
    -------
    skbio.DistanceMatrix if validate=True
        Distance matrix of variance log ratio values.
    skbio.DissimilarityMatrix if validate=False
        Dissimilarity matrix of variance log ratio values.

    Notes
    -----
    Pairwise variance log ratio transformation was described in [1]_ and [2]_.

    References
    ----------
    .. [1] V. Lovell D, Pawlowsky-Glahn V, Egozcue JJ, Marguerat S,
           Bähler J (2015) Proportionality: A Valid Alternative to
           Correlation for Relative Data. PLoS Comput Biol 11(3): e1004075.
           https://doi.org/10.1371/journal.pcbi.1004075
    .. [2] Erb, I., Notredame, C.
           How should we measure proportionality on relative gene
           expression data?. Theory Biosci. 135, 21-36 (2016).
           https://doi.org/10.1007/s12064-015-0220-8

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.stats.composition import pairwise_vlr
    >>> mat = np.array([np.exp([1, 2, 2]),
    ...                 np.exp([2, 3, 6]),
    ...                 np.exp([2, 3, 12])]).T
    >>> dism = pairwise_vlr(mat)
    >>> dism.redundant_form()
    array([[  0.,   3.,  27.],
           [  3.,   0.,  12.],
           [ 27.,  12.,   0.]])

    """
    from skbio.stats.distance import DistanceMatrix

    # Mask zeros
    mat = closure(mat.astype(np.float64))

    # Set up input and parameters
    kwargs = {
        "mat": mat,
        "ddof": ddof,
    }

    # Variance log ratio
    if robust:
        raise NotImplementedError("Pairwise version of robust VLR is not implemented.")
    else:
        vlr_data = _pairwise_vlr(**kwargs)

    # Return distance matrix
    if validate:
        vlr_data = 0.5 * (vlr_data + vlr_data.T)
        return DistanceMatrix(vlr_data, ids=ids)

    # Return dissimilarity matrix
    else:
        return DistanceMatrix(vlr_data, ids=ids, validate=False)


def tree_basis(tree):
    r"""Calculate the sparse representation of an ilr basis from a tree.

    This computes an orthonormal basis specified from a bifurcating tree.

    Parameters
    ----------
    tree : skbio.TreeNode
        Input bifurcating tree. Must be strictly bifurcating (i.e. every
        internal node needs to have exactly two children). This is used to
        specify the ilr basis.

    Returns
    -------
    scipy.sparse.coo_matrix
        The ilr basis required to perform the ilr_inv transform. This is also
        known as the sequential binary partition. Note that this matrix is
        represented in clr coordinates.
    list of str
        List of tree node names indicating the ordering in the basis.

    Raises
    ------
    ValueError
        If the tree doesn't contain two branches.

    Examples
    --------
    >>> from skbio import TreeNode
    >>> tree = u"((b,c)a, d)root;"
    >>> t = TreeNode.read([tree])
    >>> basis, nodes = tree_basis(t)
    >>> basis.toarray()
    array([[-0.40824829, -0.40824829,  0.81649658],
           [-0.70710678,  0.70710678,  0.        ]])

    """
    from scipy.sparse import coo_matrix

    # Specifies which child is numerator and denominator
    # within any given node in a tree.
    NUMERATOR = 1
    DENOMINATOR = 0

    # this is inspired by @wasade in
    # https://github.com/biocore/gneiss/pull/8
    t = tree.copy()
    D = len(list(tree.tips()))

    # calculate number of tips under each node
    for n in t.postorder(include_self=True):
        if n.is_tip():
            n._tip_count = 1
        else:
            if len(n.children) == 2:
                left, right = (
                    n.children[NUMERATOR],
                    n.children[DENOMINATOR],
                )
            else:
                raise ValueError("Not a strictly bifurcating tree.")
            n._tip_count = left._tip_count + right._tip_count

    # calculate k, r, s, t coordinate for each node
    left, right = (
        t.children[NUMERATOR],
        t.children[DENOMINATOR],
    )
    t._k, t._r, t._s, t._t = 0, left._tip_count, right._tip_count, 0
    for n in t.preorder(include_self=False):
        if n.is_tip():
            n._k, n._r, n._s, n._t = 0, 0, 0, 0

        elif n == n.parent.children[NUMERATOR]:
            n._k = n.parent._k
            n._r = n.children[NUMERATOR]._tip_count
            n._s = n.children[DENOMINATOR]._tip_count
            n._t = n.parent._s + n.parent._t
        elif n == n.parent.children[DENOMINATOR]:
            n._k = n.parent._r + n.parent._k
            n._r = n.children[NUMERATOR]._tip_count
            n._s = n.children[DENOMINATOR]._tip_count
            n._t = n.parent._t
        else:
            raise ValueError("Tree topology is not correct.")

    # navigate through tree to build the basis in a sparse matrix form
    value = []
    row, col = [], []
    nodes = []
    i = 0

    for n in t.levelorder(include_self=True):
        if n.is_tip():
            continue

        for j in range(n._k, n._k + n._r):
            row.append(i)
            # consider tips in reverse order. May want to rethink
            # this orientation in the future.
            col.append(D - 1 - j)
            A = np.sqrt(n._s / (n._r * (n._s + n._r)))

            value.append(A)

        for j in range(n._k + n._r, n._k + n._r + n._s):
            row.append(i)
            col.append(D - 1 - j)
            B = -np.sqrt(n._r / (n._s * (n._s + n._r)))

            value.append(B)
        i += 1
        nodes.append(n.name)

    basis = coo_matrix((value, (row, col)), shape=(D - 1, D))

    return basis, nodes


def _gram_schmidt_basis(n):
    """Build CLR-transformed basis derived from Gram-Schmidt orthogonalization.

    Parameters
    ----------
    n : int
        Dimension of the Aitchison simplex.

    Returns
    -------
    basis : array_like of shape (n - 1, n)
        Basis matrix.

    """
    basis = np.zeros((n, n - 1))
    for j in range(n - 1):
        i = j + 1
        e = np.array([(1 / i)] * i + [-1] + [0] * (n - i - 1)) * np.sqrt(i / (i + 1))
        basis[:, j] = e
    return basis.T


def sbp_basis(sbp):
    r"""Build an orthonormal basis from a sequential binary partition (SBP).

    A SBP is a hierarchical collection of binary divisions of compositional
    parts ([1]_). The child groups are divided again until all groups contain a
    single part. The SBP can be encoded in a :math:`(D - 1) \times D` matrix
    where, for each row, parts can be grouped by -1 and +1 tags, and 0 for
    excluded parts. The *i*-th balance is computed as follows:

    .. math::
        b_i = \sqrt{ \frac{r_i s_i}{r_i+s_i} }
        \ln \left( \frac{g(x_{r_i})}{g(x_{s_i})} \right)

    where :math:`b_i` is the *i*-th balance corresponding to the *i*-th row in
    the SBP, :math:`r_i` and :math:`s_i` and the number of respectively ``+1``
    and ``-1`` labels in the *i*-th row of the SBP and where :math:`g(x) =
    (\prod\limits_{i=1}^{D} x_i)^{1/D}` is the geometric mean of :math:`x`.

    Parameters
    ----------
    sbp : array_like of shape (n_partitions, n_features)
        A contrast matrix, also known as a sequential binary partition, where
        every row represents a partition between two groups of features. A part
        labelled ``+1`` would correspond to that feature being in the numerator
        of the given row partition, a part labelled ``-1`` would correspond to
        features being in the denominator of that given row partition, and
        ``0`` would correspond to features excluded in the row partition.

    Returns
    -------
    ndarray of shape (n_partitions, n_features)
        An orthonormal basis in the Aitchison simplex.

    Notes
    -----
    The ``sbp_basis`` method was derived from the ``gsi.buildilrBase()``
    function implemented in the R package "compositions" [2]_.

    Examples
    --------
    >>> import numpy as np
    >>> sbp = np.array([[1, 1,-1,-1,-1],
    ...                 [1,-1, 0, 0, 0],
    ...                 [0, 0, 1,-1,-1],
    ...                 [0, 0, 0, 1,-1]])
    ...
    >>> sbp_basis(sbp)
    array([[ 0.54772256,  0.54772256, -0.36514837, -0.36514837, -0.36514837],
           [ 0.70710678, -0.70710678,  0.        ,  0.        ,  0.        ],
           [ 0.        ,  0.        ,  0.81649658, -0.40824829, -0.40824829],
           [ 0.        ,  0.        ,  0.        ,  0.70710678, -0.70710678]])

    References
    ----------
    .. [1] Parent, S.É., Parent, L.E., Egozcue, J.J., Rozane, D.E.,
       Hernandes, A., Lapointe, L., Hébert-Gentile, V., Naess, K.,
       Marchand, S., Lafond, J., Mattos, D., Barlow, P., Natale, W., 2013.
       The plant ionome revisited by the nutrient balance concept.
       Front. Plant Sci. 4, 39.
    .. [2] van den Boogaart, K. Gerald, Tolosana-Delgado, Raimon and Bren,
       Matevz, 2014. `compositions`: Compositional Data Analysis. R package
       version 1.40-1. https://CRAN.R-project.org/package=compositions.

    """
    n_pos = (sbp == 1).sum(axis=1)
    n_neg = (sbp == -1).sum(axis=1)
    psi = np.zeros(sbp.shape)
    for i in range(0, sbp.shape[0]):
        psi[i, :] = sbp[i, :] * np.sqrt(
            (n_neg[i] / n_pos[i]) ** sbp[i, :] / np.sum(np.abs(sbp[i, :]))
        )
    return psi


def _check_basis(
    xp: ModuleType,
    basis: StdArray,
    orthonormal: bool = False,
    subspace_dim: int | None = None,
):
    r"""Check if basis is a valid basis for transformation.

    Parameters
    ----------
    xp : namespace
        The array API compatible namespace corresponding ``basis``.
    basis : array of shape (n_basis, n_components)
        A columns vectors for the basis.
    orthonormal : bool, optional
        If True, basis is required to be orthonormal. Default is False.
    subspace_dim : int, optional
        The dimensions of the subspace that the basis suppose to span,
        when None is give, the n_basis will be used. Default is None.

    Raises
    ------
    ValueError
        If the basis is not matching to the subspace dimension.
    ValueError
        If the basis are not orthonormal.

    """
    xp, basis = ingest_array(basis)
    if basis.ndim < 2:
        basis = basis.reshape(1, -1)
    if subspace_dim is None:
        subspace_dim = len(basis)
    elif len(basis) != subspace_dim:
        n_basis = len(basis)
        msg = f"Number of basis {n_basis} not match to the subspace dim {subspace_dim}."
        raise ValueError(msg)
    if orthonormal:
        eyes = xp.eye(subspace_dim, device=basis.device)
        if not xp.all(xp.abs(basis @ basis.T - eyes) < (1e-4 * eyes + 1e-6)):
            raise ValueError("Basis is not orthonormal.")


register_aliases(modules[__name__])
