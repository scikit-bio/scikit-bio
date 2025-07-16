r"""Composition Statistics (:mod:`skbio.stats.composition`)
=======================================================

.. currentmodule:: skbio.stats.composition

This module provides functions for compositional data analysis.

Many omics datasets are inherently compositional -- meaning that they are best
interpreted as proportions or percentages rather than absolute counts.

Formally, sample :math:`x` is a composition if :math:`\sum_{i=0}^D x_{i} = c`
and :math:`x_{i} > 0`, :math:`1 \leq i \leq D` and :math:`c` is a real-valued
constant and there are :math:`D` components (features) for this composition.
In this module :math:`c=1`. Compositional data can be analyzed using
**Aitchison geometry** [1]_.

However, in this framework, standard real Euclidean operations such as addition
and multiplication no longer apply. Only operations such as perturbation and
power can be used to manipulate this data.

This module allows two styles of manipulation of compositional data.
Compositional data can be analyzed using perturbation and power operations,
which can be useful for simulation studies. The alternative strategy is to
transform compositional data into the real space. Right now, the centre log
ratio transform (clr) and the isometric log ratio transform (ilr) [2]_ can be
used to accomplish this. This transform can be useful for performing standard
statistical methods such as parametric hypothesis testing, regression and more.

The major caveat of using this framework is dealing with zeros. In Aitchison
geometry, only compositions with non-zero components can be considered.
The multiplicative replacement technique [3]_ can be used to substitute these
zeros with small pseudocounts without introducing major distortions to the
data.


Differential abundance
----------------------

Statistical tests for the differential abundance (DA) of components among groups of
compositions.

.. autosummary::
   :toctree:

   ancom
   dirmult_ttest
   dirmult_lme


Arithmetic operations
---------------------

Manipulate compositional data within the Aitchison space.

.. autosummary::
   :toctree:

   centralize
   closure
   inner
   perturb
   perturb_inv
   power


Log-ratio transformation
------------------------

Convert compositional data into log-ratio space to enable subsequent comparison
and statistical analysis.

.. autosummary::
   :toctree:

   alr
   alr_inv
   clr
   clr_inv
   ilr
   ilr_inv

.. note::
   Arithmetic operations and log-ratio transformations support array formats compliant
   with the `Python array API standard <https://data-apis.org/array-api/latest/>`_
   without transition through NumPy. For example, they can directly consume and return
   GPU-resident PyTorch tensors.


Correlation analysis
--------------------

Measure the pairwise relationships of compositional data.

.. autosummary::
   :toctree:

   vlr
   pairwise_vlr


Zero handling
-------------

Replace zero values in compositional data with positive values, which is
necessary prior to logarithmic operations.

.. autosummary::
   :toctree:

   multi_replace


Basis construction
------------------

Generate basis vectors for compositional data via hierarchical partitioning, to
allow for decomposition and transformation, such as ilr transform.

.. autosummary::
   :toctree:

   sbp_basis
   tree_basis


References
----------
.. [1] V. Pawlowsky-Glahn, J. J. Egozcue, R. Tolosana-Delgado (2015),
   Modeling and Analysis of Compositional Data, Wiley, Chichester, UK

.. [2] J. J. Egozcue.,  "Isometric Logratio Transformations for
   Compositional Data Analysis" Mathematical Geology, 35.3 (2003)

.. [3] J. A. Martin-Fernandez,  "Dealing With Zeros and Missing Values in
   Compositional Data Sets Using Nonparametric Imputation",
   Mathematical Geology, 35.3 (2003)

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, TYPE_CHECKING
from sys import modules
from warnings import warn, catch_warnings, simplefilter
import inspect

import numpy as np
import pandas as pd

from skbio.util import get_rng
from skbio.util._decorator import aliased, register_aliases, params_aliased
from skbio.util._array import ingest_array
from skbio.table._tabular import _ingest_table


if TYPE_CHECKING:  # pragma: no cover
    from types import ModuleType
    from skbio.util._typing import ArrayLike, StdArray


def _check_composition(
    xp: "ModuleType",
    mat: "StdArray",
    axis: int = -1,
    nozero: bool = False,
    maxdim: Optional[int] = None,
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
    if not xp.all(xp.isfinite(mat)):
        raise ValueError("Input matrix cannot have infinite or NaN values.")
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


def closure(mat: "ArrayLike", axis: int = -1, validate: bool = True) -> "StdArray":
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


def _closure(xp: "ModuleType", mat: "StdArray", axis: int = -1) -> "StdArray":
    """Perform closure."""
    return mat / xp.sum(mat, axis=axis, keepdims=True)


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


def perturb(x: "ArrayLike", y: "ArrayLike", validate: bool = True) -> "StdArray":
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


def perturb_inv(x: "ArrayLike", y: "ArrayLike", validate: bool = True) -> "StdArray":
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


def power(x: "ArrayLike", a: float, validate: bool = True) -> "StdArray":
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


def inner(x: "ArrayLike", y: "ArrayLike", validate: bool = True) -> "StdArray":
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


def clr(mat: "ArrayLike", axis: int = -1, validate: bool = True) -> "StdArray":
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


def _clr(xp: "ModuleType", mat: "StdArray", axis: int) -> "StdArray":
    """Perform CLR transform."""
    return (lmat := xp.log(mat)) - xp.mean(lmat, axis=axis, keepdims=True)


def clr_inv(mat: "ArrayLike", axis: int = -1, validate: bool = True) -> "StdArray":
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


def _clr_inv(xp: "ModuleType", mat: "StdArray", axis: int) -> "StdArray":
    """Perform inverse CLR transform."""
    # for numerical stability, shift the values < 1
    diff = xp.exp(mat - xp.max(mat, axis=axis, keepdims=True))
    return _closure(xp, diff, axis)


@params_aliased([("validate", "check", "0.7.0", True)])
def ilr(
    mat: "ArrayLike",
    basis: Optional["ArrayLike"] = None,
    axis: int = -1,
    validate: bool = True,
) -> "StdArray":
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


def _ilr(xp: "ModuleType", mat: "StdArray", basis: "StdArray", axis: int) -> "StdArray":
    """Perform ILR transform."""
    mat = _clr(xp, mat, axis)
    # tensordot return's shape consists of the non-contracted axes (dimensions) of
    # the first array x1, followed by the non-contracted axes (dimensions) of the
    # second array x2
    prod = xp.tensordot(mat, basis, axes=([axis], [1]))
    return xp.permute_dims(prod, axes=_swap_axis(mat.ndim, axis))


@params_aliased([("validate", "check", "0.7.0", True)])
def ilr_inv(
    mat: "ArrayLike",
    basis: Optional["ArrayLike"] = None,
    axis: int = -1,
    validate: bool = True,
) -> "StdArray":
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


def _ilr_inv(
    xp: "ModuleType", mat: "StdArray", basis: "StdArray", axis: int
) -> "StdArray":
    """Perform ILR transform."""
    prod = xp.tensordot(mat, basis, axes=([axis], [0]))
    perm = xp.permute_dims(prod, axes=_swap_axis(mat.ndim, axis))
    return _clr_inv(xp, perm, axis)


@params_aliased([("ref_idx", "denominator_idx", "0.7.0", False)])
def alr(
    mat: "ArrayLike", ref_idx: int = 0, axis: int = -1, validate: bool = True
) -> "StdArray":
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


def _alr(xp: "ModuleType", mat: "StdArray", ref_idx: int, axis: int) -> "StdArray":
    # Given that: log(numerator / denominator) = log(numerator) - log(denominator)
    # The following code will perform logarithm on the entire matrix, then subtract
    # denominator from numerator. This is also for numerical stability.
    lmat = xp.log(mat)

    # The following code can be replaced with a single NumPy function call:
    #     numerator_matrix = xp.delete(lmat, ref_idx, axis=axis)
    # However, `delete` is not in the Python array API standard. For compatibility with
    # libraries that don't have `delete`, an arbitrary dimension slicing method is
    # is provided below.
    before = [slice(None)] * mat.ndim
    before[axis] = slice(None, ref_idx)
    before = tuple(before)
    after = [slice(None)] * mat.ndim
    after[axis] = slice(ref_idx + 1, None)
    after = tuple(after)
    numerator_matrix = xp.concat((lmat[before], lmat[after]), axis=axis)

    # The following code can be replaced with a single NumPy function call:
    #     denominator_vector = xp.take(lmat, xp.asarray([ref_idx]), axis=axis)
    # `take` is in the Python array API standard. The following code is to keep the
    # style consistent with the code above.
    column = [slice(None)] * mat.ndim
    column[axis] = slice(ref_idx, ref_idx + 1)
    column = tuple(column)
    denominator_vector = lmat[column]

    return numerator_matrix - denominator_vector


@params_aliased([("ref_idx", "denominator_idx", "0.7.0", False)])
def alr_inv(mat: "ArrayLike", ref_idx: int = 0, axis: int = -1) -> "StdArray":
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


def _alr_inv(xp: "ModuleType", mat: "StdArray", ref_idx: int, axis: int) -> "StdArray":
    # The following code can be replaced with a single NumPy function call.
    #     comp = xp.insert(emat, ref_idx, 1.0, axis=axis)
    # However, `insert` is not in the Python array API standard. For compatibility with
    # libraries that don't have `insert`, an arbitrary dimension slicing method is
    # is provided below.
    before = [slice(None)] * mat.ndim
    before[axis] = slice(None, ref_idx)
    before = tuple(before)
    after = [slice(None)] * mat.ndim
    after[axis] = slice(ref_idx, None)
    after = tuple(after)
    shape = list(mat.shape)
    shape[axis] = 1
    shape = tuple(shape)
    zeros = xp.zeros(shape, dtype=mat.dtype, device=mat.device)
    comp = xp.concat((mat[before], zeros, mat[after]), axis=axis)
    comp = xp.exp(comp - xp.max(comp, axis=axis, keepdims=True))

    return _closure(xp, comp, axis)


def centralize(mat: "ArrayLike") -> "StdArray":
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


def _check_p_adjust(name):
    r"""Construct a p-value correction function based on the method name.

    Parameters
    ----------
    name : str
        The name of the p-value correction method. This should match one of the
        method names available in :func:`statsmodels.stats.multitest.multipletests`.

    Returns
    -------
    callable, optional
        Function to correct p-values.

    """
    if name is None:
        return

    from statsmodels.stats.multitest import multipletests as sm_multipletests

    name_ = name.lower()

    # Original options are kept for backwards compatibility
    if name_ in ("holm", "holm-bonferroni"):
        name_ = "holm"
    if name_ in ("bh", "fdr_bh", "benjamini-hochberg"):
        name_ = "fdr_bh"

    def func(pvals):
        r"""Correct p-values for multiple testing problems.

        Parameters
        ----------
        pvals : ndarray of shape (n_tests,)
            Original p-values.

        Returns
        -------
        qvals : ndarray of shape (n_tests,)
            Corrected p-values.

        """
        try:
            res = sm_multipletests(pvals, alpha=0.05, method=name_)
        except ValueError as e:
            if "method not recognized" in str(e):
                raise ValueError(f'"{name}" is not an available FDR correction method.')
            else:
                raise ValueError(
                    f"Cannot perform FDR correction using the {name} method."
                )
        else:
            return res[1]

    return func


def _check_grouping(grouping, matrix, samples=None):
    """Format grouping for differential abundance analysis.

    Parameters
    ----------
    grouping : 1-D array_like
        Vector indicating the assignment of samples to groups. For example,
        these could be strings or integers denoting which group a sample
        belongs to.
    matrix : ndarray of shape (n_samples, n_features)
        Data matrix.
    samples : array_like of shape (n_samples,), optional
        Sample IDs.

    Returns
    -------
    groups : ndarray of (n_groups,)
        Class names.
    labels : ndarray of (n_samples,)
        Class indices by sample.

    Notes
    -----
    If `grouping` is indexed and `samples` is provided, `grouping` will be filtered and
    reordered to match `samples`. Otherwise, `grouping` and `matrix` must have the same
    length, with the assumption that samples are in the same order.

    """
    # match sample IDs
    if samples is not None and isinstance(grouping, pd.Series):
        try:
            grouping = grouping.loc[samples]
        except KeyError:
            raise ValueError(
                "`table` contains sample IDs that are absent in `grouping`."
            )
        else:
            grouping = grouping.to_numpy()

    # match lengths
    else:
        grouping = np.asarray(grouping)

        if grouping.ndim != 1:
            raise ValueError("`grouping` must be convertible to a 1-D vector.")

        if matrix.shape[0] != grouping.shape[0]:
            raise ValueError(
                "Sample counts in `table` and `grouping` are not consistent."
            )

    # The following code achieves what `pd.isnull` does with NumPy.
    null_errmsg = "Cannot handle missing values in `grouping`."
    if np.isdtype(grouping.dtype, "numeric"):
        if np.isnan(grouping).any():
            raise ValueError(null_errmsg)
    else:
        if (grouping != grouping).any() or np.equal(grouping, None).any():
            raise ValueError(null_errmsg)

    return np.unique(grouping, return_inverse=True)


def _check_trt_ref_groups(treatment, reference, groups, labels):
    """Extract treatment and reference group indices.

    Parameters
    ----------
    treatment : str, int or None
        Treatment group label.
    reference : str, int or None
        Reference group label.
    groups : ndarray of (n_groups,)
        Class names.
    labels : ndarray of (n_samples,)
        Class indices by sample.

    Returns
    -------
    trt_idx : ndarray of (n_samples_in_treatment,)
        Sample indices in the treatment group.
    ref_idx : ndarray of (n_samples_in_reference,)
        Sample indices in the reference group.

    Raises
    ------
    ValueError
        If treatment or reference group is not found.
    ValueError
        If treatment and reference groups are the same.
    ValueError
        If there are less than two groups.

    """
    if len(groups) < 2:
        raise ValueError("There must be at least two groups in grouping.")

    if treatment is not None:
        try:
            (trt_i,) = np.flatnonzero(groups == treatment)
        except ValueError:
            raise ValueError(f"Treatment group {treatment} is not found in grouping.")
    else:
        trt_i = 0
    trt_idx = np.flatnonzero(labels == trt_i)

    if reference is not None:
        try:
            (ref_i,) = np.flatnonzero(groups == reference)
        except ValueError:
            raise ValueError(f"Reference group {reference} is not found in grouping.")
        if trt_i == ref_i:
            raise ValueError("Treatment and reference groups must not be identical.")
        ref_idx = np.flatnonzero(labels == ref_i)
    else:
        ref_idx = np.flatnonzero(labels != trt_i)

    return trt_idx, ref_idx


def _check_metadata(metadata, matrix, samples=None):
    """Format metadata for differential abundance analysis.

    Parameters
    ----------
    metadata : dataframe_like
        Metadata table.
    matrix : ndarray of shape (n_samples, n_features)
        Data matrix.
    samples : array_like of shape (n_samples,), optional
        Sample IDs.

    Returns
    -------
    pd.DataFrame
        Validated metadata table.

    Notes
    -----
    This function resembles `_check_grouping`.

    """
    if not isinstance(metadata, pd.DataFrame):
        try:
            metadata = pd.DataFrame(metadata)
        except Exception:
            raise TypeError(
                "Metadata must be a pandas DataFrame, or a data structure that can be "
                "converted into a pandas DataFrame, such as a NumPy structured or rec "
                "array, or a dictionary."
            )

    # match lengths
    if samples is None or isinstance(metadata.index, pd.RangeIndex):
        if matrix.shape[0] != metadata.shape[0]:
            raise ValueError("Sample counts in table and metadata are not consistent.")

    # match sample IDs
    else:
        if not metadata.index.equals(pd.Index(samples)):
            try:
                metadata = metadata.loc[samples]
            except KeyError:
                raise ValueError(
                    "Metadata contains sample IDs that are absent in the table."
                )

    if metadata.isnull().values.any():
        raise ValueError("Cannot handle missing values in metadata.")

    return metadata


def _check_sig_test(test, n_groups=None):
    """Validate significance test.

    Parameters
    ----------
    test : str or callable
        Statistical testing function or its name.
    n_groups : int, optional
        Number of sample groups.

    Returns
    -------
    callable
        Statistical testing function.

    """
    if isinstance(test, str):
        import scipy.stats

        try:
            func = getattr(scipy.stats, test)
        except AttributeError:
            raise ValueError(f'Function "{test}" does not exist under scipy.stats.')
    else:
        if not callable(test):
            raise TypeError("`sig_test` must be a function or a string.")
        func = test
        test = test.__name__

    if n_groups is not None:
        sig = inspect.signature(func)
        param = next(iter(sig.parameters.values()))

        is_multi = param.kind is inspect.Parameter.VAR_POSITIONAL
        if n_groups > 2 and not is_multi:
            raise ValueError(
                f'"{test}" is a two-way statistical test whereas {n_groups} sample '
                "groups were provided."
            )

    return func


@params_aliased(
    [
        ("p_adjust", "multiple_comparisons_correction", "0.6.0", True),
        ("sig_test", "significance_test", "0.7.0", True),
    ]
)
def ancom(
    table,
    grouping,
    alpha=0.05,
    tau=0.02,
    theta=0.1,
    p_adjust="holm",
    sig_test="f_oneway",
    percentiles=None,
):
    r"""Perform differential abundance test using ANCOM.

    Analysis of composition of microbiomes (ANCOM) is done by calculating
    pairwise log ratios between all features and performing a significance
    test to determine if there is a significant difference in feature ratios
    with respect to the variable of interest.

    In an experiment with only two treatments, this tests the following
    hypothesis for feature :math:`i`:

    .. math::

        H_{0i}: \mathbb{E}[\ln(u_i^{(1)})] = \mathbb{E}[\ln(u_i^{(2)})]

    where :math:`u_i^{(1)}` is the mean abundance for feature :math:`i` in the
    first group and :math:`u_i^{(2)}` is the mean abundance for feature
    :math:`i` in the second group.

    .. versionchanged:: 0.7.0
        Computational efficiency significantly improved.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Matrix of strictly positive values (i.e. counts or proportions). See
        :ref:`supported formats <table_like>`.

        .. note::
            If the table contains zero values, one should add a pseudocount or apply
            :func:`multi_replace` to convert all values into positive numbers.

    grouping : pd.Series or 1-D array_like
        Vector indicating the assignment of samples to groups. These could be strings
        or integers denoting which group a sample belongs to. If it is a pandas Series
        and the table contains sample IDs, its index will be filtered and reordered to
        match the sample IDs. Otherwise, it must be the same length as the samples in
        the table.
    alpha : float, optional
        Significance level for each of the statistical tests. This can can be
        anywhere between 0 and 1 exclusive.
    tau : float, optional
        A constant used to determine an appropriate cutoff. A value close to
        zero indicates a conservative cutoff. This can can be anywhere between
        0 and 1 exclusive.
    theta : float, optional
        Lower bound for the proportion for the *W*-statistic. If all *W*-
        statistics are lower than theta, then no features will be detected to
        be significantly different. This can can be anywhere between 0 and 1
        exclusive.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are Holm-
        Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-
        Hochberg ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported
        by statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.
    sig_test : str or callable, optional
        A function to test for significance between classes. It must be able to
        accept at least two vectors of floats and returns a test statistic and
        a *p*-value. Functions under ``scipy.stats`` can be directly specified
        by name. The default is one-way ANOVA ("f_oneway").

        .. versionchanged:: 0.7.0
            Test funcion must accept 2-D arrays as input, perform batch testing, and
            return 1-D arrays. SciPy functions have this capability. Custom functions
            may need modification.

        .. versionchanged:: 0.6.0
            Accepts test names in addition to functions.

    percentiles : iterable of floats, optional
        Percentile abundances to return for each feature in each group. By
        default, will return the minimum, 25th percentile, median, 75th
        percentile, and maximum abundances for each feature in each group.

    Returns
    -------
    pd.DataFrame
        A table of features, their *W*-statistics and whether the null
        hypothesis is rejected.

        - ``W``: *W*-statistic, or the number of features that the current
          feature is tested to be significantly different against.

        - ``Signif``: Whether the feature is significantly differentially
          abundant across groups (``True``) or not (``False``).

        .. versionchanged:: 0.7.0
            Renamed ``Reject null hypothesis`` as ``Signif``.

    pd.DataFrame
        A table of features and their percentile abundances in each group. If
        ``percentiles`` is empty, this will be an empty ``pd.DataFrame``. The
        rows in this object will be features, and the columns will be a
        multi-index where the first index is the percentile, and the second
        index is the group.

    See Also
    --------
    multi_replace
    scipy.stats.ttest_ind
    scipy.stats.f_oneway
    scipy.stats.wilcoxon
    scipy.stats.kruskal

    Notes
    -----
    The developers of ANCOM recommend the following significance tests ([1]_,
    Supplementary File 1, top of page 11):

    - If there are two groups, use the standard parametric *t*-test
      (``ttest_ind``) or the non-parametric Mann-Whitney rank test
      (``mannwhitneyu``).

    - For paired samples, use the parametric paired *t*-test (``ttest_rel``) or
      the non-parametric Wilcoxon signed-rank test (``wilcoxon``).

    - If there are more than two groups, use the parametric one-way ANOVA
      (``f_oneway``) or the non-parametric Kruskal-Wallis test (``kruskal``).

    - If there are multiple measurements obtained from the individuals, use a
      Friedman test (``friedmanchisquare``).

    Because one-way ANOVA is equivalent to the standard *t*-test when the
    number of groups is two, we default to ``f_oneway`` here, which can be used
    when there are two or more groups.

    Users should refer to the documentation of these tests in SciPy to
    understand the assumptions made by each test.

    This method cannot handle any zero counts as input, since the logarithm
    of zero cannot be computed.  While this is an unsolved problem, many
    studies, including [1]_, have shown promising results by adding
    pseudocounts to all values in the matrix. In [1]_, a pseudocount of 0.001
    was used, though the authors note that a pseudocount of 1.0 may also be
    useful. Zero counts can also be addressed using the ``multi_replace`` method.

    References
    ----------
    .. [1] Mandal et al. "Analysis of composition of microbiomes: a novel
       method for studying microbial composition", Microbial Ecology in Health
       & Disease, (2015), 26.

    Examples
    --------
    >>> from skbio.stats.composition import ancom
    >>> import pandas as pd

    Let's load in a DataFrame with six samples and seven features (e.g., these
    may be bacterial taxa):

    >>> table = pd.DataFrame([[12, 11, 10, 10, 10, 10, 10],
    ...                       [9,  11, 12, 10, 10, 10, 10],
    ...                       [1,  11, 10, 11, 10, 5,  9],
    ...                       [22, 21, 9,  10, 10, 10, 10],
    ...                       [20, 22, 10, 10, 13, 10, 10],
    ...                       [23, 21, 14, 10, 10, 10, 10]],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'],
    ...                      columns=['b1', 'b2', 'b3', 'b4', 'b5', 'b6',
    ...                               'b7'])

    Then create a grouping vector. In this example, there is a treatment group
    and a placebo group.

    >>> grouping = pd.Series(['treatment', 'treatment', 'treatment',
    ...                       'placebo', 'placebo', 'placebo'],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'])

    Now run ``ancom`` to determine if there are any features that are
    significantly different in abundance between the treatment and the placebo
    groups. The first DataFrame that is returned contains the ANCOM test
    results, and the second contains the percentile abundance data for each
    feature in each group.

    >>> ancom_df, percentile_df = ancom(table, grouping)
    >>> ancom_df['W'] # doctest: +ELLIPSIS
    b1    0
    b2    4
    b3    0
    b4    1
    b5    1
    b6    0
    b7    1
    Name: W, dtype: ...

    The *W*-statistic is the number of features that a single feature is tested
    to be significantly different against. In this scenario, ``b2`` was
    detected to have significantly different abundances compared to four of the
    other features. To summarize the results from the *W*-statistic, let's take
    a look at the results from the hypothesis test. The ``Signif`` column in the
    table indicates whether the null hypothesis was rejected, and that a feature
    was therefore observed to be differentially abundant across the groups.

    >>> ancom_df['Signif']
    b1    False
    b2     True
    b3    False
    b4    False
    b5    False
    b6    False
    b7    False
    Name: Signif, dtype: bool

    From this we can conclude that only ``b2`` was significantly different in
    abundance between the treatment and the placebo. We still don't know, for
    example, in which group ``b2`` was more abundant. We therefore may next be
    interested in comparing the abundance of ``b2`` across the two groups.
    We can do that using the second DataFrame that was returned. Here we
    compare the median (50th percentile) abundance of ``b2`` in the treatment
    and placebo groups:

    >>> percentile_df[50.0].loc['b2']
    Group
    placebo      21.0
    treatment    11.0
    Name: b2, dtype: float64

    We can also look at a full five-number summary for ``b2`` in the treatment
    and placebo groups:

    >>> percentile_df.loc['b2'] # doctest: +NORMALIZE_WHITESPACE
    Percentile  Group
    0.0         placebo      21.0
    25.0        placebo      21.0
    50.0        placebo      21.0
    75.0        placebo      21.5
    100.0       placebo      22.0
    0.0         treatment    11.0
    25.0        treatment    11.0
    50.0        treatment    11.0
    75.0        treatment    11.0
    100.0       treatment    11.0
    Name: b2, dtype: float64

    Taken together, these data tell us that ``b2`` is present in significantly
    higher abundance in the placebo group samples than in the treatment group
    samples.

    """
    matrix, samples, features = _ingest_table(table)

    groups, labels = _check_grouping(grouping, matrix, samples)

    _check_composition(np, matrix, nozero=True)

    # validate parameters
    if not 0 < alpha < 1:
        raise ValueError("`alpha`=%f is not within 0 and 1." % alpha)

    if not 0 < tau < 1:
        raise ValueError("`tau`=%f is not within 0 and 1." % tau)

    if not 0 < theta < 1:
        raise ValueError("`theta`=%f is not within 0 and 1." % theta)

    # validate percentiles
    if percentiles is None:
        percentiles = np.arange(0, 125, 25.0)
    else:
        if not isinstance(percentiles, np.ndarray):
            percentiles = np.fromiter(percentiles, dtype=float)
        if (percentiles < 0.0).any() or (percentiles > 100.0).any():
            raise ValueError("Percentiles must be in the range [0, 100].")
        n_pcts = len(percentiles)
        percentiles = np.unique(percentiles)
        if percentiles.size != n_pcts:
            raise ValueError("Percentile values must be unique.")

    n_groups = len(groups)
    if n_groups == len(labels):
        raise ValueError(
            "All values in `grouping` are unique. This method cannot "
            "operate on a grouping vector with only unique values (e.g., "
            "there are no 'within' variance because each group of samples "
            "contains only a single sample)."
        )
    elif n_groups == 1:
        raise ValueError(
            "All values the `grouping` are the same. This method cannot "
            "operate on a grouping vector with only a single group of samples"
            "(e.g., there are no 'between' variance because there is only a "
            "single group)."
        )

    # validate significance test
    if sig_test is None:
        sig_test = "f_oneway"
    test_f = _check_sig_test(sig_test, n_groups)

    # compare log ratios
    pval_mat = _log_compare(matrix, labels, n_groups, test_f)

    # correct for multiple testing problem
    if p_adjust is not None:
        func = _check_p_adjust(p_adjust)
        pval_mat = np.apply_along_axis(func, 1, pval_mat)

    np.fill_diagonal(pval_mat, 1)

    # calculate W-statistics
    n_feats = matrix.shape[1]
    W = (pval_mat < alpha).sum(axis=1)
    c_start = W.max() / n_feats
    if c_start < theta:
        reject = np.zeros_like(W, dtype=bool)
    else:
        # Select appropriate cutoff
        cutoff = c_start - np.linspace(0.05, 0.25, 5)
        prop_cut = (W[:, None] > n_feats * cutoff).mean(axis=0)
        dels = np.abs(prop_cut - np.roll(prop_cut, -1))
        dels[-1] = 0

        if (dels[0] < tau) and (dels[1] < tau) and (dels[2] < tau):
            nu = cutoff[1]
        elif (dels[0] >= tau) and (dels[1] < tau) and (dels[2] < tau):
            nu = cutoff[2]
        elif (dels[1] >= tau) and (dels[2] < tau) and (dels[3] < tau):
            nu = cutoff[3]
        else:
            nu = cutoff[4]
        reject = W >= nu * n_feats

    ancom_df = pd.DataFrame(
        {
            "W": pd.Series(W, index=features),
            "Signif": pd.Series(reject, index=features),
        }
    )

    # calculate percentiles
    if percentiles.size == 0:
        return ancom_df, pd.DataFrame()
    data = []
    columns = []
    for i, group in enumerate(groups):
        feat_dists = matrix[labels == i]
        for percentile in percentiles:
            columns.append((percentile, group))
            data.append(np.percentile(feat_dists, percentile, axis=0))
    columns = pd.MultiIndex.from_tuples(columns, names=["Percentile", "Group"])
    percentile_df = pd.DataFrame(np.asarray(data).T, columns=columns, index=features)
    return ancom_df, percentile_df


def _log_compare(matrix, labels, n, test):
    """Compare pairwise log ratios between sample groups.

    Calculate pairwise log ratios between all features and perform a statistical test
    to determine if there is a significant difference in feature ratios with respect
    to the variable of interest.

    Parameters
    ----------
    matrix : ndarray of shape (n_samples, n_features)
        Data matrix.
    labels : ndarray of shape (n_samples,)
        Group indices (0-indexed, consecutive).
    n : int
        Number of groups.
    test : callable
        Statistical test to run.

    Returns
    -------
    ndarray of shape (n_features, n_features)
        p-value matrix.

    """
    # note: `n` can be simply computed with `labels.max()`. It is supplied instead to
    # save compute.

    # log-transform data
    log_mat = np.log(matrix)

    # divide data by sample group
    grouped = [log_mat[labels == i] for i in range(n)]

    # determine all pairs of feature indices
    m = matrix.shape[1]
    ii, jj = np.triu_indices(m, k=1)

    # calculate all log ratios (pairwise difference of log values)
    log_ratios = [x[:, ii] - x[:, jj] for x in grouped]

    # run statistical test on the 2-D arrays in a vectorized manner
    _, pvals = test(*log_ratios)

    # populate p-value matrix
    pval_mat = np.empty((m, m))
    pval_mat[ii, jj] = pval_mat[jj, ii] = pvals
    np.fill_diagonal(pval_mat, 0)

    return pval_mat


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
    xp: "ModuleType",
    basis: "StdArray",
    orthonormal: bool = False,
    subspace_dim: Optional[int] = None,
):
    r"""Check if basis is a valid basis for transformation.

    Parameters
    ----------
    xp : namespace
        The array API compatible namespace corresponding ``basis``.
    basis : array of shape (n_basis, n_components)
        A columns vetors for the basis.
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


def _dirmult_draw(matrix, rng):
    """Resample data from a Dirichlet-multinomial posterior distribution.

    See Also
    --------
    numpy.random.Generator.gamma
    numpy.random.Generator.dirichlet

    Notes
    -----
    This function uses a Gamma distribution to replace a Dirichlet distribution. The
    result is precisely identical to (and reproducible given the same seed):

    .. code-block:: python
       return clr(np.apply_along_axis(rng.dirichlet, axis=1, arr=matrix))

    A Dirichlet distribution is essentially a standard Gamma distribution normalized
    by row sums. Meanwhile, CLR is independent of scale, therefore the normalization
    step can be omitted.

    `gamma` can vectorize to a 2-D array whereas `dirichlet` cannot.

    """
    return clr(rng.gamma(shape=matrix, scale=1.0, size=matrix.shape), validate=False)


def dirmult_ttest(
    table,
    grouping,
    treatment=None,
    reference=None,
    pseudocount=0.5,
    draws=128,
    p_adjust="holm",
    seed=None,
):
    r"""*T*-test using Dirichlet-multinomial distribution.

    The Dirichlet-multinomial distribution is a compound distribution that
    combines a Dirichlet distribution over the probabilities of a multinomial
    distribution. This distribution is used to model the distribution of
    species abundances in a community.

    To perform the *t*-test, we first fit a Dirichlet-multinomial distribution
    for each sample, and then we compute the fold change and *p*-value for each
    feature. The fold change is computed as the difference between the
    samples of the two groups. *t*-tests are then performed on the posterior
    samples, drawn from each Dirichlet-multinomial distribution. The
    log-fold changes as well as their credible intervals, the *p*-values and
    the multiple comparison corrected *p*-values are reported.

    This process mirrors the approach performed by the R package "ALDEx2" [1]_.

    Additionally, this function excludes hits with a 95% confidence interval of
    fold-change crossing zero during any draw. This step further reduces false
    positive hits, especially among low-abundance features.

    .. versionchanged:: 0.7.0
        Computational efficiency significantly improved.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing count or proportional abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    grouping : pd.Series or 1-D array_like
        Vector indicating the assignment of samples to groups. These could be strings
        or integers denoting which group a sample belongs to. If it is a pandas Series
        and the table contains sample IDs, its index will be filtered and reordered to
        match the sample IDs. Otherwise, it must be the same length as the samples in
        the table.
    treatment : str, optional
        Name of the treatment group. The *t*-test is computed between the ``treatment``
        group and the ``reference`` group specified in the ``grouping`` vector. If
        omitted, the first group in the sorted order of all group names will be the
        treatment group.
    reference : str, optional
        Name of the reference group. See above. If omitted, all groups other than the
        treatment group will be combined as the reference group.

        .. versionchanged:: 0.7.0
            ``treatment`` and ``reference`` are now optional.

    pseudocount : float, optional
        A non-zero value added to the input counts to ensure that all of the
        estimated abundances are strictly greater than zero.
    draws : int, optional
        The number of draws from the Dirichilet-multinomial posterior distribution
        More draws provide higher uncertainty surrounding the estimated
        log-fold changes and *p*-values.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are Holm-
        Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-
        Hochberg ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported
        by statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance for drawing from the
        Dirichlet distribution. See :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    pd.DataFrame
        A table of features, their log-fold changes and other relevant statistics.

        - ``T-statistic``: *t*-statistic of Welch's *t*-test. The reported value is the
          average across all of the posterior draws.

        - ``Log2(FC)``: Expected log2-fold change of abundance from the reference group
          to the treatment group. The value is expressed in the center log ratio (see
          :func:`clr`) transformed coordinates. The reported value is the average of
          all of the log2-fold changes computed from each of the posterior draws.

        - ``CI(2.5)``: 2.5% quantile of the log2-fold change. The reported value is the
          minimum of all of the 2.5% quantiles computed from each of the posterior
          draws.

        - ``CI(97.5)``: 97.5% quantile of the log2-fold change. The reported value is
          the maximum of all of the 97.5% quantiles computed from each of the posterior
          draws.

        - ``pvalue``: *p*-value of Welch's *t*-test. The reported value is the average
          of all of the *p*-values computed from each of the posterior draws.

        - ``qvalue``: Corrected *p*-value of Welch's *t*-test for multiple comparisons.
          The reported value is the average of all of the *q*-values computed from each
          of the posterior draws.

        - ``Signif``: Whether feature is significantly differentially abundant between
          the treatment and reference groups. A feature marked as "True" suffice: 1)
          The *q*-value must be less than or equal to the significance level (0.05). 2)
          The confidence interval (CI(2.5)..CI(97.5)) must not overlap with zero.

        .. versionchanged:: 0.7.0
            ``df`` (degrees of freedom) was removed from the report, as this metric is
            inconsistent across draws.

        .. versionchanged:: 0.7.0
            Renamed ``T statistic`` as ``T-statistic``.
            Renamed ``Reject null hypothesis`` as ``Signif``.

    See Also
    --------
    dirmult_ttest
    scipy.stats.ttest_ind
    statsmodels.stats.weightstats.CompareMeans

    Notes
    -----
    The confidence intervals are computed using the mininum 2.5% and maximum
    97.5% bounds computed across all of the posterior draws.

    The reference frame here is the geometric mean. Extracting absolute log
    fold changes from this test assumes that the average feature abundance
    between the ``treatment`` and the ``reference`` groups are the same. If this
    assumption is violated, then the log-fold changes will be biased, and the
    *p*-values will not be reliable. However, the bias is the same across each
    feature, as a result the ordering of the log-fold changes can still be useful.

    One benefit of using the Dirichlet-multinomial distribution is that the
    statistical power increases with regards to the abundance magnitude. More counts
    per sample will shrink the size of the confidence intervals, and can result in
    lower *p*-values.

    References
    ----------
    .. [1] Fernandes, A. D., Reid, J. N., Macklaim, J. M., McMurrough, T. A., Edgell,
       D. R., & Gloor, G. B. (2014). Unifying the analysis of high-throughput
       sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and
       selective growth experiments by compositional data analysis. Microbiome, 2,
       1-13.

    Examples
    --------
    >>> import pandas as pd
    >>> from skbio.stats.composition import dirmult_ttest
    >>> table = pd.DataFrame([[20,  110, 100, 101, 100, 103, 104],
    ...                       [33,  110, 120, 100, 101, 100, 102],
    ...                       [12,  110, 100, 110, 100, 50,  90],
    ...                       [202, 201, 9,  10, 10, 11, 11],
    ...                       [200, 202, 10, 10, 13, 10, 10],
    ...                       [203, 201, 14, 10, 10, 13, 12]],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'],
    ...                      columns=['b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7'])
    >>> grouping = pd.Series(['treatment', 'treatment', 'treatment',
    ...                       'placebo', 'placebo', 'placebo'],
    ...                      index=['s1', 's2', 's3', 's4', 's5', 's6'])
    >>> result = dirmult_ttest(table, grouping, 'treatment', 'placebo', seed=0)
    >>> result
        T-statistic  Log2(FC)   CI(2.5)  CI(97.5)    pvalue    qvalue  Signif
    b1   -17.178600 -4.991987 -7.884498 -2.293463  0.003355  0.020131    True
    b2   -16.873187 -2.533729 -3.594590 -1.462339  0.001064  0.007446    True
    b3     6.942727  1.627677 -1.048219  4.750792  0.021130  0.068310   False
    b4     6.522786  1.707221 -0.467481  4.164998  0.013123  0.065613   False
    b5     6.654142  1.528243 -1.036910  3.978387  0.019360  0.068310   False
    b6     3.839520  1.182343 -0.702656  3.556061  0.045376  0.068310   False
    b7     7.600734  1.480232 -0.601277  4.043888  0.017077  0.068310   False

    """
    from statsmodels.stats.weightstats import CompareMeans

    rng = get_rng(seed)

    matrix, samples, features = _ingest_table(table)
    _check_composition(np, matrix)

    # handle zero values
    if pseudocount:
        matrix = matrix + pseudocount

    # get sample indices of treatment and reference groups
    groups, labels = _check_grouping(grouping, matrix, samples)
    trt_idx, ref_idx = _check_trt_ref_groups(treatment, reference, groups, labels)

    cm_params = dict(alternative="two-sided", usevar="unequal")

    # initiate results
    m = matrix.shape[1]
    delta = np.zeros(m)  # inter-group difference
    tstat = np.zeros(m)  # t-test statistic
    pval = np.zeros(m)  # t-test p-value
    lower = np.full(m, np.inf)  # 2.5% percentile of distribution
    upper = np.full(m, -np.inf)  # 97.5% percentile of distribution

    for i in range(draws):
        # Resample data in a Dirichlet-multinomial distribution.
        dir_mat = _dirmult_draw(matrix, rng)

        # Stratify data by group (treatment vs. reference).
        trt_mat = dir_mat[trt_idx]
        ref_mat = dir_mat[ref_idx]

        # Calculate the difference between the two means.
        delta += trt_mat.mean(axis=0) - ref_mat.mean(axis=0)

        # Create a CompareMeans object for statistical testing.
        # Welch's t-test is also available in SciPy's `ttest_ind` (with `equal_var=
        # False`). The current code uses statsmodels' `CompareMeans` instead because
        # it additionally returns confidence intervals.
        cm = CompareMeans.from_data(trt_mat, ref_mat)

        # Perform Welch's t-test to assess the significance of difference.
        tstat_, pval_, _ = cm.ttest_ind(value=0, **cm_params)
        tstat += tstat_
        pval += pval_

        # Calculate confidence intervals.
        # The final lower and upper bounds are the minimum and maximum of all lower
        # and upper bounds seen during sampling, respectively.
        lower_, upper_ = cm.tconfint_diff(alpha=0.05, **cm_params)
        np.minimum(lower, lower_, out=lower)
        np.maximum(upper, upper_, out=upper)

    # Normalize metrics to averages over all replicates.
    delta /= draws
    tstat /= draws
    pval /= draws

    # Correct p-values for multiple comparison.
    if p_adjust is not None:
        qval = _check_p_adjust(p_adjust)(pval)
    else:
        qval = pval
    reject = qval <= 0.05

    # Test if confidence interval includes 0.
    # A significant result (i.e., reject null hypothesis) must simultaneously suffice:
    # 1) q-value <= significance level. 2) confidence interval doesn't include 0.
    # This test is in addition to the original ALDEx2 method. It helps to reduce false
    # positive discoveries of low abundance.
    outer = ((lower > 0) & (upper > 0)) | ((lower < 0) & (upper < 0))
    reject &= outer

    # Convert all log fold changes to base 2.
    log2_ = np.log(2)
    delta /= log2_
    upper /= log2_
    lower /= log2_

    # construct report
    res = pd.DataFrame.from_dict(
        {
            "T-statistic": tstat,
            "Log2(FC)": delta,
            "CI(2.5)": lower,
            "CI(97.5)": upper,
            "pvalue": pval,
            "qvalue": qval,
            "Signif": reject,
        }
    )
    if features is not None:
        res.index = features

    return res


def _type_cast_to_float(df):
    """Attempt to cast all of the values in dataframe to float.

    This will try to type cast all of the series within the dataframe into floats. If a
    column cannot be type casted, it will be kept as is.

    Parameters
    ----------
    df : pd.DataFrame

    Returns
    -------
    pd.DataFrame

    """
    df = df.copy()
    for col in df.select_dtypes(exclude=["float64"]).columns:
        try:
            df[col] = df[col].astype("float64")
        except (ValueError, TypeError):
            continue
    return df


def dirmult_lme(
    table,
    metadata,
    formula,
    grouping,
    pseudocount=0.5,
    draws=128,
    p_adjust="holm",
    seed=None,
    re_formula=None,
    vc_formula=None,
    model_kwargs={},
    fit_method=None,
    fit_converge=False,
    fit_warnings=False,
    fit_kwargs={},
):
    r"""Fit a Dirichlet-multinomial linear mixed effects model.

    .. versionadded:: 0.7.0

    The Dirichlet-multinomial distribution is a compound distribution that
    combines a Dirichlet distribution over the probabilities of a multinomial
    distribution. This distribution is used to model the distribution of
    species abundances in a community.

    To fit the linear mixed effects model we first fit a Dirichlet-multinomial
    distribution for each sample, and then we compute the fold change and
    *p*-value for each feature. The fold change is computed as the slopes
    from the resulting model. Statistical tests are then performed on the posterior
    samples, drawn from each Dirichlet-multinomial distribution. The
    log-fold changes as well as their credible intervals, the *p*-values and
    the multiple comparison corrected *p*-values are reported.

    This function uses the :class:`~statsmodels.regression.mixed_linear_model.MixedLM`
    class from statsmodels.

    .. note::
        Because the analysis iteratively runs many numeric optimizations, it can take
        longer than usual to finish. Please allow extra time for completion.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing count or proportional abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    metadata : pd.DataFrame or 2-D array_like
        The metadata for the model. Rows correspond to samples and columns correspond
        to covariates in the model. Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    formula : str or generic Formula object
        The formula defining the model. Refer to `Patsy's documentation
        <https://patsy.readthedocs.io/en/latest/formulas.html>`_ on how to specify
        a formula.
    grouping : str, pd.Series or 1-D array_like
        A vector or a metadata column name indicating the assignment of samples to
        groups. Samples are independent between groups during model fitting.
    pseudocount : float, optional
        A non-zero value added to the input counts to ensure that all of the
        estimated abundances are strictly greater than zero. Default is 0.5.
    draws : int, optional
        Number of draws from the Dirichlet-multinomial posterior distribution.
        Default is 128.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are Holm-
        Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-
        Hochberg ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported
        by statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance for drawing from the
        Dirichlet distribution. See :func:`details <skbio.util.get_rng>`.
    re_formula : str, optional
        Random coefficient formula. See :meth:`MixedLM.from_formula
        <statsmodels.regression.mixed_linear_model.MixedLM.from_formula>` for details.
    vc_formula : str, optional
        Variance component formula. See ``MixedLM.from_formula`` for details.
    model_kwargs : dict, optional
        Additional keyword arguments to pass to ``MixedLM``.
    fit_method : str or list of str, optional
        Optimization method for model fitting. Can be a single method name, or a list
        of method names to be tried sequentially. See `statsmodels.optimization
        <https://www.statsmodels.org/stable/optimization.html>`_
        for available methods. If None, a default list of methods will be tried.
    fit_converge : bool, optional
        If True, model fittings that were completed but did not converge will be
        excluded from the calculation of final statistics. Default is False.
    fit_warnings : bool, optional
        Issue warnings if any during the model fitting process. Default is False.
        Warnings are usually issued when the optimization methods do not converge,
        which is common in the analysis. Default is False.
    fit_kwargs : dict, optional
        Additional keyword arguments to pass to :meth:`MixedLM.fit
        <statsmodels.regression.mixed_linear_model.MixedLM.fit>`.

    Returns
    -------
    pd.DataFrame
        A table of features and covariates, their log-fold changes and other relevant
        statistics.

        - ``FeatureID``: Feature identifier, i.e., dependent variable.

        - ``Covariate``: Covariate name, i.e., independent variable.

        - ``Reps``: Number of Dirichlet-multinomial posterior draws that supported the
          reported statistics, i.e., the number of successful model fittings on this
          feature. Max: ``draws`` (if none failed). Min: 0 (in which case all
          statistics are NaN).

        - ``Log2(FC)``: Expected log2-fold change of abundance from the reference
          category to the covariate category defined in the formula. The value is
          expressed in the center log ratio (see :func:`clr`) transformed coordinates.
          The reported value is the average of all of the log2-fold changes computed
          from each of the posterior draws.

        - ``CI(2.5)``: 2.5% quantile of the log2-fold change. The reported value is the
          minimum of all of the 2.5% quantiles computed from each of the posterior
          draws.

        - ``CI(97.5)``: 97.5% quantile of the log2-fold change. The reported value is
          the maximum of all of the 97.5% quantiles computed from each of the posterior
          draws.

        - ``pvalue``: *p*-value of the linear mixed effects model. The reported value
          is the average of all of the *p*-values computed from each of the posterior
          draws.

        - ``qvalue``: Corrected *p*-value of the linear mixed effects model for multiple
          comparisons. The reported value is the average of all of the *q*-values
          computed from each of the posterior draws.

        - ``Signif``: Whether the covariate category is significantly differentially
          abundant from the reference category. A feature-covariate pair marked as
          "True" suffice: 1) The *q*-value must be less than or equal to the
          significance level (0.05). 2) The confidence interval (CI(2.5)..CI(97.5))
          must not overlap with zero.

    See Also
    --------
    dirmult_ttest
    statsmodels.formula.api.mixedlm
    statsmodels.regression.mixed_linear_model.MixedLM

    Examples
    --------
    >>> import pandas as pd
    >>> from skbio.stats.composition import dirmult_lme
    >>> table = pd.DataFrame(
    ...     [[1.00000053, 6.09924644],
    ...      [0.99999843, 7.0000045],
    ...      [1.09999884, 8.08474053],
    ...      [1.09999758, 1.10000349],
    ...      [0.99999902, 2.00000027],
    ...      [1.09999862, 2.99998318],
    ...      [1.00000084, 2.10001257],
    ...      [0.9999991, 3.09998418],
    ...      [0.99999899, 3.9999742],
    ...      [1.10000124, 5.0001796],
    ...      [1.00000053, 6.09924644],
    ...      [1.10000173, 6.99693644]],
    ...     index=['u1', 'u2', 'u3', 'x1', 'x2', 'x3', 'y1', 'y2', 'y3', 'z1',
    ...            'z2', 'z3'],
    ...     columns=['Y1', 'Y2'])
    >>> metadata = pd.DataFrame(
    ...     {'patient': [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
    ...      'treatment': [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
    ...      'time': [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3],},
    ...     index=['x1', 'x2', 'x3', 'y1', 'y2', 'y3', 'z1', 'z2', 'z3', 'u1',
    ...            'u2', 'u3'])
    >>> result = dirmult_lme(table, metadata, formula='time + treatment',
    ...                      grouping='patient', seed=0, p_adjust='sidak')
    >>> result
      FeatureID  Covariate  Reps  Log2(FC)   CI(2.5)  CI(97.5)    pvalue  \
    0        Y1       time   128 -0.210769 -1.532255  1.122148  0.403737
    1        Y1  treatment   128 -0.744061 -3.401978  1.581917  0.252057
    2        Y2       time   128  0.210769 -1.122148  1.532255  0.403737
    3        Y2  treatment   128  0.744061 -1.581917  3.401978  0.252057
    <BLANKLINE>
         qvalue  Signif
    0  0.644470   False
    1  0.440581   False
    2  0.644470   False
    3  0.440581   False

    """
    from patsy import dmatrix
    from scipy.optimize import OptimizeWarning
    from statsmodels.regression.mixed_linear_model import MixedLM, VCSpec
    from statsmodels.tools.sm_exceptions import ConvergenceWarning

    rng = get_rng(seed)

    matrix, samples, features = _ingest_table(table)

    _check_composition(np, matrix)

    n_feats = matrix.shape[1]
    if n_feats < 2:
        raise ValueError("Table must have at least two features.")
    if features is None:
        features = np.arange(n_feats)

    # validate metadata
    metadata = _check_metadata(metadata, matrix, samples)

    # cast metadata to numbers where applicable
    metadata = _type_cast_to_float(metadata)

    # Instead of directly calling `MixedLM.from_formula` on merged table + metadata,
    # the following code converts metadata into a design matrix based on the formula
    # (as well as re_formula and vc_formula, if applicable), and calls `MixedLM`.
    # This is because the design matrix (independent variable) is always the same
    # whereas the table (dependent variable) is resampled in every replicate. Fixing
    # the design matrix can save conversion overheads.

    # Create a design matrix based on metadata and formula.
    dmat = dmatrix(formula, metadata, return_type="matrix")

    # Obtain the list of covariates by selecting the relevant columns
    covars = dmat.design_info.column_names

    # Remove intercept since it is not a covariate, and is included by default.
    # Then determine the range of rows to be extracted from the model fitting result.
    # (See also `result.model.k_fe`, number of fixed effects.)
    if covars[0] == "Intercept":
        covars = covars[1:]
        n_covars = len(covars)
        covar_range = slice(1, n_covars + 1)
    else:
        n_covars = len(covars)
        covar_range = slice(0, n_covars)

    exog_mat = np.asarray(dmat)

    # parse grouping
    if isinstance(grouping, str):
        try:
            grouping = metadata[grouping].to_numpy()
        except KeyError:
            raise ValueError("Grouping is not a column in the metadata.")
        uniq, grouping = np.unique(grouping, return_inverse=True)
    else:
        uniq, grouping = _check_grouping(grouping, matrix, samples)
    n_groups = len(uniq)

    # random effects matrix
    if re_formula is not None:
        exog_re = np.asarray(dmatrix(re_formula, metadata, return_type="matrix"))
    else:
        exog_re = None

    # variance component matrices
    # see: https://www.statsmodels.org/v0.12.2/examples/notebooks/generated/
    # variance_components.html
    if vc_formula is not None:
        metas = [metadata.iloc[grouping == x] for x in range(n_groups)]
        names, cols, mats = [], [], []
        for name, formula in vc_formula.items():
            names.append(name)
            dmats = [dmatrix(formula, x, return_type="matrix") for x in metas]
            cols.append([x.design_info.column_names for x in dmats])
            mats.append([np.asarray(x) for x in dmats])
        exog_vc = VCSpec(names, cols, mats)
    else:
        exog_vc = None

    # handle zero values
    if pseudocount:
        matrix = matrix + pseudocount

    # initiate results
    shape = (n_feats, n_covars)
    coef = np.zeros(shape)  # coefficient (fold change)
    pval = np.zeros(shape)  # p-value
    lower = np.full(shape, np.inf)  # 2.5% CI
    upper = np.full(shape, -np.inf)  # 97.5% CI

    # number of replicates (draws) LME fitting is successful for each feature
    fitted = np.zeros(n_feats, dtype=int)

    fit_fail_msg = "LME fit failed for feature {} in replicate {}, outputting NaNs."

    with catch_warnings():
        if not fit_warnings:
            simplefilter("ignore", UserWarning)
            simplefilter("ignore", ConvergenceWarning)
            simplefilter("ignore", OptimizeWarning)
            # This is temporary because statsmodels calls scipy in a deprecated way as
            # of v0.14.4.
            simplefilter("ignore", DeprecationWarning)

        for i in range(draws):
            # Resample data in a Dirichlet-multinomial distribution.
            dir_mat = _dirmult_draw(matrix, rng)

            # Fit a linear mixed effects (LME) model for each feature.
            for j in range(n_feats):
                model = MixedLM(
                    dir_mat[:, j],
                    exog_mat,
                    grouping,
                    exog_re=exog_re,
                    exog_vc=exog_vc,
                    **model_kwargs,
                )

                # model fitting (computationally expensive)
                try:
                    result = model.fit(method=fit_method, **fit_kwargs)

                # There are many ways model fitting may fail. Examples are LinAlgError,
                # RuntimeError, OverflowError, and ZeroDivisionError. If any error
                # occurs, the function will still proceed but the current run will be
                # discarded.
                except Exception:
                    warn(fit_fail_msg.format(features[j], i), UserWarning)
                    continue

                # It is common that model fitting successfully finished (no error) but
                # the optimizer did not converge, making the calculated statistics less
                # reliable. The `fit_converge` flag can discard these runs.
                if fit_converge and not result.converged:
                    warn(fit_fail_msg.format(features[j], i), UserWarning)
                    continue

                # update results
                coef[j] += result.params[covar_range]
                pval[j] += result.pvalues[covar_range]

                # calculate confidence interval and update results
                ci = result.conf_int()
                np.minimum(lower[j], ci[covar_range, 0], out=lower[j])
                np.maximum(upper[j], ci[covar_range, 1], out=upper[j])

                fitted[j] += 1

    # deal with fitting failures
    all_fail_msg = "LME fit failed for {} features in all replicates, reporting NaNs."
    mask = fitted > 0
    n_failed = n_feats - mask.sum()
    # all succeeded
    if n_failed == 0:
        mask = slice(None)
    # some failed
    elif n_failed < n_feats:
        warn(all_fail_msg.format(n_failed), UserWarning)
        for x in (coef, pval, lower, upper):
            x[~mask] = np.nan
    # all failed
    else:
        raise ValueError("LME fit failed for all features in all replicates.")

    # normalize metrics to averages over all successful replicates
    fitted_ = fitted.reshape(-1, 1)[mask]
    for x in (coef, pval):
        x[mask] /= fitted_

    # convert all log fold changes to base 2
    log2_ = np.log(2)
    for x in (coef, lower, upper):
        x[mask] /= log2_

    # correct p-values for multiple comparison
    # (only valid replicates are included)
    if p_adjust is not None:
        func = _check_p_adjust(p_adjust)
        qval = np.full(shape, np.nan)
        qval[mask] = np.apply_along_axis(func, 1, pval[mask])
    else:
        qval = pval

    # get significant results (q-value <= 0.05 and CI doesn't cross 0)
    # see `dirmult_ttest`
    reject = np.full(shape, np.nan)
    ii = np.where(mask)[0] if n_failed else np.arange(n_feats)
    for i in ii:
        outer = ((lower[i] > 0) & (upper[i] > 0)) | ((lower[i] < 0) & (upper[i] < 0))
        reject[i] = (qval[i] <= 0.05) & outer

    # construct report
    res = pd.DataFrame.from_dict(
        {
            "FeatureID": [x for x in features for _ in range(n_covars)],
            "Covariate": list(covars) * n_feats,
            "Reps": np.repeat(fitted, n_covars),
            "Log2(FC)": coef.ravel(),
            "CI(2.5)": lower.ravel(),
            "CI(97.5)": upper.ravel(),
            "pvalue": pval.ravel(),
            "qvalue": qval.ravel(),
        }
    )
    # pandas' nullable boolean type
    res["Signif"] = pd.Series(reject.ravel(), dtype="boolean")
    return res


register_aliases(modules[__name__])
