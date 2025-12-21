# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from skbio.util import get_rng
from skbio.stats.composition import closure
from skbio.table._tabular import _ingest_table

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Sequence
    from numpy.typing import ArrayLike, NDArray
    from skbio.tree import TreeNode
    from skbio.util._typing import TableLike, SeedLike


def _validate_labels(  # type: ignore[return]
    labels: NDArray, n: int
) -> tuple[NDArray, NDArray]:
    r"""Ensure that the provided label is appropriate for the provided matrix.

    Parameters
    ----------
    labels : array_like of shape (n_samples,) or (n_samples, n_classes)
        Class labels for the data.
    n : int
        Number of samples.

    Returns
    -------
    labels_index : ndarray of shape (n_samples,)
        The class labels in 1-D format. Contains indices from 0 to n_classes - 1. If
        the input was already 1-D, this is the same as the input. If the input was
        one-hot encoded, this is the argmax reconstruction.
    labels_one_hot : ndarray of shape (n_samples, n_classes)
        The class labels in one-hot encoded format. Each row contains exactly one 1 and
        the rest 0s, indicating the class membership for that sample.

    Raises
    ------
    ValueError
        If labels are not 1-D or 2-D.
    ValueError
        If labels' sample count is not ``n``.
    ValueError
        If index labels contain non-integer values.
    ValueError
        If index labels do not start from zero.
    ValueError
        If index labels are not consecutive integers.
    ValueError
        If one-hot encoded labels are not valid.

    """
    labels = np.asarray(labels)

    if labels.ndim not in (1, 2):
        raise ValueError(
            f"Labels should be 1-D or 2-D, but got {labels.ndim} dimensions instead."
        )

    if labels.shape[0] != n:
        raise ValueError(
            f"Number of labels ({labels.shape[0]}) does not match number of samples "
            f"in the data ({n})."
        )

    # input as indexed labels
    if labels.ndim == 1:
        # if labels is not integer-type, make sure it contains only whole numbers,
        # then convert to int
        if not np.issubdtype(labels.dtype, np.integer):
            if not (np.mod(labels, 1) == 0).all():
                raise ValueError(f"Labels must only contain integer values.")
            labels = labels.astype(int)

        # check that labels are 0-indexed
        unique_labels = np.unique(labels)
        if unique_labels.min() != 0:
            raise ValueError("Labels must be zero-indexed. Minimum value must be 0.")
        n_classes = len(unique_labels)
        exp_labels = np.arange(n_classes)
        # check that label is consecutive integers starting at 0
        if not np.array_equal(unique_labels, exp_labels):
            raise ValueError(
                "Labels must be consecutive integers from 0 to n_classes - 1."
            )
        # perform one-hot encoding
        labels_one_hot = np.eye(n_classes, dtype=int)[labels]
        return labels, labels_one_hot

    # input as one-hot encoded labels
    else:
        # all rows should sum to 1
        if not np.allclose(labels.sum(axis=1), 1):
            raise ValueError("Labels are not properly one-hot encoded.")
        # generated label indices
        labels_index = labels.argmax(axis=1)
        return labels_index, labels


def _normalize_matrix(matrix: NDArray) -> NDArray:
    r"""Normalize a data matrix if needed such that each row sums to 1.

    Parameters
    ----------
    matrix : ndarray of shape (n_samples, n_features)
        Original matrix.

    Returns
    -------
    ndarray of shape (n_samples, n_features)
        Normalized matrix.

    """
    if np.allclose(matrix.sum(axis=1), 1):
        return matrix
    return closure(matrix)


def _all_pairs(n: int) -> NDArray:
    r"""Get all pairs of sample indices.

    Parameters
    ----------
    n : int
        Number of samples.

    Returns
    -------
    ndarray of shape (n_pairs, 2), dtype=int
        Pairs of sample indices.

    """
    return np.column_stack(np.triu_indices(n, k=1))


def _intra_class_pairs(labels: ArrayLike) -> NDArray:
    r"""Get pairs of sample indices within each class.

    Parameters
    ----------
    labels : array_like of shape (n_samples,), dtype=int
        Class labels of samples.

    Returns
    -------
    ndarray of shape (n_pairs, 2), dtype=int
        Pairs of sample indices.

    """
    pairs = []
    for label in np.unique(labels):
        idx = np.where(labels == label)[0]
        if idx.size > 1:
            i, j = np.triu_indices(idx.size, k=1)
            pairs.append(np.column_stack((idx[i], idx[j])))
    if pairs:
        return np.vstack(pairs)
    else:
        return np.empty((0, 2), dtype=np.intp)


def _format_input(
    table: TableLike,
    labels: Sequence | NDArray,
    intra_class: bool | None = False,
    normalize: bool | None = False,
    taxa: Sequence | NDArray | None = None,
) -> tuple[
    NDArray,
    NDArray | None,
    NDArray,
    Sequence | NDArray | None,
]:
    """Format input data for augmentation.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table.
    labels : array_like of shape (n_samples,) or (n_samples, n_classes), optional
        Class labels for the data.
    intra_class : bool, optional
        If True, pair samples within each class.
    normalize : bool, optional
        If True, normalize the data.
    taxa : array_like of shape (n_features,), optional
        Feature IDs corresponding to tip names in a tree.

    Returns
    -------
    matrix : ndarray of shape (n_samples, n_features)
        Data matrix (X).
    labels : ndarray of shape (n_samples, n_classes), optional
        Class labels (one-hot encoded) (Y).
    pairs : ndarray of shape (n_pairs, 2)
        Pairs of sample indices to be mixed up.
    taxa : array_like of shape (n_features,), optional
        Feature IDs corresponding to tip names in a tree.

    """
    matrix, _, taxa = _ingest_table(table, feature_ids=taxa)
    if normalize:
        matrix = _normalize_matrix(matrix)
    n_samples = matrix.shape[0]
    if labels is None:
        pairs = _all_pairs(n_samples)
        labels_2d = None
    else:
        labels_1d, labels_2d = _validate_labels(labels, n_samples)
        if intra_class:
            pairs = _intra_class_pairs(labels_1d)
        else:
            pairs = _all_pairs(n_samples)
    if len(pairs) == 0:
        raise ValueError("Cannot find a pair of samples to mix.")
    return matrix, labels_2d, pairs, taxa


def _make_aug_labels(
    labels: NDArray | None,
    pairs: NDArray,
    lams: NDArray,
    intra_class: bool,
) -> NDArray | None:
    r"""Generate class labels for synthetic samples.

    Parameters
    ----------
    labels : ndarray of shape (n, n_classes), or None
        Original class labels in one-hot encoded format.
    pairs : ndarray of shape (n_pairs, 2)
        Pairs of sample indices that were mixed up to create synthetic samples.
    lams : ndarray of shape (n, 1)
        Lambda parameters (probabilities) used during mixing-up.
    intra_class : bool
        If True, samples were paired within each class.

    Returns
    -------
    aug_labels : ndarray of shape (n, n_classes), or None
        Synthetic class labels in one-hot encoded format.

    """
    if labels is None:
        return None
    elif intra_class:
        return labels[pairs[:, 0]]
    else:
        return lams * labels[pairs[:, 0]] + (1.0 - lams) * labels[pairs[:, 1]]


def _format_output(
    aug_matrix: NDArray,
    aug_labels: NDArray | None,
    matrix: NDArray,
    labels: NDArray | None,
    append: bool = False,
) -> tuple[NDArray, NDArray]:
    """Format output data of augmentation.

    Parameters
    ----------
    aug_matrix : array_like of shape (n, n_features)
        Augmented data matrix.
    aug_labels : array_like of shape (n, n_classes), optional
        Augmented class labels in one-hot encoded format.
    matrix : ndarray of shape (n_samples, n_features), optional
        Original data matrix.
    labels : ndarray of shape (n_samples, n_classes), optional
        Original class labels in one-hot encoded format.

    Returns
    -------
    aug_matrix : ndarray of shape ([n_samples + ]n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape ([n_samples + ]n, n_classes), optional
        Augmented class labels in one-hot encoded format.

    """
    # concatenate original and synthetic data
    if append:
        aug_matrix = np.concatenate([matrix, aug_matrix])
        if aug_labels is not None:
            aug_labels = np.concatenate([labels, aug_labels])

    # not necessary because upstream code already enforced it; but for safety
    else:
        aug_matrix = np.asarray(aug_matrix)
        if aug_labels is not None:
            aug_labels = np.asarray(aug_labels)

    return aug_matrix, aug_labels


def mixup(
    table: TableLike,
    n: int,
    labels: NDArray | None = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    append: bool = False,
    seed: SeedLike | None = None,
) -> tuple[NDArray, NDArray | None]:
    r"""Data augmentation by vanilla mixup.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table to be augmented. See
        :ref:`supported formats <table_like>`.
    n : int
        Number of synthetic samples to generate.
    labels : array_like of shape (n_samples,) or (n_samples, n_classes), optional
        Class labels for the data. Accepts either indices (1-D) or one-hot encoded
        labels (2-D).
    intra_class : bool, optional
        If True, synthetic samples will be created by mixing samples within each class.
        If False (default), any samples regardless of class can be mixed.
    alpha : float, optional
        Shape parameter of the beta distribution.
    append : bool, optional
        If True, the returned data include both the original and synthetic samples. If
        False (default), only the synthetic samples are returned.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    See Also
    --------
    phylomix
    aitchison_mixup
    compos_cutmix

    Returns
    -------
    aug_matrix : ndarray of shape (n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape (n, n_classes), optional
        Augmented class labels in one-hot encoded format. Available if ``labels`` are
        provided. One can call ``aug_labels.argmax(axis=1)`` to get class indices.

    Notes
    -----
    The mixup method is based on [1]_. It randomly selects two samples :math:`s_1` and
    :math:`s_2` from the data table, and generates a new sample :math:`s` by a linear
    combination of them, as follows:

    .. math::

       s = \lambda \cdot s_1 + (1 - \lambda) \cdot s_2

    where :math:`\lambda` is a mixing coefficient drawn from a beta distribution:

    .. math::

       \lambda \sim \mathrm{Beta}(\alpha, \alpha)

    The label :math:`y` is computed as the linear combination of the labels of the two
    samples (:math:`y_1` and :math:`y_2`):

    .. math::

       y = \lambda \cdot y_1 + (1 - \lambda) \cdot y_2

    This function shares the same core concept as PyTorch's
    `MixUp <https://pytorch.org/vision/
    main/generated/torchvision.transforms.v2.MixUp.html>`_ class. There are some key
    differences:

    1. This implementation returns synthetic samples and class labels from a dataset,
       while PyTorch's MixUp is applied on-the-fly during training to batches of data.

    2. This implementation randomly selects pairs of samples from the entire dataset,
       while PyTorch's implementation typically mixes consecutive samples in a batch
       (requiring prior shuffling).

    3. This implementation is tailored for biological omic data, while PyTorch's is
       primarily for image data.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.table import mixup
    >>> matrix = np.arange(40).reshape(4, 10)
    >>> labels = np.array([0, 1, 0, 1])
    >>> aug_matrix, aug_labels = mixup(matrix, n=5, labels=labels)
    >>> print(aug_matrix.shape)
    (5, 10)
    >>> print(aug_labels.shape)
    (5, 2)

    References
    ----------
    .. [1] Zhang, H., Cisse, M., Dauphin, Y. N., & Lopez-Paz, D. (2017). mixup: Beyond
       Empirical Risk Minimization. arXiv preprint arXiv:1710.09412.

    """
    rng = get_rng(seed)
    matrix, labels, pairs, _ = _format_input(table, labels, intra_class)

    # Draw pairs of samples to mix.
    pairs_sel = pairs[rng.integers(0, len(pairs), size=n)]

    # Draw mixing coefficients (lambda) from a beta distribution.
    lams = rng.beta(alpha, alpha, size=(n, 1))

    # Collect paired sample data.
    x1 = matrix[pairs_sel[:, 0]]
    x2 = matrix[pairs_sel[:, 1]]

    # Mix paired sample data.
    aug_matrix = lams * x1 + (1.0 - lams) * x2

    # Generate synthetic class labels.
    aug_labels = _make_aug_labels(labels, pairs_sel, lams, intra_class)

    return _format_output(aug_matrix, aug_labels, matrix, labels, append)


def _aitchison_add(x: NDArray, y: NDArray) -> NDArray:
    r"""Perform Aitchison addition on two samples or sample groups.

    Parameters
    ----------
    x : ndarray of shape (..., n_features)
        The first sample(s).
    y : ndarray of shape (..., n_features)
        The second sample(s).

    Returns
    -------
    ndarray of shape (..., n_features)
        The result of Aitchison addition.

    """
    return (prod := x * y) / prod.sum(axis=-1, keepdims=True)


def _aitchison_multiply(x: NDArray, p: float) -> NDArray:
    r"""Perform Aitchison multiplication on a sample (group) with a scalar.

    Parameters
    ----------
    x : ndarray of shape (..., n_features)
        The sample to multiply.
    p : float or ndarray of shape (...,)
        The scalar to multiply the sample by.

    Returns
    -------
    ndarray of shape (..., n_features)
        The result of Aitchison multiplication.

    """
    return (exp := x**p) / exp.sum(axis=-1, keepdims=True)


def aitchison_mixup(
    table: TableLike,
    n: int,
    labels: NDArray | None = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    normalize: bool = True,
    append: bool = False,
    seed: SeedLike | None = None,
) -> tuple[NDArray, NDArray | None]:
    r"""Data augmentation by Aitchison mixup.

    This function requires the data to be compositional (values per sample sum to one).
    If not, the function will automatically normalize them prior to augmentation.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table to be augmented. See
        :ref:`supported formats <table_like>`.
    n : int
        Number of synthetic samples to generate.
    labels : array_like of shape (n_samples,) or (n_samples, n_classes), optional
        Class labels for the data. Accepts either indices (1-D) or one-hot encoded
        labels (2-D).
    intra_class : bool, optional
        If ``True``, synthetic samples will be created by mixing samples within each
        class. If ``False`` (Default), any samples regardless of class can be mixed.
    alpha : float, optional
        Shape parameter of the beta distribution.
    normalize : bool, optional
        If True (default), and the input is not already compositional, scikit-bio's
        :func:`~skbio.stats.composition.closure` function will be called, ensuring
        values for each sample add up to 1.
    append : bool, optional
        If True, the returned data include both the original and synthetic samples. If
        False (default), only the synthetic samples are returned.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    aug_matrix : ndarray of shape (n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape (n, n_classes), optional
        Augmented class labels in one-hot encoded format. Available if ``labels`` are
        provided. One can call ``aug_labels.argmax(axis=1)`` to get class indices.

    See Also
    --------
    mixup
    compos_cutmix

    Notes
    -----
    The algorithm is based on [1]_, and leverages the Aitchison geometry to guide the
    augmentation of compositional data. It is essentially the vanilla mixup method in
    the Aitchison space.

    This method only works on compositional data, where a set of data points live in
    the simplex: :math:`x_i > 0`, and :math:`\sum_{i=1}^{p} x_i = 1`.

    An augmented sample :math:`s` is computed as the linear combination of two samples
    :math:`s_1` and :math:`s_2` in the Aitchison space:

    .. math::

       s = (\lambda \otimes  s_1) \oplus ((1 - \lambda) \otimes s_2)

    where :math:`\otimes` is the Aitchison scalar multiplication, defined as:

    .. math::

       \lambda \otimes s =
       \frac{1}{\sum_{i=1}^{n} s_i^{\lambda}}
       (s_1^{\lambda}, s_2^{\lambda}, ..., s_n^{\lambda})

    :math:`\oplus` is the Aitchison addition, defined as:

    .. math::

       s \oplus t =
       \frac{1}{\sum_{i=1}^{n} s_i t_i}
       (s_1 t_1, s_2 t_2, ..., s_n t_n)

    :math:`\lambda` is a mixing coefficient drawn from a beta distribution:

    .. math::

       \lambda \sim \mathrm{Beta}(\alpha, \alpha)

    The label :math:`y` is computed as the linear combination of the labels of the two
    samples (:math:`y_1` and :math:`y_2`):

    .. math::

       y = \lambda \cdot y_1 + (1 - \lambda) \cdot y_2

    By mixing the counts of two samples, Aitchison mixup preserves the compositional
    nature of the data, and the sum-to-one property.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.table import aitchison_mixup
    >>> matrix = np.arange(40).reshape(4, 10)
    >>> labels = np.array([0, 1, 0, 1])
    >>> aug_matrix, aug_labels = aitchison_mixup(matrix, n=5, labels=labels)
    >>> print(aug_matrix.shape)
    (5, 10)
    >>> print(aug_labels.shape)
    (5, 2)

    References
    ----------
    .. [1] Gordon-Rodriguez, E., Quinn, T., & Cunningham, J. P. (2022). Data
       augmentation for compositional data: Advancing predictive models of the
       microbiome. Advances in Neural Information Processing Systems, 35,
       20551-20565.

    """
    rng = get_rng(seed)
    matrix, labels, pairs, _ = _format_input(table, labels, intra_class, normalize)
    pairs_sel = pairs[rng.integers(0, len(pairs), size=n)]

    # Draw mixing coefficients (lambda) from a beta distribution.
    lams = rng.beta(alpha, alpha, size=(n, 1))

    # Collect paired sample data to mix.
    x1 = matrix[pairs_sel[:, 0]]
    x2 = matrix[pairs_sel[:, 1]]

    # Perform Aitchison scalar multiplication on each side of mixture.
    s1 = _aitchison_multiply(x1, lams)
    s2 = _aitchison_multiply(x2, 1.0 - lams)

    # Perform Aitchison addition to mix the two sides.
    aug_matrix = _aitchison_add(s1, s2)

    aug_labels = _make_aug_labels(labels, pairs_sel, lams, intra_class)
    return _format_output(aug_matrix, aug_labels, matrix, labels, append)


def compos_cutmix(
    table: TableLike,
    n: int,
    labels: NDArray | None = None,
    normalize: bool = True,
    append: bool = False,
    seed: SeedLike | None = None,
) -> tuple[NDArray, NDArray | None]:
    r"""Data augmentation by compositional cutmix.

    This function requires the data to be compositional (values per sample sum to one).
    If not, the function will automatically normalize them prior to augmentation.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table to be augmented. See
        :ref:`supported formats <table_like>`.
    n : int
        Number of synthetic samples to generate.
    labels : array_like of shape (n_samples,) or (n_samples, n_classes), optional
        Class labels for the data. Accepts either indices (1-D) or one-hot encoded
        labels (2-D).
    normalize : bool, optional
        If True (default), and the input is not already compositional, scikit-bio's
        :func:`~skbio.stats.composition.closure` function will be called, ensuring
        values for each sample add up to 1.
    append : bool, optional
        If True, the returned data include both the original and synthetic samples. If
        False (default), only the synthetic samples are returned.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

        .. note::
           This function does not have the ``intra_class`` parameter, as it always
           operates in intra-class mode in order to preserve the compositional
           structure within classes.

    See Also
    --------
    mixup
    aitchison_mixup

    Returns
    -------
    aug_matrix : ndarray of shape (n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape (n, n_classes), optional
        Augmented class labels in one-hot encoded format. Available if ``labels`` are
        provided. One can call ``aug_labels.argmax(axis=1)`` to get class indices.

    Notes
    -----
    The compositional cutmix method was described in [1]_.

    This method randomly selects values from one of a pair of samples to generate
    a new sample. It has four steps:

    1. Draw a mixing coefficient :math:`\lambda` from a uniform distribution:

    .. math::

       \lambda \sim U(0, 1)

    2. Draw a binary selector :math:`I` for each feature from a Bernoulli distribution:

    .. math::

       I \sim \mathrm{Bernoulli}(\lambda)

    3. For the :math:`i`-th feature, set the augmented value :math:`x_i` as from sample
       1 if :math:`I_i = 0` or from sample 2 if :math:`I_i = 1`.

    4. Normalize the augment sample such that it is compositional (sum-to-one).

    .. math::

       s = \frac{1}{\sum_{i=1}^{n} x_i} (x_1, x_2, ..., x_n)

    This method is applied separately to samples of each class. If ``labels`` is None,
    all samples will be considered as the same class, and ``aug_labels`` will be
    returned as None.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.table import compos_cutmix
    >>> matrix = np.arange(40).reshape(4, 10)
    >>> labels = np.array([0, 1, 0, 1])
    >>> aug_matrix, aug_labels = compos_cutmix(matrix, n=5, labels=labels)
    >>> print(aug_matrix.shape)
    (5, 10)
    >>> print(aug_labels.shape)
    (5, 2)

    References
    ----------
    .. [1] Gordon-Rodriguez, E., Quinn, T., & Cunningham, J. P. (2022). Data
       augmentation for compositional data: Advancing predictive models of the
       microbiome. Advances in Neural Information Processing Systems, 35,
       20551-20565.

    """
    rng = get_rng(seed)
    matrix, labels, pairs, _ = _format_input(table, labels, True, normalize)
    pairs_sel = pairs[rng.integers(0, len(pairs), size=n)]

    # Draw binary gates to control the choice between samples 1 and 2 per feature.
    probs = rng.uniform(0, 1, size=(n, 1))
    gates = rng.binomial(1, probs, size=(n, matrix.shape[1])).astype(bool)

    # Mix sample pairs based on binary gates.
    samples_sel = matrix[pairs_sel]
    aug_matrix = np.where(gates, samples_sel[:, 0], samples_sel[:, 1])

    # Normalize synthetic samples (although original samples are compositional,
    # synthetic one at this point are not, so they need to be normalized)
    aug_matrix = closure(aug_matrix)

    aug_labels = labels[pairs_sel[:, 0]] if labels is not None else None
    return _format_output(aug_matrix, aug_labels, matrix, labels, append)


def _indices_under_nodes(
    tree: TreeNode, taxa: Sequence | NDArray
) -> list[dict[int, None]]:
    """Pre-compute feature indices descending from each node.

    Parameters
    ----------
    tree : TreeNode
        Reference tree.
    taxa : sequence of str
        Taxa (tip names) corresponding to feature indices.

    Returns
    -------
    list of dict of {int: None}
        Feature indices descending from each node.

    Raises
    ------
    ValueError
        If there are duplicate taxa.
    ValueError
        If some taxa are not present as tip names in the tree.

    See Also
    --------
    skbio.tree._utils._validate_taxa_and_tree

    Notes
    -----
    Only internal nodes with at least one descending feature are included.

    This function does not ensure all tip names in the tree are unique.

    The output is a list of dicts instead of just list of lists. The reason is
    explained in the inline comments of `phylomix`.

    """
    n_taxa = len(taxa)
    taxon_map = {taxon: i for i, taxon in enumerate(taxa)}
    if len(taxon_map) < n_taxa:
        raise ValueError("All taxa must be unique.")

    seen: set[str] = set()
    seen_add = seen.add
    res: list[dict[int, None]] = []
    res_append = res.append
    dict_fromkeys = dict.fromkeys

    for node in tree.postorder(include_self=True):
        if node.children:
            lst = []
            for child in node.children:
                lst.extend(child._taxa)
                del child._taxa
            node._taxa = lst
            if lst:
                res_append(dict_fromkeys(lst))
        else:
            if (taxon := node.name) in taxon_map:
                node._taxa = [taxon_map[taxon]]
                seen_add(taxon)
            else:
                node._taxa = []
    del tree._taxa

    if missing := n_taxa - len(seen):
        raise ValueError(f"{missing} taxa are not present as tip names in the tree.")

    return res


def phylomix(
    table: TableLike,
    n: int,
    tree: TreeNode,
    taxa: ArrayLike | None = None,
    labels: NDArray | None = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    append: bool = False,
    seed: SeedLike | None = None,
) -> tuple[NDArray, NDArray | None]:
    r"""Data augmentation by PhyloMix.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table to be augmented. See
        :ref:`supported formats <table_like>`.
    n : int
        Number of synthetic samples to generate.
    tree : :class:`~skbio.tree.TreeNode`
        Tree structure modeling the relationships between features.
    taxa : array_like of shape (n_features,), optional
        Taxa (tip names) in ``tree`` corresponding to individual features. Can be
        omitted if ``table`` already contains feature IDs that are taxa. Otherwise
        they need to be explicitly provided. Should be a subset of taxa in the tree.
    labels : array_like of shape (n_samples,) or (n_samples, n_classes), optional
        Class labels for the data. Accepts either indices (1-D) or one-hot encoded
        labels (2-D).
    intra_class : bool, optional
        If ``True``, synthetic samples will be created by mixing samples within each
        class. If ``False`` (Default), any samples regardless of class can be mixed.
    alpha : float, optional
        Shape parameter of the beta distribution.
    append : bool, optional
        If True, the returned data include both the original and synthetic samples. If
        False (default), only the synthetic samples are returned.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    aug_matrix : ndarray of shape (n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape (n, n_classes), optional
        Augmented class labels in one-hot encoded format. Available if ``labels`` are
        provided. One can call ``aug_labels.argmax(axis=1)`` to get class indices.

    Raises
    ------
    ValueError
        If taxa are unavailable.

    See Also
    --------
    mixup

    Notes
    -----
    The Phylomix method was described in [1]_.

    This method leverages phylogenetic relationships to guide data augmentation in
    microbiome and other omic data. By mixing the abundances of evolutionarily related
    taxa (tips of a selected node), Phylomix preserves the biological structure while
    introducing new synthetic samples.

    The selection of nodes follows a random sampling approach, where a subset of taxa
    is chosen based on a Beta-distributed mixing coefficient. This ensures that the
    augmented data maintains biologically meaningful compositional relationships.

    In the original paper, the authors assumed a bifurcated phylogenetic tree, but this
    implementation works with any tree structure. If desired, the user can bifurcate
    the tree using :meth:`~skbio.tree.TreeNode.bifurcate` before augmentation.

    Phylomix is particularly valuable for microbiome-trait association studies, where
    preserving phylogenetic similarity between related taxa is crucial for accurate
    downstream predictions. This approach helps address the common challenge of limited
    sample sizes in omic data studies.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.table import phylomix
    >>> from skbio.tree import TreeNode
    >>> matrix = np.arange(20).reshape(4, 5)
    >>> labels = np.array([0, 1, 0, 1])
    >>> tree = TreeNode.read(['(((a,b),c),(d,e));'])
    >>> taxa = ['a', 'b', 'c', 'd', 'e']
    >>> aug_matrix, aug_labels = phylomix(
    ...     matrix, n=5, tree=tree, taxa=taxa, labels=labels)
    >>> print(aug_matrix.shape)
    (5, 5)
    >>> print(aug_labels.shape)
    (5, 2)

    References
    ----------
    .. [1] Jiang, Y., Liao, D., Zhu, Q., & Lu, Y. Y. (2025). PhyloMix: Enhancing
       microbiome-trait association prediction through phylogeny-mixing augmentation.
       Bioinformatics, btaf014.

    """
    rng = get_rng(seed)
    matrix, labels, pairs, taxa = _format_input(table, labels, intra_class, taxa=taxa)
    if taxa is None:
        raise ValueError("Taxa must be included in table or explicitly provided.")
    n_taxa = len(taxa)

    pairs_sel = pairs[rng.integers(0, len(pairs), size=n)]

    # Pre-calculate feature indices descending from each node.
    node_indices = _indices_under_nodes(tree, taxa)
    n_nodes = len(node_indices)

    # Draw mixing coefficients (lambda) from a beta distribution.
    lams = rng.beta(alpha, alpha, size=n)

    # Randomly shuffle about half of the pairs. This is because each pair of samples
    # are not treated symmetrically in the algorithm. Meanwhile, pairs were generated
    # with the smaller index on the left. If shuffling isn't performed, the algorithm
    # will bias toward one side of the list of samples.
    gates = rng.random(n) < 0.5
    pairs_sel[gates] = pairs_sel[gates][:, ::-1]

    aug_matrix = []
    for (idx1, idx2), lam in zip(pairs_sel, lams):
        x1, x2 = matrix[idx1], matrix[idx2]

        # Number of features to select (ratio: 1 - lambda).
        n_selected = int(np.ceil((1.0 - lam) * n_taxa))

        # Shuffle nodes.
        it = iter(rng.permutation(n_nodes))

        # Randomly sample nodes and add their descending feature indices to selection,
        # until it reaches a given size.
        # Note 1: Nodes may overlap each other, thus descending feature indices may be
        # redundant.
        # Note 2: Dictionary is used instead of set because the order of elements in a
        # dictionary is stable (unlike set, which is unstable).
        # See `_indices_under_nodes` for how it was done.
        # Note 3: One cannot predict how many nodes need to be sampled in order to get
        # n_selected features. Otherwise one can do the following:
        #
        #     nodes_sel = rng.choice(n_nodes, size=num_nodes_to_sample)
        #     selected = np.unique(np.fromiter(itertools.chain.from_iterable(
        #         node_indices[i].keys() for i in nodes_self), dtype=np.intp))
        #
        # Note 4: It it guaranteed that n_selected <= n_taxa <= n_taxa_in_tree,
        # therefore the `while` loop won't last infinitely.
        selected: dict[str, None] = {}
        while len(selected) < n_selected:
            selected.update(node_indices[next(it)])

        # Drop out extra indices.
        # Permutation is used instead of sorting because we just want stable selection
        # of indices whereas order doesn't matter. `matrix[selected]` should always be
        # the same regardless of order.
        selected = np.fromiter(selected.keys(), dtype=np.intp)
        selected = rng.permutation(selected)[:n_selected]

        # Take selected features from sample 2, scale them to match sample 1, then
        # replace those in sample 1.
        # Note that samples 1 and 2 are not treated equally. The synthetic sample will
        # will have the same scale as sample 1.
        x1_sel, x2_sel = x1[selected], x2[selected]
        total1, total2 = x1_sel.sum(), x2_sel.sum()
        if total1 > 0 and total2 > 0:
            # removed astype(int) so phylomix could handle compositional input
            # (total1 / total2) is guaranteed to be float, therefore
            # leaf_counts2 * (total1 / total2) is also guaranteed to be float
            x_mix = x1.copy()
            x_mix[selected] = x2_sel * (total1 / total2)
            aug_matrix.append(x_mix)

        # If either is zero, just use sample 1.
        else:
            aug_matrix.append(x1)

    aug_matrix = np.vstack(aug_matrix)

    aug_labels = _make_aug_labels(labels, pairs_sel, lams.reshape(-1, 1), intra_class)
    return _format_output(aug_matrix, aug_labels, matrix, labels, append)
