# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from numpy.random import RandomState, Generator
    from numpy.typing import ArrayLike, NDArray
    from skbio.tree import TreeNode

import numpy as np
from skbio.util import get_rng
from skbio.table._tabular import _ingest_table, _create_table


def _normalize_matrix(matrix):
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
    from skbio.stats.composition import closure

    return closure(matrix)


def _validate_labels(  # type: ignore[return]
    labels: "NDArray", n: int
) -> tuple["NDArray", "NDArray"]:
    r"""Ensure that the provided label is appropriate for the provided matrix.

    Parameters
    ----------
    labels : array_like of shape (n_samples,) or (n_samples, n_classes)
        Class labels for the data.
    n : int
        Number of samples.

    Returns
    -------
    label_index : ndarray of shape (n_samples,)
        The class labels in 1-D format. Contains indices from 0 to n_classes - 1. If
        the input was already 1-D, this is the same as the input. If the input was
        one-hot encoded, this is the argmax reconstruction.
    label_one_hot : ndarray of shape (n_samples, n_classes)
        The class labels in one-hot encoded format. Each row contains exactly one 1 and
        the rest 0s, indicating the class membership for that sample.

    """
    labels = np.asarray(labels)

    if labels.ndim not in (1, 2):
        raise ValueError(
            f"Label should have shape (n_samples,) or (n_samples, n_classes)"
            f"but got {labels.shape} instead."
        )

    if labels.shape[0] != n:
        raise ValueError(
            f"Number of elements in label ({labels.shape[0]}) does not match "
            f"number of samples in input data ({n})."
        )

    # make sure labels are only whole numbers, then convert to int
    if not np.all(np.equal(np.mod(labels, 1), 0)):
        raise TypeError(f"Label must only contain integer values.")

    # input as indexed labels
    if labels.ndim == 1:
        # check that labels are 0-indexed
        unique_labels = np.unique(labels)
        if unique_labels.min() != 0:
            raise ValueError("Label must be zero-indexed. Minimum value must be 0.")
        n_classes = len(unique_labels)
        exp_labels = np.arange(n_classes)
        # check that label is consecutive integers starting at 0
        if not np.array_equal(unique_labels, exp_labels):
            raise ValueError(
                "Label must be consecutive integers from 0 to n_classes - 1."
            )
        # perform one-hot encoding
        labels_one_hot = np.eye(n_classes, dtype=int)[labels]
        return labels, labels_one_hot

    # input as one-hot encoded labels
    elif labels.ndim == 2:
        # all rows should sum to 1
        msg = "Labels are not properly one-hot encoded."
        if not all(x == 1 for x in labels.sum(axis=1)):
            raise ValueError(
                msg + " Rows (samples) were found with more than one label."
            )
        # sum of all values in label should match number of samples
        if labels.sum() != n:
            raise ValueError(
                msg + " Sum of all values does not match number of samples."
            )
        # generated label indices
        labels_index = labels.argmax(axis=1)
        return labels_index, labels


def _all_pairs(n: int) -> "NDArray":
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


def _intra_class_pairs(labels: "ArrayLike") -> "NDArray":
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


def _format_input(table, labels, intra_class=False, normalize=False, taxa=None):
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
        return matrix, None, pairs, taxa
    labels_1d, labels_2d = _validate_labels(labels, n_samples)
    if intra_class:
        pairs = _intra_class_pairs(labels_1d)
    else:
        pairs = _all_pairs(n_samples)
    if len(pairs) == 0:
        raise ValueError("Cannot find a pair of samples to mix.")
    return matrix, labels_2d, pairs, taxa


def _make_aug_labels(labels, pairs, lams, intra_class):
    r"""Generate class labels for synthetic samples.

    Parameters
    ----------
    labels : ndarray of shape (n, n_classes), or None
        Original class labels in one-hot encoded format.
    pairs : ndarray of shape (n_pairs, 2)
        Pairs of sample indices that were mixed up to create synthetic samples.
    lams : ndarray of shape (n,) or (n, 1)
        Lambda parameters (probabilities) used during mixing-up.
    intra_class : bool
        If True, samples were paired within each class.

    Returns
    -------
    aug_labels : ndarray of shape (n, n_classes), or None
        Synthetic class labels in one-hot encoded format.

    """
    if labels is None:
        return
    if intra_class:
        return labels[pairs[:, 0]]
    if lams.ndim == 1:
        lams = lams.reshape(-1, 1)
    return lams * labels[pairs[:, 0]] + (1 - lams) * labels[pairs[:, 1]]


def _format_output(aug_matrix, aug_labels, matrix, labels, append=False):
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
    if append:
        aug_matrix = np.concatenate([matrix, aug_matrix])
        if aug_labels is not None:
            aug_labels = np.concatenate([labels, aug_labels])
    else:
        aug_matrix = np.asarray(aug_matrix)
        if aug_labels is not None:
            aug_labels = np.asarray(aug_labels)
    return aug_matrix, aug_labels


def mixup(
    table,
    n: int,
    labels: Optional["NDArray"] = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    append: bool = False,
    seed: Optional[Union[int, "Generator", "RandomState"]] = None,
) -> tuple:
    r"""Data augmentation by vanilla mixup.

    Randomly select two samples :math:`s_1` and :math:`s_2` from the OTU table,
    and generate a new sample :math:`s` by a linear combination
    of :math:`s_1` and :math:`s_2`, as follows:

    .. math::

        s = \lambda \cdot s_1 + (1 - \lambda) \cdot s_2

    where :math:`\lambda` is a random number sampled from a beta distribution
    with parameters :math:`\alpha` and :math:`\alpha`.
    The label is computed as the linear combination of
    the two labels of the two samples:

    .. math::

        y = \lambda \cdot y_1 + (1 - \lambda) \cdot y_2

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table to be augmented. See :doc:`../articles/table_like` for
        supported formats.
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
        False (default) If False, only the synthetic samples are returned.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    aug_matrix : ndarray of shape (n_samples + n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape (n_samples + n, n_classes), optional
        Augmented class labels in one-hot encoded format. Available if ``labels`` are
        provided. One can call ``aug_labels.argmax(axis=1)`` to get class indices.

    Notes
    -----
    The mixup method is based on [1]_, and shares the same core concept as PyTorch's
    `MixUp <https://pytorch.org/vision/
    main/generated/torchvision.transforms.v2.MixUp.html>`_ class. There are some key
    differences:

    1. This implementation generates new samples to augment a dataset, while PyTorch's
       MixUp is applied on-the-fly during training to batches of data.

    2. This implementation randomly selects pairs of samples from the entire dataset,
       while PyTorch's implementation typically mixes consecutive samples in a batch
       (requiring prior shuffling).

    3. This implementation returns an augmented dataset with both original and new
       new samples, while PyTorch's implementation transforms a batch in-place.

    4. This implementation is tailored for biological omic data tables, while PyTorch's
       is primarily for image data.

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
    matrix, labels, pairs, _ = _format_input(table, labels, intra_class)

    rng = get_rng(seed)
    pairs_sel = pairs[rng.integers(0, len(pairs), size=n)]

    idx1 = pairs_sel[:, 0]
    idx2 = pairs_sel[:, 1]
    lams = rng.beta(alpha, alpha, size=n).reshape(-1, 1)

    aug_matrix = lams * matrix[idx1] + (1 - lams) * matrix[idx2]

    aug_labels = _make_aug_labels(labels, pairs_sel, lams, intra_class)

    return _format_output(aug_matrix, aug_labels, matrix, labels, append)


def _aitchison_add(x: "NDArray", y: "NDArray") -> "NDArray":
    r"""Perform Aitchison addition on two samples x and y.

    Parameters
    ----------
    x : ndarray of shape (n_features,)
        The first sample.
    y : ndarray of shape (n_features,)
        The second sample.

    Returns
    -------
    ndarray of shape (n_features,)
        The result of Aitchison addition.
    """
    return (xy := x * y) / xy.sum()


def _aitchison_scalar_multiply(lam: float, x: "NDArray") -> "NDArray":
    r"""Perform Aitchison multiplication on sample x, with scalar lambda.

    Parameters
    ----------
    lam : float
        The scalar to multiply the sample by.
    x : ndarray of shape (n_features,)
        The sample to multiply.

    Returns
    -------
    ndarray of shape (n_features,)
        The result of Aitchison multiplication.

    """
    return (x_to_lam := x**lam) / x_to_lam.sum()


def aitchison_mixup(
    table,
    n: int,
    labels: Optional["NDArray"] = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    normalize: bool = True,
    append: bool = False,
    seed: Optional[Union[int, "Generator", "RandomState"]] = None,
) -> tuple:
    r"""Data augmentation by Aitchison mixup.

    This function requires the data to be compositional. If the
    table is not normalized, it will be normalized first.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table to be augmented. See :doc:`../articles/table_like` for
        supported formats.
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
        False (default) If False, only the synthetic samples are returned.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    aug_matrix : ndarray of shape (n_samples + n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape (n_samples + n, n_classes), optional
        Augmented class labels in one-hot encoded format. Available if ``labels`` are
        provided. One can call ``aug_labels.argmax(axis=1)`` to get class indices.

    See Also
    --------
    mixup

    Notes
    -----
    The algorithm is based on [1]_, and leverages the Aitchison geometry to guide the
    augmentation of compositional data. It is essentially the vanilla mixup method in
    the Aitchison space.

    This method only works on compositional data, where a set of datapoints are living
    in the simplex: :math:`x_i > 0`, and :math:`\sum_{i=1}^{p} x_i = 1`. The augmented
    sample is computed as the linear combination of the two samples in the Aitchison
    space. The Aitchison addition and scalar multiplication are defined as:

    .. math::

        \lambda \otimes s =
        \frac{1}{\sum_{j=1}^{p} s_j^{\lambda}}
        (x_1^{\lambda}, x_2^{\lambda}, ..., x_p^{\lambda})

    .. math::

        s \oplus t =
        \frac{1}{\sum_{j=1}^{p} s_j t_j}
        (s_1 t_1, s_2 t_2, ..., s_p t_p)

    .. math::

        s = (\lambda \otimes  s_1) \oplus ((1 - \lambda) \otimes s_2)

    The label is computed as the linear combination of the labels of the two samples.

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
    matrix, labels, pairs, _ = _format_input(table, labels, intra_class, normalize)

    rng = get_rng(seed)
    pairs_sel = pairs[rng.integers(0, len(pairs), size=n)]

    lams = rng.beta(alpha, alpha, size=n)

    # TODO: this can be optimized
    aug_matrix = []
    aug_labels = [] if labels is not None else None
    for (idx1, idx2), lam in zip(pairs_sel, lams):
        aug_matrix.append(
            _aitchison_add(
                _aitchison_scalar_multiply(lam, matrix[idx1]),
                _aitchison_scalar_multiply(1 - lam, matrix[idx2]),
            )
        )

    aug_labels = _make_aug_labels(labels, pairs_sel, lams, intra_class)

    return _format_output(aug_matrix, aug_labels, matrix, labels, append)


def compositional_cutmix(
    table,
    n: int,
    labels: Optional["NDArray"] = None,
    normalize: bool = True,
    append: bool = False,
    seed: Optional[Union[int, "Generator", "RandomState"]] = None,
) -> tuple:
    r"""Data augmentation by compositional cutmix.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table to be augmented. See :doc:`../articles/table_like` for
        supported formats.
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
        False (default) If False, only the synthetic samples are returned.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

        .. notes::
           This function does not have the ``intra_class`` parameter. It always
           operates in intra-class mode.

    Returns
    -------
    aug_matrix : ndarray of shape (n_samples + n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape (n_samples + n, n_classes), optional
        Augmented class labels in one-hot encoded format. Available if ``labels`` are
        provided. One can call ``aug_labels.argmax(axis=1)`` to get class indices.

    Notes
    -----
    The compositional cutmix method was described in [1]_.

    This method randomly selects counts from one of a pair of samples to generate
    a new sample. It has four steps:

    1. Draw a class :math:`c` from the class prior and draw
    :math:`\lambda \sim Uniform(0, 1)`

    2. Draw two training points :math:`i_1, i_2` from the training set such that
    :math:`y_{i_1} = y_{i_2} = c`, uniformly at random.

    3. For each :math:`j \in \{1, ..., p\}`, draw :math:`I_j \sim Binomial(\lambda)`
    and set :math:`\tilde{x}_j = x_{i_1j}` if :math:`I_j = 1`,
    and :math:`\tilde{x}_j = x_{i_2j}` if :math:`I_j = 0`

    4. Set :math:`\tilde{y} = c`

    This method is applied separately to samples of each class. Therefore, it expects
    class labels (``label``) as part of the input. If ``label`` is None, all samples
    will be considered as belonging to the same class, and ``augmented_label`` will be
    returned as None.

    This algorithm currently only works with binary classification problems (i.e., two
    classes), as it requires intra-class generation of possible sample pairs.

    Examples
    --------
    >>> import numpy as np
    >>> from skbio.table import compositional_cutmix
    >>> matrix = np.arange(40).reshape(4, 10)
    >>> labels = np.array([0, 1, 0, 1])
    >>> aug_matrix, aug_labels = compositional_cutmix(matrix, n=5, labels=labels)
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
    matrix, labels, pairs, _ = _format_input(table, labels, True, normalize)

    rng = get_rng(seed)
    pairs_sel = pairs[rng.integers(0, len(pairs), size=n)]

    # Sample binary gates to control the choice between samples 1 and 2 per feature.
    probs = rng.uniform(0, 1, size=(n, 1))
    gates = rng.binomial(1, probs, size=(n, matrix.shape[1])).astype(bool)

    samples_sel = matrix[pairs_sel]
    aug_matrix = np.where(gates, samples_sel[:, 0], samples_sel[:, 1])

    aug_labels = labels[pairs_sel[:, 0]] if labels is not None else None

    return _format_output(aug_matrix, aug_labels, matrix, labels, append)


def _indices_under_nodes(
    tree: "TreeNode", taxon_map: dict[str, int]
) -> list[list[int]]:
    """Generate a nested list representing feature indices descending from each node.

    `taxon_map` is a mapping of taxa (tip names) to feature indices.

    Only internal nodes with at least one descending feature are included.

    """
    res = []
    res_append = res.append
    for node in tree.postorder(include_self=True):
        if node.children:
            lst = []
            for child in node.children:
                lst.extend(child._taxa)
                del child._taxa
            node._taxa = lst
            if lst:
                res_append(lst)
        else:
            if node.name in taxon_map:
                node._taxa = [taxon_map[node.name]]
            else:
                node._taxa = []
    del tree._taxa
    return res


def phylomix(
    table,
    n: int,
    tree: "TreeNode",
    taxa: Optional["ArrayLike"] = None,
    labels: Optional["NDArray"] = None,
    intra_class: bool = False,
    alpha: float = 2.0,
    append: bool = False,
    seed: Optional[Union[int, "Generator", "RandomState"]] = None,
) -> tuple:
    r"""Data augmentation by PhyloMix.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Input data table to be augmented. See :doc:`../articles/table_like` for
        supported formats.
    n : int
        Number of synthetic samples to generate.
    tree : :class:`~skbio.tree.TreeNode`
        Tree structure modeling the relationships between features.
    taxa : array_like of shape (n_features,), optional
        Feature IDs corresponding to tip names (taxa) in ``tree``. Can be omitted if
        ``table`` already contains relevant feature IDs. Otherwise they need to be
        explicitly provided.
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
        False (default) If False, only the synthetic samples are returned.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    aug_matrix : ndarray of shape (n_samples + n, n_features)
        Augmented data matrix.
    aug_labels : ndarray of shape (n_samples + n, n_classes), optional
        Augmented class labels in one-hot encoded format. Available if ``labels`` are
        provided. One can call ``aug_labels.argmax(axis=1)`` to get class indices.

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

    References
    ----------
    .. [1] Jiang, Y., Liao, D., Zhu, Q., & Lu, Y. Y. (2025). PhyloMix: Enhancing
       microbiome-trait association prediction through phylogeny-mixing augmentation.
       Bioinformatics, btaf014.

    """
    matrix, labels, pairs, taxa = _format_input(table, labels, intra_class, taxa=taxa)

    # Map taxa (tip names) to feature indices.
    from skbio.tree._utils import _validate_taxa_and_tree

    _validate_taxa_and_tree(taxa, tree)
    n_taxa = len(taxa)
    taxon_map = {taxon: i for i, taxon in enumerate(taxa)}

    # Pre-calculate feature indices descending from each node.
    node_indices = _indices_under_nodes(tree, taxon_map)
    n_nodes = len(node_indices)

    rng = get_rng(seed)
    pairs_sel = pairs[rng.integers(0, len(pairs), size=n)]

    # Sample probabilities from a beta distribution.
    lams = rng.beta(alpha, alpha, size=n)

    # Randomly shuffle about half of the pairs. This is because each pair of samples
    # are not treated symmetrically in the algorithm. Meanwhile, pairs were generated
    # with the smaller index on the left. If shuffling isn't performed, the algorithm
    # will bias toward one side of the list of samples.
    gates = rng.random(n) < 0.5
    pairs_sel[gates] = pairs_sel[gates][:, ::-1]

    aug_matrix = []
    aug_labels = [] if labels is not None else None
    for (idx1, idx2), lam in zip(pairs_sel, lams):
        x1, x2 = matrix[idx1], matrix[idx2]

        # Number of features to select (ratio: 1 - lambda).
        n_selected = int(np.ceil((1 - lam) * n_taxa))

        # Randomly sample nodes and add their descending feature indices to selection,
        # until it reaches a given size.
        # Note 1: Nodes may overlap each other, thus descending feature indices may be
        # redundant.
        # Note 2: One cannot predict how many nodes need to be sampled in order to get
        # n_selected features. Otherwise one can do the following:
        #     nodes_sel = rng.choice(n_nodes, size=num_nodes_to_sample)
        #     selected = np.unique(np.fromiter(itertools.chain.from_iterable(
        #         node_indices[i] for i in nodes_self), dtype=np.intp))
        # Note 3: A dictionary is used because the order of elements is stable (unlike
        # set, which is unstable and can impact reproducibility.
        # Note 4: Permutation is used instead of sorting because we just want stable
        # selection of indices whereas order doesn't matter. `matrix[indices]` should
        # always be the same.
        selected = {}
        for i in rng.permutation(n_nodes):
            selected.update(dict.fromkeys(node_indices[i]))
            if len(selected) >= n_selected:
                break
        # Note: len(selected) < n_selected is impossible, therefore a check is omitted.
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

    aug_labels = _make_aug_labels(labels, pairs_sel, lams, intra_class)

    return _format_output(aug_matrix, aug_labels, matrix, labels, append)
