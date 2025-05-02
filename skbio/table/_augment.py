# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import warnings
import numpy as np
import pandas as pd
from skbio.table import Table
from skbio.tree import TreeNode
from skbio._base import SkbioObject
from skbio.stats.composition import closure
from skbio.util import get_rng
from skbio.util.config._dispatcher import _ingest_array


def _validate_tree(tree):
    if tree is not None and not isinstance(tree, TreeNode):
        raise TypeError("`tree` must be a skbio.tree.TreeNode object.")


def _validate_label(label, matrix):
    if not isinstance(label, np.ndarray):
        raise ValueError(
            f"label must be a numpy.ndarray, but got {type(label)} instead."
        )
    if label.ndim not in [1, 2]:
        raise ValueError(
            f"labels should have shape (n_samples,) or (n_samples, n_classes)"
            f"but got {label.shape} instead."
        )
    if label.ndim == 1:  # and num_classes is None:
        if label.dtype != int:
            raise TypeError(f"label must only contain integer values.")
        # check that labels is 0 indexed
        unique_labels = np.unique(label)
        if min(unique_labels) != 0:
            raise ValueError("Labels must be zero-indexed. Minimum value must " "be 0.")
        num_classes = len(unique_labels)
        exp_labels = np.arange(num_classes)
        # check that label is consecutive integers starting at 0
        if not np.array_equal(unique_labels, exp_labels):
            raise ValueError(
                "Labels must be consecutive integers from 0 " "to num_classes - 1."
            )
        one_hot_label = np.eye(num_classes, dtype=int)[label]
        return label, one_hot_label
    if label.ndim == 2:
        # sanity checks to ensure valid one hot encoding.
        # all rows should sum to 1
        if not all(x == 1 for x in label.sum(axis=1)):
            raise ValueError(
                (
                    "label is not properly one hot encoded. Rows "
                    "(samples) were found with more than one label."
                )
            )
        # sum of all values in label should match number of samples
        if not (label_samples := label.sum(axis=0).sum()) == matrix.shape[1]:
            raise ValueError(
                (
                    f"The number of samples represented by label "
                    f"{label_samples} does not match the number of "
                    f"samples in data {matrix.shape[1]}"
                )
            )
        return None, label


def _aitchison_addition(x, v):
    r"""Perform Aitchison addition on two samples x and v.

    Parameters
    ----------
    x : numpy.ndarray
        The first sample.
    v : numpy.ndarray
        The second sample.

    Returns
    -------
    numpy.ndarray
        The result of Aitchison addition.
    """
    return (xv := x * v) / np.sum(xv)


def _aitchison_scalar_multiplication(lam, x):
    r"""Perform Aitchison multiplication on sample x, with scalar lambda.

    Parameters
    ----------
    lam : float
        The scalar to multiply the sample by.
    x : numpy.ndarray
        The sample to multiply.

    Returns
    -------
    numpy.ndarray
        The result of Aitchison multiplication.

    """
    return (x_to_lam := x**lam) / np.sum(x_to_lam)


def _get_all_possible_pairs(matrix, label=None, intra_class=False):
    r"""Get all possible pairs of samples that can be used for augmentation.

    Parameters
    ----------
    intra_class : bool
        If ``True``, only return pairs of samples within the same class. This
        functionality is only available for binary classification.
        If ``False``, return all possible pairs of samples.

    Returns
    -------
    numpy.ndarray
        An array of all possible pairs of samples.

    """
    possible_pairs = []
    if intra_class:
        if label is None:
            raise ValueError("label is required for intra-class augmentation")
        # so is this only for binary?!
        matrix_cls0_indices = np.where(label == 0)[0]
        matrix_cls1_indices = np.where(label == 1)[0]
        for idx1 in range(matrix_cls0_indices.shape[0]):
            for idx2 in range(idx1 + 1, matrix_cls0_indices.shape[0]):
                possible_pairs.append(
                    (matrix_cls0_indices[idx1], matrix_cls0_indices[idx2])
                )
        for idx1 in range(matrix_cls1_indices.shape[0]):
            for idx2 in range(idx1 + 1, matrix_cls1_indices.shape[0]):
                possible_pairs.append(
                    (matrix_cls1_indices[idx1], matrix_cls1_indices[idx2])
                )
    else:
        n_samples = matrix.shape[0]
        for idx1 in range(n_samples):
            for idx2 in range(idx1 + 1, n_samples):
                possible_pairs.append((idx1, idx2))
    return np.array(possible_pairs)


def mixup(table, n_samples, label=None, alpha=2, seed=None):
    r"""Data Augmentation by vanilla mixup.

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
    n_samples : int
        The number of new samples to generate.
    alpha : float
        The alpha parameter of the beta distribution.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Examples
    --------
    >>> from skbio.table import Table
    >>> from skbio.table import Augmentation
    >>> data = np.arange(40).reshape(10, 4)
    >>> sample_ids = ['S%d' % i for i in range(4)]
    >>> feature_ids = ['O%d' % i for i in range(10)]
    >>> table = Table(data, feature_ids, sample_ids)
    >>> label = np.random.randint(0, 2, size=table.shape[1])
    >>> augmentation = Augmentation(table, label)
    >>> aug_matrix, aug_label = augmentation.mixup(n_samples=5)
    >>> print(aug_matrix.shape)
    (9, 10)
    >>> print(aug_label.shape)
    (9, 2)

    Returns
    -------
    augmented_matrix : numpy.ndarray
        The augmented matrix.
    augmented_label : numpy.ndarray
        The augmented label, in one-hot encoding.
        If the user wants to use the augmented label for regression,
        users can simply call ``np.argmax(aug_label, axis=1)``
        to get the discrete labels.

    Notes
    -----
    The mixup is based on [1]_, and shares the same core concept as PyTorch's
    `MixUp <https://pytorch.org/vision/
    main/generated/torchvision.transforms.v2.MixUp.html>`_.
    there are key differences:

    1. This implementation generates new samples to augment a dataset,
        while PyTorch's MixUp is applied on-the-fly during training
        to batches of data.

    2. This implementation randomly selects pairs of samples from the entire
        dataset, while PyTorch's implementation typically mixes consecutive
        samples in a batch (requiring prior shuffling).

    3. This implementation returns an augmented dataset with both original and
        new samples, while PyTorch's implementation transforms a batch in-place.

    4. This implementation is designed for omic data tables,
        while PyTorch's is primarily for image data.
        And this implementation is mainly based on the Numpy Library.

    References
    ----------
    .. [1] Zhang, H., Cisse, M., Dauphin, Y. N., & Lopez-Paz, D. (2017).
        mixup: Beyond Empirical Risk Minimization.
        arXiv preprint arXiv:1710.09412.

    """
    matrix, row_ids, col_ids = _ingest_array(table)
    label, one_hot_label = _validate_label(label, matrix)
    rng = get_rng(seed)
    possible_pairs = _get_all_possible_pairs(matrix=matrix, label=label)

    selected_pairs = possible_pairs[
        rng.integers(0, len(possible_pairs), size=n_samples)
    ]

    indices1 = selected_pairs[:, 0]
    indices2 = selected_pairs[:, 1]
    lambdas = rng.beta(alpha, alpha, size=(n_samples, 1))
    augmented_x = lambdas * matrix[indices1] + (1 - lambdas) * matrix[indices2]

    augmented_matrix = np.concatenate([matrix, augmented_x], axis=0)
    if label is not None:
        augmented_y = (
            lambdas * one_hot_label[indices1] + (1 - lambdas) * one_hot_label[indices2]
        )
        augmented_label = np.concatenate([one_hot_label, augmented_y])
    else:
        augmented_label = None

    return augmented_matrix, augmented_label


def aitchison_mixup(table, n_samples, label=None, alpha=2, seed=None):
    r"""Data Augmentation by Aitchison mixup.

    This function requires the data to be compositional. If the
    table is not normalized, it will be normalized first.

    Parameters
    ----------
    n_samples : int
        The number of new samples to generate.
    alpha : float
        The alpha parameter of the beta distribution.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    augmented_matrix : numpy.ndarray
        The augmented matrix.
    augmented_label : numpy.ndarray
        The augmented label, in one-hot encoding.
        if the user want to use the augmented label for regression,
        users can simply call ``np.argmax(aug_label, axis=1)``
        to get the discrete labels.

    Examples
    --------
    >>> from skbio.table import Table
    >>> from skbio.table import Augmentation
    >>> data = np.arange(40).reshape(10, 4)
    >>> sample_ids = ['S%d' % i for i in range(4)]
    >>> feature_ids = ['O%d' % i for i in range(10)]
    >>> table = Table(data, feature_ids, sample_ids)
    >>> table_compositional = table.norm(axis="sample")
    >>> label = np.random.randint(0, 2, size=table.shape[1])
    >>> augmentation = Augmentation(table_compositional, label)
    >>> aug_matrix, aug_label = augmentation.aitchison_mixup(n_samples=5)
    >>> print(aug_matrix.shape)
    (9, 10)
    >>> print(aug_label.shape)
    (9, 2)

    Notes
    -----
    The algorithm is based on [1]_, and leverages the Aitchison geometry
    to guide data augmentation in compositional data,
    this is essentially the vanilla mixup in the Aitchison geometry.
    This mixup method only works on the Compositional data.
    where a set of datapoints are living in the simplex:
    :math:`x_i > 0`, and :math:`\sum_{i=1}^{p} x_i = 1`.
    The augmented sample is computed as the linear combination of
    the two samples in the Aitchison geometry. In the Aitchision
    Geometry, we define the addition and scalar multiplication as:

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

    The label is computed as the linear combination of
    the two labels of the two samples

    .. math::

        y = \lambda \cdot y_1 + (1 - \lambda) \cdot y_2

    By mixing the counts of two samples, Aitchison mixup preserves the
    compositional nature of the data, and the sum-to-one property.

    References
    ----------
    .. [1] Gordon-Rodriguez, E., Quinn, T., & Cunningham, J. P. (2022).
        Data augmentation for compositional data: Advancing predictive
        models of the microbiome. Advances in Neural Information Processing
        Systems, 35, 20551-20565.

    """
    matrix, row_ids, col_ids = _ingest_array(table)
    label, one_hot_label = _validate_label(label, matrix)

    if not np.allclose(np.sum(matrix, axis=1), 1):
        matrix = closure(matrix)

    rng = get_rng(seed)
    possible_pairs = _get_all_possible_pairs(matrix)
    selected_pairs = possible_pairs[
        rng.integers(0, len(possible_pairs), size=n_samples)
    ]

    augmented_matrix = []
    augmented_label = []
    for idx1, idx2 in selected_pairs:
        _lambda = rng.beta(alpha, alpha)
        augmented_x = _aitchison_addition(
            _aitchison_scalar_multiplication(_lambda, matrix[idx1]),
            _aitchison_scalar_multiplication(1 - _lambda, matrix[idx2]),
        )
        if label is not None:
            augmented_y = (
                _lambda * one_hot_label[idx1] + (1 - _lambda) * one_hot_label[idx2]
            )
            augmented_label.append(augmented_y)
        augmented_matrix.append(augmented_x)
    augmented_matrix = np.concatenate([matrix, np.array(augmented_matrix)], axis=0)
    if label is not None:
        augmented_label = np.concatenate([one_hot_label, augmented_label])
    else:
        augmented_label = None

    return augmented_matrix, augmented_label


def compositional_cutmix(table, n_samples, label=None, seed=None):
    r"""Data Augmentation by compositional cutmix.

    Parameters
    ----------
    n_samples : int
        The number of new samples to generate.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    augmented_matrix : numpy.ndarray
        The augmented matrix.
    augmented_label : numpy.ndarray
        The augmented label, the label is 1D array.
        User can use the 1D label for both classification and regression.

    Examples
    --------
    >>> from skbio.table import Table
    >>> from skbio.table import Augmentation
    >>> data = np.arange(40).reshape(10, 4)
    >>> sample_ids = ['S%d' % i for i in range(4)]
    >>> feature_ids = ['O%d' % i for i in range(10)]
    >>> table = Table(data, feature_ids, sample_ids)
    >>> label = np.random.randint(0, 2, size=4)
    >>> augmentation = Augmentation(table, label)
    >>> aug_matrix, aug_label = augmentation.compositional_cutmix(n_samples=5)
    >>> print(aug_matrix.shape)
    (9, 10)
    >>> print(aug_label.shape)
    (9,)

    Notes
    -----

    The algorithm is described in [1]_,
    This method needs to do cutmix on compositional data in the same class.
    by randomly selecting counts from one of two samples to generate
    a new sample. For this method to work, the label must be provided.
    The algorithm has 4 steps:

    1. Draw a class :math:`c` from the class prior
    and draw :math:`\lambda \sim Uniform(0, 1)`

    2. Draw two training points :math:`i_1, i_2` from the training set
    such that :math:`y_{i_1} = y_{i_2} = c`, uniformly at random

    3. For each :math:`j \in \{1, ..., p\}`, draw :math:`I_j \sim Binomial(\lambda)`
    and set :math:`\tilde{x}_j = x_{i_1j}` if :math:`I_j = 1`,
    and :math:`\tilde{x}_j = x_{i_2j}` if :math:`I_j = 0`

    4. Set :math:`\tilde{y} = c`

    References
    ----------
    .. [1] Gordon-Rodriguez, E., Quinn, T., & Cunningham, J. P. (2022).
        Data augmentation for compositional data: Advancing predictive
        models of the microbiome. Advances in Neural Information Processing
        Systems, 35, 20551-20565.

    """

    rng = get_rng(seed)

    matrix, row_ids, col_ids = _ingest_array(table)
    label, one_hot_label = _validate_label(label, matrix)

    if not np.allclose(np.sum(matrix, axis=1), 1):
        matrix = closure(matrix)

    possible_pairs = _get_all_possible_pairs(matrix, label=label, intra_class=True)
    selected_pairs = possible_pairs[
        rng.integers(0, len(possible_pairs), size=n_samples)
    ]

    augmented_matrix = []
    augmented_label = []
    for idx1, idx2 in selected_pairs:
        x1, x2 = matrix[idx1], matrix[idx2]
        _lambda = rng.uniform(0, 1)
        indicator_binomial = rng.binomial(1, _lambda, size=matrix.shape[1])
        augmented_x = x1 * indicator_binomial + x2 * (1 - indicator_binomial)
        augmented_matrix.append(augmented_x)
        label_ = label[idx1]
        augmented_label.append(label_)

    augmented_matrix = np.array(augmented_matrix)
    augmented_label = np.array(augmented_label)
    augmented_matrix = np.concatenate([matrix, augmented_matrix], axis=0)
    augmented_label = np.concatenate([label, augmented_label])
    return augmented_matrix, augmented_label


def phylomix(
    table, tree, tip_to_obs_mapping, n_samples, label=None, alpha=2, seed=None
):
    r"""Data augmentation by phylomix.

    Parameters
    ----------
    tip_to_obs_mapping : dict
        A dictionary mapping tips to feature indices.
    n_samples : int
        The number of new samples to generate.
    alpha : float
        The alpha parameter of the beta distribution.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.


    Returns
    -------
    augmented_matrix : numpy.ndarray
        The augmented matrix.
    augmented_label : numpy.ndarray
        The augmented label, in one-hot encoding.
        if the user want to use the augmented label for regression,
        users can simply call ``np.argmax(aug_label, axis=1)``
        to get the discrete labels.

    Examples
    --------
    >>> from skbio.table import Table
    >>> from skbio.table import Augmentation
    >>> data = np.arange(10).reshape(5, 2)
    >>> sample_ids = ['S%d' % i for i in range(2)]
    >>> feature_ids = ['O%d' % i for i in range(5)]
    >>> tree = TreeNode.read(["(((a,b)int1,c)int2,(x,y)int3);"])
    >>> table = Table(data, feature_ids, sample_ids)
    >>> label = np.random.randint(0, 2, size=2)
    >>> aug = Augmentation(table, label, tree=tree)
    >>> tip_to_obs_mapping = {'a': 0, 'b': 1, 'c': 2, 'x': 3, 'y': 4}
    >>> aug_matrix, aug_label = aug.phylomix(tip_to_obs_mapping, n_samples=5)
    >>> print(aug_matrix.shape)
    (7, 5)
    >>> print(aug_label.shape)
    (7, 2)

    Notes
    -----
    The algorithm is based on [1]_, and leverages phylogenetic
    relationships to guide data augmentation in microbiome and other omic data.
    By mixing the abundances of phylogenetically related
    taxa (leaves of a selected node), Phylomix preserves the biological
    structure while introducing new synthetic samples.

    The selection of nodes follows a random sampling approach,
    where a subset of taxa is chosen based on a
    Beta-distributed mixing coefficient. This ensures that the augmented
    data maintains biologically meaningful compositional relationships.

    In the original paper, the authors assumed a bifurcated phylogenetic tree,
    but this implementation works with any tree structure. If desired,
    users can bifurcate their tree using ``skbio.tree.TreeNode.bifurcate()``
    before augmentation.

    Phylomix is particularly valuable for microbiome-trait association studies,
    where preserving phylogenetic similarity between related taxa is crucial for
    accurate downstream predictions. This approach helps address the
    common challenge of limited sample sizes in omic data studies.

    The method assumes that all tips in the phylogenetic tree
    are represented in the ``tip_to_obs_mapping`` dictionary.

    References
    ----------
    .. [1] Jiang, Y., Liao, D., Zhu, Q., & Lu, Y. Y. (2025).
        PhyloMix: Enhancing microbiome-trait association prediction through
        phylogeny-mixing augmentation. Bioinformatics, btaf014.

    """
    rng = get_rng(seed)
    matrix, row_ids, col_ids = _ingest_array(table)
    label, one_hot_label = _validate_label(label, matrix)
    _validate_tree(tree)

    leave_names = [tip.name for tip in tree.tips()]

    if set(tip_to_obs_mapping.keys()) != set(leave_names):
        raise ValueError("tip_to_obs_mapping must contain all tips in the tree")

    # Convert nodes to indices for random selection
    all_nodes = [node for node in tree.levelorder()]
    node_indices = np.arange(len(all_nodes))
    num_leaves = len(leave_names)

    possible_pairs = _get_all_possible_pairs(matrix)
    selected_pairs = possible_pairs[
        rng.integers(0, len(possible_pairs), size=n_samples)
    ]
    feature_dict = {feature_name: idx for idx, feature_name in enumerate(tree.tips())}
    augmented_matrix = []
    augmented_label = []
    for pair in selected_pairs:
        x1, x2 = matrix[pair[0]], matrix[pair[1]]
        _lambda = rng.beta(alpha, alpha)
        n_leaves = int(np.ceil((1 - _lambda) * num_leaves))
        selected_index = set()
        mixed_x = x1.copy()

        while len(selected_index) < n_leaves:
            # Select a random node using index
            node_idx = rng.choice(node_indices)
            available_node = all_nodes[node_idx]
            leaf_idx = [feature_dict[leaf] for leaf in available_node.tips()]
            obs_idx = [tip_to_obs_mapping[leaf.name] for leaf in available_node.tips()]
            selected_index.update(obs_idx)

        selected_index = rng.choice(list(selected_index), n_leaves, replace=False)

        leaf_counts1, leaf_counts2 = (
            x1[selected_index].astype(np.float32),
            x2[selected_index].astype(np.float32),
        )

        total1, total2 = leaf_counts1.sum(), leaf_counts2.sum()
        if total1 > 0 and total2 > 0:
            leaf_counts2_normalized = leaf_counts2 / total2
            new_counts = (total1 * leaf_counts2_normalized).astype(int)
            mixed_x[selected_index] = new_counts
        else:
            mixed_x[selected_index] = leaf_counts1

        if label is not None:
            augment_label = (
                _lambda * one_hot_label[pair[0]]
                + (1 - _lambda) * one_hot_label[pair[1]]
            )
            augmented_label.append(augment_label)
        augmented_matrix.append(mixed_x)

    if label is not None:
        augmented_matrix = np.array(augmented_matrix)
        augmented_label = np.array(augmented_label)
        augmented_matrix = np.concatenate([matrix, augmented_matrix], axis=0)
        augmented_label = np.concatenate([one_hot_label, augmented_label])
    else:
        augmented_matrix = np.concatenate([matrix, np.array(augmented_matrix)], axis=0)
        augmented_label = None

    return augmented_matrix, augmented_label


class Augmentation(SkbioObject):
    """Data Augmentation of a omic data table, for predictive models on omic data.

    The Augmentation class is mainly designed to enhance the performance of
    classification models on omic data.

    Parameters
    ----------
    table : skbio.table.Table
        The table to augment.
    label : numpy.ndarray of int, optional
        The label of the table. The label is expected to has a shape of ``(n_samples,)``
        or ``(n_samples, n_classes)`` if one hot encoding has been performed. Can be
        none if the label is not needed. Values in label must be integer, and
        correspond to class labels for model training and prediction.
    tree : skbio.tree.TreeNode, optional
        The tree to use to augment the table. Only required for method ``phylomix``.
    """

    def __init__(self, table, label=None, tree=None):
        if table is not None and not isinstance(table, Table):
            raise ValueError("table must be a skbio.table.Table")
        self.table = table

        if tree is not None and not isinstance(tree, TreeNode):
            raise ValueError("tree must be a skbio.tree.TreeNode")
        self.tree = tree

        self.dataframe = self.table.to_dataframe()
        self.matrix = self.dataframe.values.T
        self.label = label
        if label is not None:
            if not isinstance(label, np.ndarray):
                raise ValueError(
                    f"label must be a numpy.ndarray, but got {type(label)} instead."
                )
            if label.ndim not in [1, 2]:
                raise ValueError(
                    f"labels should have shape (n_samples,) or (n_samples, n_classes)"
                    f"but got {label.shape} instead."
                )
            if label.ndim == 1:  # and num_classes is None:
                if label.dtype != int:
                    raise TypeError(f"label must only contain integer values.")
                # check that labels is 0 indexed
                unique_labels = np.unique(label)
                if min(unique_labels) != 0:
                    raise ValueError(
                        "Labels must be zero-indexed. Minimum value must " "be 0."
                    )
                num_classes = len(unique_labels)
                exp_labels = np.arange(num_classes)
                # check that label is consecutive integers starting at 0
                if not np.array_equal(unique_labels, exp_labels):
                    raise ValueError(
                        "Labels must be consecutive integers from 0 "
                        "to num_classes - 1."
                    )
                self.one_hot_label = np.eye(num_classes, dtype=int)[self.label]
            if label.ndim == 2:
                # sanity checks to ensure valid one hot encoding.
                # all rows should sum to 1
                if not all(x == 1 for x in label.sum(axis=1)):
                    raise ValueError(
                        (
                            "label is not properly one hot encoded. Rows "
                            "(samples) were found with more than one label."
                        )
                    )
                # sum of all values in label should match number of samples
                if (
                    not (label_samples := label.sum(axis=0).sum())
                    == self.matrix.shape[1]
                ):
                    raise ValueError(
                        (
                            f"The number of samples represented by label "
                            f"{label_samples} does not match the number of "
                            f"samples in data {self.matrix.shape[1]}"
                        )
                    )
                self.one_hot_label = label

    def __str__(self):
        return f"Augmentation(shape={self.matrix.shape})"

    def _aitchison_addition(self, x, v):
        r"""Perform Aitchison addition on two samples x and v.

        Parameters
        ----------
        x : numpy.ndarray
            The first sample.
        v : numpy.ndarray
            The second sample.

        Returns
        -------
        numpy.ndarray
            The result of Aitchison addition.
        """
        return (xv := x * v) / np.sum(xv)

    def _aitchison_scalar_multiplication(self, lam, x):
        r"""Perform Aitchison multiplication on sample x, with scalar lambda.

        Parameters
        ----------
        lam : float
            The scalar to multiply the sample by.
        x : numpy.ndarray
            The sample to multiply.

        Returns
        -------
        numpy.ndarray
            The result of Aitchison multiplication.

        """
        return (x_to_lam := x**lam) / np.sum(x_to_lam)

    def _get_all_possible_pairs(self, intra_class=False):
        r"""Get all possible pairs of samples that can be used for augmentation.

        Parameters
        ----------
        intra_class : bool
            If ``True``, only return pairs of samples within the same class.
            If ``False``, return all possible pairs of samples.

        Returns
        -------
        numpy.ndarray
            An array of all possible pairs of samples.

        """
        possible_pairs = []
        if intra_class:
            if self.label is None:
                raise ValueError("label is required for intra-class augmentation")
            matrix_cls0_indices = np.where(self.label == 0)[0]
            matrix_cls1_indices = np.where(self.label == 1)[0]
            for idx1 in range(matrix_cls0_indices.shape[0]):
                for idx2 in range(idx1 + 1, matrix_cls0_indices.shape[0]):
                    possible_pairs.append(
                        (matrix_cls0_indices[idx1], matrix_cls0_indices[idx2])
                    )
            for idx1 in range(matrix_cls1_indices.shape[0]):
                for idx2 in range(idx1 + 1, matrix_cls1_indices.shape[0]):
                    possible_pairs.append(
                        (matrix_cls1_indices[idx1], matrix_cls1_indices[idx2])
                    )
        else:
            for idx1 in range(self.table.shape[1]):
                for idx2 in range(idx1 + 1, self.table.shape[1]):
                    possible_pairs.append((idx1, idx2))
        return np.array(possible_pairs)

    def mixup(self, n_samples, alpha=2, seed=None):
        r"""Data Augmentation by vanilla mixup.

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
        n_samples : int
            The number of new samples to generate.
        alpha : float
            The alpha parameter of the beta distribution.
        seed : int, Generator or RandomState, optional
            A user-provided random seed or random generator instance. See
            :func:`details <skbio.util.get_rng>`.

        Examples
        --------
        >>> from skbio.table import Table
        >>> from skbio.table import Augmentation
        >>> data = np.arange(40).reshape(10, 4)
        >>> sample_ids = ['S%d' % i for i in range(4)]
        >>> feature_ids = ['O%d' % i for i in range(10)]
        >>> table = Table(data, feature_ids, sample_ids)
        >>> label = np.random.randint(0, 2, size=table.shape[1])
        >>> augmentation = Augmentation(table, label)
        >>> aug_matrix, aug_label = augmentation.mixup(n_samples=5)
        >>> print(aug_matrix.shape)
        (9, 10)
        >>> print(aug_label.shape)
        (9, 2)

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.
            If the user wants to use the augmented label for regression,
            users can simply call ``np.argmax(aug_label, axis=1)``
            to get the discrete labels.

        Notes
        -----
        The mixup is based on [1]_, and shares the same core concept as PyTorch's
        `MixUp <https://pytorch.org/vision/
        main/generated/torchvision.transforms.v2.MixUp.html>`_.
        there are key differences:

        1. This implementation generates new samples to augment a dataset,
           while PyTorch's MixUp is applied on-the-fly during training
           to batches of data.

        2. This implementation randomly selects pairs of samples from the entire
           dataset, while PyTorch's implementation typically mixes consecutive
           samples in a batch (requiring prior shuffling).

        3. This implementation returns an augmented dataset with both original and
           new samples, while PyTorch's implementation transforms a batch in-place.

        4. This implementation is designed for omic data tables,
           while PyTorch's is primarily for image data.
           And this implementation is mainly based on the Numpy Library.

        References
        ----------
        .. [1] Zhang, H., Cisse, M., Dauphin, Y. N., & Lopez-Paz, D. (2017).
            mixup: Beyond Empirical Risk Minimization.
            arXiv preprint arXiv:1710.09412.

        """
        rng = get_rng(seed)
        possible_pairs = self._get_all_possible_pairs()
        selected_pairs = possible_pairs[
            rng.integers(0, len(possible_pairs), size=n_samples)
        ]

        indices1 = selected_pairs[:, 0]
        indices2 = selected_pairs[:, 1]
        lambdas = rng.beta(alpha, alpha, size=(n_samples, 1))
        augmented_x = (
            lambdas * self.matrix[indices1] + (1 - lambdas) * self.matrix[indices2]
        )

        augmented_matrix = np.concatenate([self.matrix, augmented_x], axis=0)
        if self.label is not None:
            augmented_y = (
                lambdas * self.one_hot_label[indices1]
                + (1 - lambdas) * self.one_hot_label[indices2]
            )
            augmented_label = np.concatenate([self.one_hot_label, augmented_y])
        else:
            augmented_label = None

        return augmented_matrix, augmented_label

    def aitchison_mixup(self, n_samples, alpha=2, seed=None):
        r"""Data Augmentation by Aitchison mixup.

        This function requires the data to be compositional. If the
        table is not normalized, it will be normalized first.

        Parameters
        ----------
        n_samples : int
            The number of new samples to generate.
        alpha : float
            The alpha parameter of the beta distribution.
        seed : int, Generator or RandomState, optional
            A user-provided random seed or random generator instance. See
            :func:`details <skbio.util.get_rng>`.

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.
            if the user want to use the augmented label for regression,
            users can simply call ``np.argmax(aug_label, axis=1)``
            to get the discrete labels.

        Examples
        --------
        >>> from skbio.table import Table
        >>> from skbio.table import Augmentation
        >>> data = np.arange(40).reshape(10, 4)
        >>> sample_ids = ['S%d' % i for i in range(4)]
        >>> feature_ids = ['O%d' % i for i in range(10)]
        >>> table = Table(data, feature_ids, sample_ids)
        >>> table_compositional = table.norm(axis="sample")
        >>> label = np.random.randint(0, 2, size=table.shape[1])
        >>> augmentation = Augmentation(table_compositional, label)
        >>> aug_matrix, aug_label = augmentation.aitchison_mixup(n_samples=5)
        >>> print(aug_matrix.shape)
        (9, 10)
        >>> print(aug_label.shape)
        (9, 2)

        Notes
        -----
        The algorithm is based on [1]_, and leverages the Aitchison geometry
        to guide data augmentation in compositional data,
        this is essentially the vanilla mixup in the Aitchison geometry.
        This mixup method only works on the Compositional data.
        where a set of datapoints are living in the simplex:
        :math:`x_i > 0`, and :math:`\sum_{i=1}^{p} x_i = 1`.
        The augmented sample is computed as the linear combination of
        the two samples in the Aitchison geometry. In the Aitchision
        Geometry, we define the addition and scalar multiplication as:

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

        The label is computed as the linear combination of
        the two labels of the two samples

        .. math::

            y = \lambda \cdot y_1 + (1 - \lambda) \cdot y_2

        By mixing the counts of two samples, Aitchison mixup preserves the
        compositional nature of the data, and the sum-to-one property.

        References
        ----------
        .. [1] Gordon-Rodriguez, E., Quinn, T., & Cunningham, J. P. (2022).
            Data augmentation for compositional data: Advancing predictive
            models of the microbiome. Advances in Neural Information Processing
            Systems, 35, 20551-20565.

        """
        if not np.allclose(np.sum(self.matrix, axis=1), 1):
            self.matrix = closure(self.matrix)

        rng = get_rng(seed)
        possible_pairs = self._get_all_possible_pairs()
        selected_pairs = possible_pairs[
            rng.integers(0, len(possible_pairs), size=n_samples)
        ]

        augmented_matrix = []
        augmented_label = []
        for idx1, idx2 in selected_pairs:
            _lambda = rng.beta(alpha, alpha)
            augmented_x = self._aitchison_addition(
                self._aitchison_scalar_multiplication(_lambda, self.matrix[idx1]),
                self._aitchison_scalar_multiplication(1 - _lambda, self.matrix[idx2]),
            )
            if self.label is not None:
                augmented_y = (
                    _lambda * self.one_hot_label[idx1]
                    + (1 - _lambda) * self.one_hot_label[idx2]
                )
                augmented_label.append(augmented_y)
            augmented_matrix.append(augmented_x)
        augmented_matrix = np.concatenate(
            [self.matrix, np.array(augmented_matrix)], axis=0
        )
        if self.label is not None:
            augmented_label = np.concatenate([self.one_hot_label, augmented_label])
        else:
            augmented_label = None

        return augmented_matrix, augmented_label

    def compositional_cutmix(self, n_samples, seed=None):
        r"""Data Augmentation by compositional cutmix.

        Parameters
        ----------
        n_samples : int
            The number of new samples to generate.
        seed : int, Generator or RandomState, optional
            A user-provided random seed or random generator instance. See
            :func:`details <skbio.util.get_rng>`.

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, the label is 1D array.
            User can use the 1D label for both classification and regression.

        Examples
        --------
        >>> from skbio.table import Table
        >>> from skbio.table import Augmentation
        >>> data = np.arange(40).reshape(10, 4)
        >>> sample_ids = ['S%d' % i for i in range(4)]
        >>> feature_ids = ['O%d' % i for i in range(10)]
        >>> table = Table(data, feature_ids, sample_ids)
        >>> label = np.random.randint(0, 2, size=4)
        >>> augmentation = Augmentation(table, label)
        >>> aug_matrix, aug_label = augmentation.compositional_cutmix(n_samples=5)
        >>> print(aug_matrix.shape)
        (9, 10)
        >>> print(aug_label.shape)
        (9,)

        Notes
        -----

        The algorithm is described in [1]_,
        This method needs to do cutmix on compositional data in the same class.
        by randomly selecting counts from one of two samples to generate
        a new sample. For this method to work, the label must be provided.
        The algorithm has 4 steps:

        1. Draw a class :math:`c` from the class prior
        and draw :math:`\lambda \sim Uniform(0, 1)`

        2. Draw two training points :math:`i_1, i_2` from the training set
        such that :math:`y_{i_1} = y_{i_2} = c`, uniformly at random

        3. For each :math:`j \in \{1, ..., p\}`, draw :math:`I_j \sim Binomial(\lambda)`
        and set :math:`\tilde{x}_j = x_{i_1j}` if :math:`I_j = 1`,
        and :math:`\tilde{x}_j = x_{i_2j}` if :math:`I_j = 0`

        4. Set :math:`\tilde{y} = c`

        References
        ----------
        .. [1] Gordon-Rodriguez, E., Quinn, T., & Cunningham, J. P. (2022).
            Data augmentation for compositional data: Advancing predictive
            models of the microbiome. Advances in Neural Information Processing
            Systems, 35, 20551-20565.

        """

        rng = get_rng(seed)

        if not np.allclose(np.sum(self.matrix, axis=1), 1):
            self.matrix = closure(self.matrix)

        possible_pairs = self._get_all_possible_pairs(intra_class=True)
        selected_pairs = possible_pairs[
            rng.integers(0, len(possible_pairs), size=n_samples)
        ]

        augmented_matrix = []
        augmented_label = []
        for idx1, idx2 in selected_pairs:
            x1, x2 = self.matrix[idx1], self.matrix[idx2]
            _lambda = rng.uniform(0, 1)
            indicator_binomial = rng.binomial(1, _lambda, size=self.matrix.shape[1])
            augmented_x = x1 * indicator_binomial + x2 * (1 - indicator_binomial)
            augmented_matrix.append(augmented_x)
            label = self.label[idx1]
            augmented_label.append(label)

        augmented_matrix = np.array(augmented_matrix)
        augmented_label = np.array(augmented_label)
        augmented_matrix = np.concatenate([self.matrix, augmented_matrix], axis=0)
        augmented_label = np.concatenate([self.label, augmented_label])
        return augmented_matrix, augmented_label

    def phylomix(self, tip_to_obs_mapping, n_samples, alpha=2, seed=None):
        r"""Data augmentation by phylomix.

        Parameters
        ----------
        tip_to_obs_mapping : dict
            A dictionary mapping tips to feature indices.
        n_samples : int
            The number of new samples to generate.
        alpha : float
            The alpha parameter of the beta distribution.
        seed : int, Generator or RandomState, optional
            A user-provided random seed or random generator instance. See
            :func:`details <skbio.util.get_rng>`.


        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.
            if the user want to use the augmented label for regression,
            users can simply call ``np.argmax(aug_label, axis=1)``
            to get the discrete labels.

        Examples
        --------
        >>> from skbio.table import Table
        >>> from skbio.table import Augmentation
        >>> data = np.arange(10).reshape(5, 2)
        >>> sample_ids = ['S%d' % i for i in range(2)]
        >>> feature_ids = ['O%d' % i for i in range(5)]
        >>> tree = TreeNode.read(["(((a,b)int1,c)int2,(x,y)int3);"])
        >>> table = Table(data, feature_ids, sample_ids)
        >>> label = np.random.randint(0, 2, size=2)
        >>> aug = Augmentation(table, label, tree=tree)
        >>> tip_to_obs_mapping = {'a': 0, 'b': 1, 'c': 2, 'x': 3, 'y': 4}
        >>> aug_matrix, aug_label = aug.phylomix(tip_to_obs_mapping, n_samples=5)
        >>> print(aug_matrix.shape)
        (7, 5)
        >>> print(aug_label.shape)
        (7, 2)

        Notes
        -----
        The algorithm is based on [1]_, and leverages phylogenetic
        relationships to guide data augmentation in microbiome and other omic data.
        By mixing the abundances of phylogenetically related
        taxa (leaves of a selected node), Phylomix preserves the biological
        structure while introducing new synthetic samples.

        The selection of nodes follows a random sampling approach,
        where a subset of taxa is chosen based on a
        Beta-distributed mixing coefficient. This ensures that the augmented
        data maintains biologically meaningful compositional relationships.

        In the original paper, the authors assumed a bifurcated phylogenetic tree,
        but this implementation works with any tree structure. If desired,
        users can bifurcate their tree using ``skbio.tree.TreeNode.bifurcate()``
        before augmentation.

        Phylomix is particularly valuable for microbiome-trait association studies,
        where preserving phylogenetic similarity between related taxa is crucial for
        accurate downstream predictions. This approach helps address the
        common challenge of limited sample sizes in omic data studies.

        The method assumes that all tips in the phylogenetic tree
        are represented in the ``tip_to_obs_mapping`` dictionary.

        References
        ----------
        .. [1] Jiang, Y., Liao, D., Zhu, Q., & Lu, Y. Y. (2025).
            PhyloMix: Enhancing microbiome-trait association prediction through
            phylogeny-mixing augmentation. Bioinformatics, btaf014.

        """
        rng = get_rng(seed)
        if self.tree is None:
            raise ValueError("tree is required for phylomix augmentation")

        leave_names = [tip.name for tip in self.tree.tips()]

        if set(tip_to_obs_mapping.keys()) != set(leave_names):
            raise ValueError("tip_to_obs_mapping must contain all tips in the tree")

        # Convert nodes to indices for random selection
        all_nodes = [node for node in self.tree.levelorder()]
        node_indices = np.arange(len(all_nodes))
        num_leaves = len(leave_names)

        possible_pairs = self._get_all_possible_pairs()
        selected_pairs = possible_pairs[
            rng.integers(0, len(possible_pairs), size=n_samples)
        ]
        feature_dict = {
            feature_name: idx for idx, feature_name in enumerate(self.tree.tips())
        }
        augmented_matrix = []
        augmented_label = []
        for pair in selected_pairs:
            x1, x2 = self.matrix[pair[0]], self.matrix[pair[1]]
            _lambda = rng.beta(alpha, alpha)
            n_leaves = int(np.ceil((1 - _lambda) * num_leaves))
            selected_index = set()
            mixed_x = x1.copy()

            while len(selected_index) < n_leaves:
                # Select a random node using index
                node_idx = rng.choice(node_indices)
                available_node = all_nodes[node_idx]
                leaf_idx = [feature_dict[leaf] for leaf in available_node.tips()]
                obs_idx = [
                    tip_to_obs_mapping[leaf.name] for leaf in available_node.tips()
                ]
                selected_index.update(obs_idx)

            selected_index = rng.choice(list(selected_index), n_leaves, replace=False)

            leaf_counts1, leaf_counts2 = (
                x1[selected_index].astype(np.float32),
                x2[selected_index].astype(np.float32),
            )

            total1, total2 = leaf_counts1.sum(), leaf_counts2.sum()
            if total1 > 0 and total2 > 0:
                leaf_counts2_normalized = leaf_counts2 / total2
                new_counts = (total1 * leaf_counts2_normalized).astype(int)
                mixed_x[selected_index] = new_counts
            else:
                mixed_x[selected_index] = leaf_counts1

            if self.label is not None:
                augment_label = (
                    _lambda * self.one_hot_label[pair[0]]
                    + (1 - _lambda) * self.one_hot_label[pair[1]]
                )
                augmented_label.append(augment_label)
            augmented_matrix.append(mixed_x)

        if self.label is not None:
            augmented_matrix = np.array(augmented_matrix)
            augmented_label = np.array(augmented_label)
            augmented_matrix = np.concatenate([self.matrix, augmented_matrix], axis=0)
            augmented_label = np.concatenate([self.one_hot_label, augmented_label])
        else:
            augmented_matrix = np.concatenate(
                [self.matrix, np.array(augmented_matrix)], axis=0
            )
            augmented_label = None

        return augmented_matrix, augmented_label
