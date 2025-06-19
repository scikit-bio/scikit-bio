# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, Union, TYPE_CHECKING
from itertools import chain, combinations

if TYPE_CHECKING:
    from numpy.random import RandomState, Generator
    from numpy.typing import NDArray

import numpy as np
from skbio.tree import TreeNode
from skbio.stats.composition import closure
from skbio.util import get_rng
from skbio.util.config._dispatcher import _ingest_array, _create_table


def _validate_tree(tree: TreeNode) -> None:
    """Ensure that tree is a TreeNode object."""
    if not isinstance(tree, TreeNode):
        raise TypeError("`tree` must be a skbio.tree.TreeNode object.")


def _validate_label(  # type: ignore[return]
    label: "NDArray", matrix: "NDArray"
) -> tuple["NDArray", "NDArray"]:
    """Ensure that the provided label is appropriate for the provided matrix.

    Parameters
    ----------
    label : ndarray
        The class labels for the data. The label is expected to have shape
        ``(samples,)`` or ``(samples, n_classes)``.
    matrix : ndarray
        The data matrix of shape ``(samples, n_features).``

    Returns
    -------
    label_1d : ndarray
        The class labels in 1D format with shape ``(samples,)``. Contains
        integer class labels from 0 to n_classes-1. If the input was already
        1D, this is the same as the input. If the input was one-hot encoded,
        this is the argmax reconstruction.
    label_one_hot : ndarray
        The class labels in one-hot encoded format with shape
        ``(samples, n_classes)``. Each row contains exactly one 1 and the
        rest 0s, indicating the class membership for that sample.

    """
    if not isinstance(label, np.ndarray):
        raise ValueError(
            f"label must be a numpy.ndarray, but got {type(label)} instead."
        )
    if label.ndim not in [1, 2]:
        raise ValueError(
            f"labels should have shape (samples,) or (samples, n_classes)"
            f"but got {label.shape} instead."
        )
    if not label.shape[0] == matrix.shape[0]:
        raise ValueError(
            f"Number of elements in label ({label.shape[0]}) does not "
            f"match number of samples in input data ({matrix.shape[0]})"
        )
    # make sure label is only whole numbers, then convert to int
    if not np.all(np.equal(np.mod(label, 1), 0)):
        raise TypeError(f"label must only contain integer values.")
    if label.ndim == 1:
        # check that labels is 0 indexed
        unique_labels = np.unique(label)
        if min(unique_labels) != 0:
            raise ValueError("Labels must be zero-indexed. Minimum value must be 0.")
        num_classes = len(unique_labels)
        exp_labels = np.arange(num_classes)
        # check that label is consecutive integers starting at 0
        if not np.array_equal(unique_labels, exp_labels):
            raise ValueError(
                "Labels must be consecutive integers from 0 to num_classes - 1."
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
        if not (label_samples := label.sum(axis=0).sum()) == matrix.shape[0]:
            raise ValueError(
                (
                    f"The number of samples represented by label "
                    f"({label_samples}) does not match the number of "
                    f"samples in data ({matrix.shape[0]})"
                )
            )
        # do this for consistent return values whether a one hot encoded or 1D label
        # is provided
        label_1d = np.argmax(label, axis=1)
        return label_1d, label


def _aitchison_add(x: "NDArray", v: "NDArray") -> "NDArray":
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
    return (xv := x * v) / xv.sum()


def _aitchison_scalar_multiply(lam: float, x: "NDArray") -> "NDArray":
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
    return (x_to_lam := x**lam) / x_to_lam.sum()


def _get_all_possible_pairs(
    matrix: "NDArray", label: Optional["NDArray"] = None, intra_class: bool = False
) -> np.ndarray:
    r"""Get all possible pairs of samples that can be used for augmentation.

    Parameters
    ----------
    matrix : ndarray
        Sample by features array of the input data.
    label : ndarray, optional
        The 1D class labels for the data.
    intra_class : bool
        If ``True``, only return pairs of samples within the same class. This
        functionality is only available for binary classification.
        If ``False``, return all possible pairs of samples.

    Returns
    -------
    ndarray
        An array of all possible pairs of samples.

    """
    possible_pairs = []
    # Intra-class pair generation is currently only possible for binary classification.
    if intra_class:
        if label is None:
            raise ValueError("Label is required for intra-class augmentation.")
        possible_pairs = list(
            chain(
                combinations(np.where(label == 0)[0], 2),
                combinations(np.where(label == 1)[0], 2),
            )
        )
    else:
        possible_pairs = list(combinations(np.arange(matrix.shape[0]), 2))
    return np.array(possible_pairs)


def mixup(
    table,
    samples: int,
    label: Optional["NDArray"] = None,
    alpha: float = 2,
    normalize: bool = False,
    seed: Optional[Union[int, "Generator", "RandomState"]] = None,
    output_format: Optional[str] = None,
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
    table : table_like
        Samples by features table (n, m). See the `DataTable <https://scikit.bio/
        docs/dev/generated/skbio.util.config.html#the-datatable-type>`_ type
        documentation for details.
    samples : int
        The number of new samples to generate.
    label : ndarray
        The label of the table. The label is expected to has a shape of ``(samples,)``
        or ``(samples, n_classes)``.
    alpha : float
        The alpha parameter of the beta distribution.
    normalize : bool, optional
        If ``True`` and the input is not already compositional, scikit-bio's
        :func:`~skbio.stats.composition.closure` function will be called, ensuring
        values for each sample add up to 1. Defaults to ``False``.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.
    output_format : str, optional
        Standard ``DataTable`` parameter. See the `DataTable <https://scikit.bio/
        docs/dev/generated/skbio.util.config.html#the-datatable-type>`_ type
        documentation for details.

    Examples
    --------
    >>> from skbio.table import mixup
    >>> data = np.arange(40).reshape(4, 10)
    >>> label = np.array([0, 1, 0, 1])
    >>> aug_matrix, aug_label = mixup(data, label=label, samples=5)
    >>> print(aug_matrix.shape)
    (9, 10)
    >>> print(aug_label.shape)
    (9, 2)

    Returns
    -------
    augmented_matrix : table_like
        The augmented matrix.
    augmented_label : table_like
        The augmented label, in one-hot encoding.
        If the user wants to use the augmented label for regression,
        users can simply call ``np.argmax(aug_label, axis=1)``
        to get the discrete labels.

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

    References
    ----------
    .. [1] Zhang, H., Cisse, M., Dauphin, Y. N., & Lopez-Paz, D. (2017). mixup: Beyond
       Empirical Risk Minimization. arXiv preprint arXiv:1710.09412.

    """
    matrix, row_ids, col_ids = _ingest_array(table)
    if normalize:
        if not np.allclose(np.sum(matrix, axis=1), 1):
            matrix = closure(matrix)
    if label is not None:
        label, one_hot_label = _validate_label(label, matrix)
    rng = get_rng(seed)
    possible_pairs = _get_all_possible_pairs(matrix=matrix, label=label)

    selected_pairs = possible_pairs[rng.integers(0, len(possible_pairs), size=samples)]

    indices1 = selected_pairs[:, 0]
    indices2 = selected_pairs[:, 1]
    lambdas = rng.beta(alpha, alpha, size=(samples, 1))
    augmented_x = lambdas * matrix[indices1] + (1 - lambdas) * matrix[indices2]

    augmented_matrix = np.concatenate([matrix, augmented_x], axis=0)
    if label is not None:
        augmented_y = (
            lambdas * one_hot_label[indices1] + (1 - lambdas) * one_hot_label[indices2]
        )
        augmented_label = np.concatenate([one_hot_label, augmented_y])
    else:
        augmented_label = None

    return _create_table(augmented_matrix, backend=output_format), _create_table(
        augmented_label, backend=output_format
    )


def aitchison_mixup(
    table,
    samples: int,
    label: Optional["NDArray"] = None,
    alpha: float = 2,
    normalize: bool = True,
    seed: Optional[Union[int, "Generator", "RandomState"]] = None,
    output_format: Optional[str] = None,
) -> tuple:
    r"""Data augmentation by Aitchison mixup.

    This function requires the data to be compositional. If the
    table is not normalized, it will be normalized first.

    Parameters
    ----------
    table : table_like
        Samples by features table (n, m). See the `DataTable <https://scikit.bio/
        docs/dev/generated/skbio.util.config.html#the-datatable-type>`_ type
        documentation for details.
    samples : int
        The number of new samples to generate.
    label : ndarray
        The label of the table. The label is expected to has a shape of ``(samples,)``
        or ``(samples, n_classes)``.
    alpha : float
        The alpha parameter of the beta distribution.
    normalize : bool, optional
        If ``True`` and the input is not already compositional, scikit-bio's
        :func:`~skbio.stats.composition.closure` function will be called, ensuring
        values for each sample add up to 1. Defaults to ``True``.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.
    output_format : str, optional
        Standard ``DataTable`` parameter. See the `DataTable <https://scikit.bio/
        docs/dev/generated/skbio.util.config.html#the-datatable-type>`_ type
        documentation for details.

    Returns
    -------
    augmented_matrix : table_like
        The augmented matrix.
    augmented_label : table_like
        The augmented label, in one-hot encoding.
        if the user want to use the augmented label for regression,
        users can simply call ``np.argmax(aug_label, axis=1)``
        to get the discrete labels.

    See Also
    --------
    mixup

    Examples
    --------
    >>> from skbio.table import aitchison_mixup
    >>> data = np.arange(40).reshape(4, 10)
    >>> label = np.array([0, 1, 0, 1])
    >>> aug_matrix, aug_label = aitchison_mixup(data, label=label, samples=5)
    >>> print(aug_matrix.shape)
    (9, 10)
    >>> print(aug_label.shape)
    (9, 2)

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

    References
    ----------
    .. [1] Gordon-Rodriguez, E., Quinn, T., & Cunningham, J. P. (2022). Data
       augmentation for compositional data: Advancing predictive models of the
       microbiome. Advances in Neural Information Processing Systems, 35,
       20551-20565.

    """
    matrix, _, _ = _ingest_array(table)
    if label is not None:
        label, one_hot_label = _validate_label(label, matrix)

    if normalize:
        if not np.allclose(np.sum(matrix, axis=1), 1):
            matrix = closure(matrix)

    rng = get_rng(seed)
    possible_pairs = _get_all_possible_pairs(matrix)
    selected_pairs = possible_pairs[rng.integers(0, len(possible_pairs), size=samples)]

    augmented_matrix = []
    augmented_label = []
    for idx1, idx2 in selected_pairs:
        _lambda = rng.beta(alpha, alpha)
        augmented_x = _aitchison_add(
            _aitchison_scalar_multiply(_lambda, matrix[idx1]),
            _aitchison_scalar_multiply(1 - _lambda, matrix[idx2]),
        )
        if label is not None:
            augmented_y = (
                _lambda * one_hot_label[idx1] + (1 - _lambda) * one_hot_label[idx2]
            )
            augmented_label.append(augmented_y)
        augmented_matrix.append(augmented_x)
    augmented_matrix = np.concatenate([matrix, np.array(augmented_matrix)], axis=0)
    if label is not None:
        augmented_label_ = np.concatenate([one_hot_label, augmented_label])
    else:
        augmented_label_ = None
    return _create_table(augmented_matrix, backend=output_format), _create_table(
        augmented_label_, backend=output_format
    )


def compositional_cutmix(
    table,
    samples: int,
    label: Optional["NDArray"] = None,
    normalize: bool = True,
    seed: Optional[Union[int, "Generator", "RandomState"]] = None,
    output_format: Optional[str] = None,
) -> tuple:
    r"""Data augmentation by compositional cutmix.

    Parameters
    ----------
    table : table_like
        Samples by features table (n, m). See the `DataTable <https://scikit.bio/
        docs/dev/generated/skbio.util.config.html#the-datatable-type>`_ type
        documentation for details.
    samples : int
        The number of new samples to generate.
    label : ndarray, optional
        The label of the table. The label is expected to has a shape of
        ``(samples,)`` or ``(samples, n_classes)``.
    normalize : bool, optional
        If ``True`` and the input is not already compositional, scikit-bio's
        :func:`~skbio.stats.composition.closure` function will be called, ensuring
        values for each sample add up to 1. Defaults to ``True``.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.
    output_format : str, optional
        Standard ``DataTable`` parameter. See the `DataTable <https://scikit.bio/
        docs/dev/generated/skbio.util.config.html#the-datatable-type>`_ type
        documentation for details.

    Returns
    -------
    augmented_matrix : table_like
        The augmented matrix.
    augmented_label : table_like
        The augmented label, the label is 1D array.
        User can use the 1D label for both classification and regression.

    Examples
    --------
    >>> from skbio.table import compositional_cutmix
    >>> data = np.arange(40).reshape(4, 10)
    >>> label = np.array([0, 1, 0, 1])
    >>> aug_matrix, aug_label = compositional_cutmix(data, label=label, samples=5)
    >>> print(aug_matrix.shape)
    (9, 10)
    >>> print(aug_label.shape)
    (9, 2)

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

    References
    ----------
    .. [1] Gordon-Rodriguez, E., Quinn, T., & Cunningham, J. P. (2022). Data
       augmentation for compositional data: Advancing predictive models of the
       microbiome. Advances in Neural Information Processing Systems, 35,
       20551-20565.

    """
    rng = get_rng(seed)

    matrix, _, _ = _ingest_array(table)
    # If label isn't provided, assume all samples have same class label
    no_out_label = False
    if label is None:
        no_out_label = True
        label = np.zeros(matrix.shape[0], dtype=int)
    label, one_hot_label = _validate_label(label, matrix)

    if normalize:
        if not np.allclose(np.sum(matrix, axis=1), 1):
            matrix = closure(matrix)
    possible_pairs = _get_all_possible_pairs(matrix, label=label, intra_class=True)
    selected_pairs = possible_pairs[rng.integers(0, len(possible_pairs), size=samples)]

    augmented_matrix = []
    augmented_label = []
    for idx1, idx2 in selected_pairs:
        x1, x2 = matrix[idx1], matrix[idx2]
        _lambda = rng.uniform(0, 1)
        indicator_binomial = rng.binomial(1, _lambda, size=matrix.shape[1])
        augmented_x = x1 * indicator_binomial + x2 * (1 - indicator_binomial)
        augmented_matrix.append(augmented_x)
        label_ = one_hot_label[idx1]
        augmented_label.append(label_)

    augmented_matrix_ = np.array(augmented_matrix)
    augmented_matrix = np.concatenate([matrix, augmented_matrix_], axis=0)
    if no_out_label:
        augmented_label = None
    else:
        augmented_label_ = np.array(augmented_label)
        augmented_label = np.concatenate([one_hot_label, augmented_label_])
    return _create_table(augmented_matrix, backend=output_format), _create_table(
        augmented_label, backend=output_format
    )


def phylomix(
    table,
    tree: TreeNode,
    tip_to_obs_mapping: dict[str, int],
    samples: int,
    label: Optional["NDArray"] = None,
    alpha: float = 2,
    normalize: bool = False,
    seed: Optional[Union[int, "Generator", "RandomState"]] = None,
    output_format: Optional[str] = None,
) -> tuple:
    r"""Data augmentation by phylomix.

    Parameters
    ----------
    table : table_like
        Samples by features table (n, m). See the `DataTable <https://scikit.bio/
        docs/dev/generated/skbio.util.config.html#the-datatable-type>`_ type
        documentation for details.
    tree : skbio.tree.TreeNode
        The tree to use to augment the table.
    tip_to_obs_mapping : dict
        A dictionary mapping tips to feature indices.
    samples : int
        The number of new samples to generate.
    label : ndarray
        The label of the table. The label is expected to has a shape of ``(samples,)``
        or ``(samples, n_classes)``.
    alpha : float
        The alpha parameter of the beta distribution.
    normalize : bool, optional
        If ``True`` and the input is not already compositional, scikit-bio's
        :func:`~skbio.stats.composition.closure` function will be called, ensuring
        values for each sample add up to 1. Defaults to ``False``.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.
    output_format : str, optional
        Standard ``DataTable`` parameter. See the `DataTable <https://scikit.bio/
        docs/dev/generated/skbio.util.config.html#the-datatable-type>`_ type
        documentation for details.

    Returns
    -------
    augmented_matrix : table_like
        The augmented matrix.
    augmented_label : table_like
        The augmented label, in one-hot encoding.
        if the user want to use the augmented label for regression,
        users can simply call ``np.argmax(aug_label, axis=1)``
        to get the discrete labels.

    Examples
    --------
    >>> from skbio.table import phylomix
    >>> data = np.arange(10).reshape(2, 5)
    >>> tree = TreeNode.read(["(((a,b)int1,c)int2,(x,y)int3);"])
    >>> label = np.array([0, 1])
    >>> tip_to_obs_mapping = {'a': 0, 'b': 1, 'c': 2, 'x': 3, 'y': 4}
    >>> aug_matrix, aug_label = phylomix(data,
    ...                                  tree,
    ...                                  tip_to_obs_mapping,
    ...                                  label=label,
    ...                                  samples=5)
    >>> print(aug_matrix.shape)
    (7, 5)
    >>> print(aug_label.shape)
    (7, 2)

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

    The method assumes that all tips in the phylogenetic tree are represented in the
    ``tip_to_obs_mapping`` dictionary.

    References
    ----------
    .. [1] Jiang, Y., Liao, D., Zhu, Q., & Lu, Y. Y. (2025). PhyloMix: Enhancing
       microbiome-trait association prediction through phylogeny-mixing augmentation.
       Bioinformatics, btaf014.

    """
    rng = get_rng(seed)
    matrix, row_ids, col_ids = _ingest_array(table)

    if normalize:
        if not np.allclose(np.sum(matrix, axis=1), 1):
            matrix = closure(matrix)

    if label is not None:
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
    selected_pairs = possible_pairs[rng.integers(0, len(possible_pairs), size=samples)]
    feature_dict = {feature_name: idx for idx, feature_name in enumerate(tree.tips())}
    augmented_matrix = []
    augmented_label = []
    for pair in selected_pairs:
        x1, x2 = matrix[pair[0]], matrix[pair[1]]
        _lambda = rng.beta(alpha, alpha)
        n_leaves = int(np.ceil((1 - _lambda) * num_leaves))
        selected_index: set[int] = set()
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
            # removed astype(int) so phylomix could handle compositional input
            new_counts = total1 * leaf_counts2_normalized
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
        augmented_matrix_ = np.array(augmented_matrix)
        augmented_label_ = np.array(augmented_label)
        augmented_matrix = np.concatenate([matrix, augmented_matrix_], axis=0)
        augmented_label = np.concatenate([one_hot_label, augmented_label_])
    else:
        augmented_matrix = np.concatenate([matrix, np.array(augmented_matrix)], axis=0)
        augmented_label = None  # type: ignore[assignment]

    return _create_table(augmented_matrix, backend=output_format), _create_table(
        augmented_label, backend=output_format
    )
