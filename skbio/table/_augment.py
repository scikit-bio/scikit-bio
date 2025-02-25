# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.table import Table
from skbio.tree import TreeNode
from skbio._base import SkbioObject
from skbio.util import get_rng
import numpy as np
import pandas as pd
import warnings


class Augmentation(SkbioObject):
    """Data Augmentation of a table,
    for predictive microbiome models.

    Parameters
    ----------
    table : skbio.table.Table
        The table to augment.
    label : numpy.ndarray, optional
        The label of the table.
    tree : skbio.tree.TreeNode, optional
        The tree to use to augment the table.
    """

    def __init__(self, table, label=None, tree=None):
        if table is not None and not isinstance(table, Table):
            raise ValueError("table must be a skbio.table.Table")
        self.table = table

        if tree is not None and not isinstance(tree, TreeNode):
            raise ValueError("tree must be a skbio.tree.TreeNode")
        self.tree = tree

        self.dataframe = self.table.to_dataframe()
        self.matrix = self.dataframe.values
        self.label = label
        self.one_hot_label = np.eye(2)[self.label]

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
        r"""Get all possible pairs of samples that can be used for augmentation
        if intra_class is True, only return pairs of samples within the same class
        if intra_class is False, return all possible pairs of samples.

        Parameters
        ----------
        intra_class : bool
            Whether to get all possible pairs of samples within the same class.

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
        the two labels of the two samples

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
        >>> observ_ids = ['O%d' % i for i in range(10)]
        >>> table = Table(data, observ_ids, sample_ids)
        >>> label = np.random.randint(0, 2, size=4)
        >>> augmentation = Augmentation(table, label)
        >>> aug_matrix, aug_label = augmentation.mixup(n_samples=5)
        >>> print(aug_matrix.shape)
        (10, 9)
        >>> print(aug_label.shape)
        (9, 2)

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.

        """
        rng = get_rng(seed)
        possible_pairs = self._get_all_possible_pairs()
        selected_pairs = possible_pairs[
            rng.integers(0, len(possible_pairs), size=n_samples)
        ]

        indices1 = selected_pairs[:, 0]
        indices2 = selected_pairs[:, 1]
        lambdas = rng.beta(alpha, alpha, size=(1, n_samples))
        augmented_x = (
            lambdas * self.matrix[:, indices1]
            + (1 - lambdas) * self.matrix[:, indices2]
        )
        if self.label is not None:
            augmented_y = (
                lambdas.T * self.one_hot_label[indices1]
                + (1 - lambdas).T * self.one_hot_label[indices2]
            )
            augmented_matrix = np.concatenate([self.matrix, augmented_x], axis=1)
            augmented_label = np.concatenate([self.one_hot_label, augmented_y])
        else:
            augmented_matrix = np.concatenate([self.matrix, augmented_x], axis=1)
            augmented_label = None

        return augmented_matrix, augmented_label

    def aitchison_mixup(self, n_samples, alpha=2, seed=None):
        r"""Data Augmentation by Aitchison mixup.
        it requires the data to be compositional, if the
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

            .. versionadded:: 0.6.3

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.

        Examples
        --------
        >>> from skbio.table import Table
        >>> from skbio.table import Augmentation
        >>> data = np.arange(40).reshape(10, 4)
        >>> sample_ids = ['S%d' % i for i in range(4)]
        >>> observ_ids = ['O%d' % i for i in range(10)]
        >>> table = Table(data, observ_ids, sample_ids)
        >>> table_compositional = table.norm(axis="sample")
        >>> label = np.random.randint(0, 2, size=4)
        >>> augmentation = Augmentation(table_compositional, label)
        >>> aug_matrix, aug_label = augmentation.aitchison_mixup(n_samples=5)
        >>> print(aug_matrix.shape)
        (10, 9)
        >>> print(aug_label.shape)
        (9, 2)

        Notes
        -----
        The algorithm is based on [1]_, and leverages the Aitchison geometry
        to guide data augmentation in compositional microbiome data,
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
        if not np.allclose(np.sum(self.matrix, axis=0), 1):
            warnings.warn("The data is not compositional, it will be normalized first.")
            if np.any(np.sum(self.matrix, axis=0) == 0):
                self.matrix = self.matrix + 1e-5
            self.matrix = self.matrix / np.sum(self.matrix, axis=0)

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
                self._aitchison_scalar_multiplication(_lambda, self.matrix[:, idx1]),
                self._aitchison_scalar_multiplication(
                    1 - _lambda, self.matrix[:, idx2]
                ),
            )
            if self.label is not None:
                augmented_y = (
                    _lambda * self.one_hot_label[idx1]
                    + (1 - _lambda) * self.one_hot_label[idx2]
                )
            else:
                augmented_y = None
            augmented_matrix.append(augmented_x)
            augmented_label.append(augmented_y)
        augmented_matrix = np.array(augmented_matrix).T
        if self.label is not None:
            augmented_matrix = np.concatenate([self.matrix, augmented_matrix], axis=1)
            augmented_label = np.concatenate([self.one_hot_label, augmented_label])
        else:
            augmented_matrix = np.concatenate([self.matrix, augmented_matrix], axis=1)
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
            The augmented label, the label is in 0 and 1.

        Examples
        --------
        >>> from skbio.table import Table
        >>> from skbio.table import Augmentation
        >>> data = np.arange(40).reshape(10, 4)
        >>> sample_ids = ['S%d' % i for i in range(4)]
        >>> observ_ids = ['O%d' % i for i in range(10)]
        >>> table = Table(data, observ_ids, sample_ids)
        >>> label = np.random.randint(0, 2, size=4)
        >>> augmentation = Augmentation(table, label)
        >>> aug_matrix, aug_label = augmentation.compositional_cutmix(n_samples=5)
        >>> print(aug_matrix.shape)
        (10, 9)
        >>> print(aug_label.shape)
        (9,)

        Notes
        -----

        The algorithm is described in [1]_,
        This method needs to do cutmix on data in the same class.
        by randomly select count from one of two samples to generate
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
        possible_pairs = self._get_all_possible_pairs(intra_class=True)
        selected_pairs = possible_pairs[
            rng.integers(0, len(possible_pairs), size=n_samples)
        ]

        augmented_matrix = []
        augmented_label = []
        for idx1, idx2 in selected_pairs:
            x1, x2 = self.matrix[:, idx1], self.matrix[:, idx2]
            _lambda = rng.uniform(0, 1)
            indicator_binomial = rng.binomial(1, _lambda, size=self.matrix.shape[0])
            augmented_x = x1 * indicator_binomial + x2 * (1 - indicator_binomial)
            augmented_matrix.append(augmented_x)
            label = self.label[idx1]
            augmented_label.append(label)

        augmented_matrix = np.array(augmented_matrix).T
        augmented_label = np.array(augmented_label)
        augmented_matrix = np.concatenate([self.matrix, augmented_matrix], axis=1)
        augmented_label = np.concatenate([self.label, augmented_label])
        return augmented_matrix, augmented_label

    def phylomix(self, tip_to_obs_mapping, n_samples, alpha=2, seed=None):
        r"""Data Augmentation by phylomix.

        Parameters
        ----------
        tip_to_obs_mapping : dict
            A dictionary mapping tips to observation indices.
        n_samples : int
            The number of new samples to generate.
        alpha : float
            The alpha parameter of the beta distribution.
        seed : int, Generator or RandomState, optional
            A user-provided random seed or random generator instance. See
            :func:`details <skbio.util.get_rng>`.

            .. versionadded:: 0.6.3

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.

        Examples
        --------
        >>> from skbio.table import Table
        >>> from skbio.table import Augmentation
        >>> data = np.arange(10).reshape(5, 2)
        >>> sample_ids = ['S%d' % i for i in range(2)]
        >>> observ_ids = ['O%d' % i for i in range(5)]
        >>> tree = TreeNode.read(["(((a,b)int1,c)int2,(x,y)int3);"])
        >>> table = Table(data, observ_ids, sample_ids)
        >>> label = np.random.randint(0, 2, size=2)
        >>> aug = Augmentation(table, label, tree)
        >>> tip_to_obs_mapping = {'a': 0, 'b': 1, 'c': 2, 'x': 3, 'y': 4}
        >>> aug_matrix, aug_label = aug.phylomix(tip_to_obs_mapping, n_samples=5)
        >>> print(aug_matrix.shape)
        (5, 7)
        >>> print(aug_label.shape)
        (7, 2)

        Notes
        -----
        The algorithm is based on [1]_, and leverages phylogenetic
        relationships to guide data augmentation in compositional microbiome data.
        By mixing the counts of the leaves of a selected node, Phylomix preserves
        phylogenetic structure while introducing new synthetic samples.

        The selection of nodes follows a random sampling approach,
        where a subset of leaves is chosen based on a
        Beta-distributed mixing coefficient. This ensures that the augmented
        data maintains meaningful compositional relationships.

        Phylomix is particularly useful for microbiome-trait association studies,
        where the preservation of phylogenetic similarity is crucial for
        accurate downstream predictions.

        The method assumes that all tips in the phylogenetic tree
        are represented in the ``tip_to_obs_mapping`` dictionary,
        and that the tree is bifurcated before augmentation.

        References
        ----------
        .. [1] Jiang, Y., Liao, D., Zhu, Q., & Lu, Y. Y. (2025).
            PhyloMix: Enhancing microbiome-trait association prediction through
            phylogeny-mixing augmentation. Bioinformatics, btaf014.

        """
        rng = get_rng(seed)
        if self.tree is None:
            raise ValueError("tree is required for phylomix augmentation")

        self.tree.bifurcate()
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
            x1, x2 = self.matrix[:, pair[0]], self.matrix[:, pair[1]]
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
            augment_label = (
                _lambda * self.one_hot_label[pair[0]]
                + (1 - _lambda) * self.one_hot_label[pair[1]]
            )

            augmented_matrix.append(mixed_x)
            augmented_label.append(augment_label)

        if self.label is not None:
            augmented_matrix = np.array(augmented_matrix).T
            augmented_label = np.array(augmented_label)
            augmented_matrix = np.concatenate([self.matrix, augmented_matrix], axis=1)
            augmented_label = np.concatenate([self.one_hot_label, augmented_label])
        else:
            augmented_matrix = np.concatenate(
                [self.matrix, np.array(augmented_matrix).T], axis=1
            )
            augmented_label = None

        return augmented_matrix, augmented_label
