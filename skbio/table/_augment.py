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
import numpy as np
import pandas as pd
import random


class Augmentation(SkbioObject):
    """Augmentation of a table with metadata.

    Parameters
    ----------
    table : skbio.table.Table
        The table to augment.
    method : str
        The method to use to augment the table.
    tree [Optional[skbio.tree.TreeNode]]: skbio.tree.TreeNode
        The tree to use to augment the table.
    """

    def __init__(self, table, method, label=None, tree=None):
        self.table = table
        self.method = method
        self.tree = tree
        self.dataframe = self.table.to_dataframe()
        self.matrix = self.dataframe.values
        self.label = label
        if self.label is not None and len(self.label.shape) == 1:
            self.one_hot_label = pd.get_dummies(self.label).values
        else:
            self.one_hot_label = self.label

    def __str__(self):
        return f"Augmentation(method={self.method})"

    def _aitchison_addition(self, x, v):
        r"""Perform Aitchison addition on two samples x and v

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
        sum_xv = np.sum(x * v)
        return (x * v) / sum_xv

    def _aitchison_scalar_multiplication(self, lam, x):
        r"""Perform Aitchison multiplication on sample x, with scal ar lambda

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
        sum_x_times_lam = np.sum(x * lam)
        return (x**lam) / sum_x_times_lam

    def _get_all_possible_pairs(self):
        r"""Get all possible pairs of samples that can be used for augmentation

        Returns
        -------
        numpy.ndarray
            An array of all possible pairs of samples.
        """
        possible_pairs = []
        for idx1 in range(self.table.shape[1]):
            for idx2 in range(idx1 + 1, self.table.shape[1]):
                possible_pairs.append((idx1, idx2))
        return np.array(possible_pairs)

    def mixup(self, n_samples, alpha=2):
        r"""Data Augmentation by vanilla mixup
        Randomly select two samples s1 and s2 from the OTU table,
        and generate a new sample s by a linear combination of s1 and s2.
        s = lambda * s1 + (1 - lambda) * s2
        where lambda is a random number sampled from a beta distribution
        with parameters alpha and alpha.


        Parameters
        ----------
        n_samples : int
            The number of new samples to generate.
        alpha : float
            The alpha parameter of the beta distribution.

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.
        """
        possible_pairs = self._get_all_possible_pairs()
        selected_pairs = possible_pairs[
            np.random.randint(0, len(possible_pairs), size=n_samples)
        ]

        indices1 = selected_pairs[:, 0]
        indices2 = selected_pairs[:, 1]
        lambdas = np.random.beta(alpha, alpha, size=(1, n_samples))
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

    def aitchison_mixup(self, n_samples, alpha=2):
        r"""Implementation of Aitchison mixup,
        which is essentially the vanilla mixup in the Aitchison geometry.
        This mixup method only works on the Compositional data.
        where a set of datapoints are living in the simplex:
        x_i > 0, and \sum_{i=1}^{p} x_i = 1
        See paper:
        Data Augmentation for Compositional Data:
        Advancing Predictive Models of the Microbiome

        Parameters
        ----------
        n_samples : int
            The number of new samples to generate.
        alpha : float
            The alpha parameter of the beta distribution.

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.
        """
        possible_pairs = self._get_all_possible_pairs()
        selected_pairs = possible_pairs[
            np.random.randint(0, len(possible_pairs), size=n_samples)
        ]

        augmented_matrix = []
        augmented_label = []
        for idx1, idx2 in selected_pairs:
            _lambda = np.random.beta(alpha, alpha)
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

    def phylomix(self, tip_to_obs_mapping, n_samples, alpha=2):
        r"""Data Augmentation by phylomix

        Jiang, Y., Liao, D., Zhu, Q., & Lu, Y. Y. (2025).
        PhyloMix: Enhancing microbiome-trait association prediction through
        phylogeny-mixing augmentation. Bioinformatics, btaf014.
        Phylomix use the phylogenetic tree to generate new samples
        by mixing the counts of the leaves of a selected node.

        Parameters
        ----------
        tip_to_obs_mapping : dict
            A dictionary mapping tips to observation indices.
        n_samples : int
            The number of new samples to generate.
        alpha : float
            The alpha parameter of the beta distribution.

        Returns
        -------
        augmented_matrix : numpy.ndarray
            The augmented matrix.
        augmented_label : numpy.ndarray
            The augmented label, in one-hot encoding.
        """

        if self.tree is None:
            raise ValueError("tree is required for phylomix augmentation")

        self.tree.bifurcate()
        leave_names = [tip.name for tip in self.tree.tips()]

        if set(tip_to_obs_mapping.keys()) != set(leave_names):
            raise ValueError("tip_to_obs_mapping must contain all tips in the tree")

        all_nodes = [node for node in self.tree.levelorder()]
        num_leaves = len(leave_names)

        possible_pairs = self._get_all_possible_pairs()
        selected_pairs = possible_pairs[
            np.random.randint(0, len(possible_pairs), size=n_samples)
        ]
        feature_dict = {
            feature_name: idx for idx, feature_name in enumerate(self.tree.tips())
        }
        augmented_matrix = []
        augmented_label = []
        for pair in selected_pairs:
            x1, x2 = self.matrix[:, pair[0]], self.matrix[:, pair[1]]
            _lambda = np.random.beta(alpha, alpha)
            n_leaves = int((1 - _lambda) * num_leaves)
            selected_index = set()
            mixed_x = x1.copy()

            while len(selected_index) < n_leaves:
                available_node = random.choice(all_nodes)
                leaf_idx = [feature_dict[leaf] for leaf in available_node.tips()]
                obs_idx = [
                    tip_to_obs_mapping[leaf.name] for leaf in available_node.tips()
                ]
                selected_index.update(obs_idx)

            selected_index = random.sample(list(selected_index), n_leaves)
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
            augment_label = _lambda * self.one_hot_label[pair[0]]
            +(1 - _lambda) * self.one_hot_label[pair[1]]

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
