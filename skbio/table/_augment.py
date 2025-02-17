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


def biom_to_dataframe(biom_table):
    """Convert a biom table to pandas DataFrame

    Parameters
    ----------
    biom_table : biom.Table
        The biom table to convert to a pandas DataFrame.

    Returns
    -------
    pandas.DataFrame
    """
    return biom_table.to_dataframe()


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

    def __str__(self):
        return f"Augmentation(method={self.method})"

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
            The augmented label.
        """
        possible_pairs = self.get_all_possible_pairs()
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
                lambdas * self.label[indices1] + (1 - lambdas) * self.label[indices2]
            )
            augmented_matrix = np.concatenate([self.matrix, augmented_x], axis=1)
            augmented_label = np.concatenate([self.label, augmented_y])
        else:
            augmented_matrix = np.concatenate([self.matrix, augmented_x], axis=1)
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
            The augmented label.
        """

        if self.tree is None:
            raise ValueError("tree is required for phylomix augmentation")

        self.tree.bifurcate()
        leave_names = [tip.name for tip in self.tree.tips()]

        if set(tip_to_obs_mapping.keys()) != set(leave_names):
            raise ValueError("tip_to_obs_mapping must contain all tips in the tree")

        all_nodes = [node for node in self.tree.levelorder()]
        num_leaves = len(leave_names)

        possible_pairs = self.get_all_possible_pairs()
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

            augmented_matrix.append(mixed_x)
            augmented_label.append(pair[0])

        if self.label is not None:
            augmented_matrix = np.array(augmented_matrix).T
            augmented_label = np.array(augmented_label)
            augmented_matrix = np.concatenate([self.matrix, augmented_matrix], axis=1)
            augmented_label = np.concatenate([self.label, augmented_label])
        else:
            augmented_matrix = np.concatenate(
                [self.matrix, np.array(augmented_matrix).T], axis=1
            )
            augmented_label = None

        return augmented_matrix, augmented_label


if __name__ == "__main__":
    from biom import example_table
    from skbio.tree import TreeNode
    from collections import defaultdict

    # construnct 4 x 4 table
    data = np.arange(16).reshape(4, 4)
    sample_ids = ["S%d" % i for i in range(4)]
    observ_ids = ["O%d" % i for i in range(4)]
    sample_metadata = [
        {"environment": "A"},
        {"environment": "B"},
        {"environment": "A"},
        {"environment": "B"},
    ]
    observ_metadata = [
        {"phylogeny": "a"},
        {"phylogeny": "b"},
        {"phylogeny": "d"},
        {"phylogeny": "e"},
    ]
    table = Table(
        data, observ_ids, sample_ids, observ_metadata, sample_metadata, observ_metadata
    )

    # Create a simple tree
    tree = TreeNode.read(["((a,b)c,(d, e)f)root;"])
    tree.bifurcate()
    tree_tips = {tip.name for tip in tree.tips()}
    tips_to_obs_mapping = {}
    for idx, metadata in enumerate(table.metadata(axis="observation")):
        if metadata and "phylogeny" in metadata:
            phylogeny_label = metadata["phylogeny"]
            if phylogeny_label in tree_tips:
                tips_to_obs_mapping[phylogeny_label] = idx

    print("tips to Observation index mapping:")
    print(tips_to_obs_mapping)

    # Test phylomix method
    print(table)
    augmentation_phylomix = Augmentation(table, "phylomix", tree=tree)
    augmented_matrix_phylomix, augmented_label_phylomix = (
        augmentation_phylomix.phylomix(tips_to_obs_mapping, n_samples=3)
    )
    print(augmented_matrix_phylomix)
    print(augmented_label_phylomix)

    # Test mixup method
    augmentation_mixup = Augmentation(table, "mixup")
    augmented_matrix_mixup, augmented_label_mixup = augmentation_mixup.mixup(
        n_samples=3
    )
    print(augmented_matrix_mixup)
    print(augmented_label_mixup)
