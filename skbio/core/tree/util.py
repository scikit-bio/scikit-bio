from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from random import shuffle as shuffle_list

from future.builtins import zip


def shuffle(tree, n=None, names=None, shuffle_f=shuffle_list):
    """Yield trees with shuffled tip names

    Parameters
    ----------
    tree : TreeNode
        The tree to shuffle
    n : int, optional
        The number of tips to shuffle. If n is not `None`, n tips are randomly
        selected, and only those names will be shuffled.
    names : list, optional
        The specific tip names to shuffle. n and names cannot be specified at
        the same time.
    shuffle_f : func
        Shuffle method, this function must accept a list and modify inplace

    Notes
    -----
    Tip names are shuffled inplace.

    Returns
    -------
    GeneratorType
        Yielding TreeNode

    Raises
    ------
    ValueError
        If n is < 2
    ValueError
        If n and names are specified
    MissingNodeError
        If names are specified and a node name cannot be found

    Examples
    --------
    Alternate the names on two of the tips, 'a', and 'b', and do this 5 times.

    >>> from skbio.core.tree import TreeNode
    >>> tree = TreeNode.from_newick("((a,b),(c,d))")
    >>> rev = lambda items: items[::-1]
    >>> shuffler = shuffle(tree, names=['a', 'b'], shuffle_f=rev)
    >>> for idx, shuffled_tree in zip(range(5), shuffler):
    ...     print(shuffled_tree.to_newick())
    ((b,a),(c,d));
    ((a,b),(c,d));
    ((b,a),(c,d));
    ((a,b),(c,d));
    ((b,a),(c,d));

    """
    if n is not None and n < 2:
        raise ValueError("n must be None or >= 2")
    if n is not None and names is not None:
        raise ValueError("n and names cannot be specified at the sametime")

    tree.assign_ids()

    if names is None:
        all_tips = list(tree.tips())

        if n is None:
            n = len(all_tips)

        shuffle_f(all_tips)
        names = [tip.name for tip in all_tips[:n]]

    nodes = [tree.find(name) for name in names]

    while True:
        shuffle_f(names)
        for node, name in zip(nodes, names):
            node.name = name

        yield tree
