from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from random import shuffle as shuffle_list_inplace

from future.builtins import zip


def _shuffle_list(items):
    """Python's shuffle is inplace, wrap it to return"""
    shuffle_list_inplace(items)
    return items


def shuffle(tree, n=None, names=None, inplace=False, shuffle_f=_shuffle_list):
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
    inplace : bool
        If `True`, the names are shuffled on the tree inplace. If `False`, a
        copy of the tree is made on each iteration.
    shuffle_f : func
        Shuffle method, this function must accept a list and return a list.

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

    original_tree = tree
    tree.assign_ids()
    all_tips = shuffle_f(list(tree.tips()))

    if names is None:
        if n is None:
            n = len(all_tips)

        ids_to_shuffle = [tip.id for tip in all_tips[:n]]
        names_to_shuffle = [tip.name for tip in all_tips[:n]]
    else:
        ids_to_shuffle = [tip.id for tip in all_tips if tip.name in set(names)]

        if len(ids_to_shuffle) != len(names):
            raise ValueError("Cannot find all the names to shuffle!")

        names_to_shuffle = [tree.find_by_id(i).name for i in ids_to_shuffle]

    while True:
        if not inplace:
            tree = original_tree.copy()

        names_to_shuffle = shuffle_f(names_to_shuffle)
        for id_, name in zip(ids_to_shuffle, names_to_shuffle):
            tree.find_by_id(id_).name = name

        yield tree
