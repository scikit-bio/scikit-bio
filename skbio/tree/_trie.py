# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import viewitems


class _CompressedNode(object):
    """Represents a node in the compressed trie

    Parameters
    ----------
    key : string
        the key attached to the node

    values : list of objects, optional
        the values attached to this node

    Attributes
    ----------
    values : list of objects
        the values attached to this node
    key : string
        the key attached to the node
    children : dict of {string: _CompressedNode}
        the children nodes below this node
    """

    def __init__(self, key, values=None):
        self.values = values or []
        self.key = key
        self.children = {}

    def __nonzero__(self):
        return (self.key != "" or len(self.values) > 0 or
                len(self.children.keys()) > 0)

    def __len__(self):
        """Returns the number of values attached to the node

        .. warning:: This method is recursive
        """
        return sum(len(n) for n in self.children.values()) + len(self.values)

    @property
    def size(self):
        """int with the number of nodes below the node

        .. warning:: This method is recursive
        """
        return sum(n.size for n in self.children.values()) + 1

    @property
    def prefix_map(self):
        """Dict with the prefix map

        Dictionary of {values: list of values} containing the prefix map
            of this node
        """
        mapping = {}

        if len(self.children) == 0:
            # we have a leaf
            mapping = {self.values[0]: self.values[1:]}
        else:
            # we are at an internal node
            for child in self.children.values():
                mapping.update(child.prefix_map)
            # get largest group
            n = -1
            key_largest = None
            for key, value in viewitems(mapping):
                if len(value) > n:
                    n = len(value)
                    key_largest = key
            # append this node's values
            mapping[key_largest].extend(self.values)

        return mapping

    def insert(self, key, value):
        """Inserts key with value in the node

        Parameters
        ----------
        key : string
            The string key attached to the value

        value : object
            Object to attach to the key
        """
        node_key_len = len(self.key)
        length = min(node_key_len, len(key))
        # Follow the key into the tree
        split_node = False
        index = 0
        while index < length and not split_node:
            split_node = key[index] != self.key[index]
            index += 1

        if split_node:
            # Index has been incremented after split_node was set to true,
            # decrement it to make it work
            index -= 1
            # We need to split up the node pointed by index
            # Get the key for the new node
            new_key_node = _CompressedNode(key[index:], [value])
            # Get a new node for the old key node
            old_key_node = _CompressedNode(self.key[index:], self.values)
            old_key_node.children = self.children
            self.children = {key[index]: new_key_node,
                             self.key[index]: old_key_node}
            self.key = self.key[:index]
            self.values = []
        elif index == len(self.key) and index == len(key):
            # The new key matches node key exactly
            self.values.append(value)
        elif index < node_key_len:
            # Key shorter than node key
            lower_node = _CompressedNode(self.key[index:], self.values)
            lower_node.children = self.children
            self.children = {self.key[index]: lower_node}
            self.key = key
            self.values = [value]
        else:
            # New key longer than current node key
            node = self.children.get(key[index])
            if node:
                # insert into next node
                node.insert(key[index:], value)
            else:
                # Create new node
                new_node = _CompressedNode(key[index:], [value])
                self.children[key[index]] = new_node

    def find(self, key):
        """Searches for key and returns values stored for the key.

        Parameters
        ----------
        key : string
            The key of the value to search for

        Returns
        -------
        object
            The value attached to the key
        """
        # key exhausted
        if len(key) == 0:
            return self.values

        # find matching part of key and node_key
        min_length = min(len(key), len(self.key))
        keys_diff = False
        index = 0
        while index < min_length and not keys_diff:
            keys_diff = key[index] != self.key[index]
            index += 1

        if keys_diff:
            return []
        elif index == len(key):
            # key and node_key match exactly
            return self.values
        else:
            node = self.children.get(key[index])
            if node:
                # descend to next node
                return node.find(key[index:])
        return []


class CompressedTrie(object):
    """ A compressed Trie for a list of (key, value) pairs

    Parameters
    ----------
    pair_list : list of tuples, optional
        List of (key, value) pairs to initialize the Trie

    Attributes
    ----------
    size
    prefix_map
    """

    def __init__(self, pair_list=None):
        self._root = _CompressedNode("")
        if pair_list:
            for key, value in pair_list:
                self.insert(key, value)

    def __nonzero__(self):
        return bool(self._root)

    def __len__(self):
        return len(self._root)

    @property
    def size(self):
        """int with the number of nodes in the Trie"""
        return self._root.size

    @property
    def prefix_map(self):
        """Dict with the prefix map

        Dictionary of {values: list of values} containing the prefix map
        """
        return self._root.prefix_map

    def insert(self, key, value):
        """Inserts key with value in Trie

        Parameters
        ----------
        key : string
            The string key attached to the value

        value : object
            Object to attach to the key
        """
        self._root.insert(key, value)

    def find(self, key):
        """Searches for key and returns values stored for the key.

        Parameters
        ----------
        key : string


        Returns
        -------
        object
            The value attached to the key
        """
        return self._root.find(key)


def fasta_to_pairlist(seqs):
    """Yields (key, value) pairs, useful for populating a Trie object

    Parameters
    ----------
    seqs : Iterable
        tuples of the form ``(label, seq)``

    Returns
    -------
    GeneratorType
        yields tuples of the form ``(seq, label)``
    """
    for label, seq in seqs:
        yield seq, label
