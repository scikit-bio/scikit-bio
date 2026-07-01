r"""Balanced parentheses (BP) tree backend (:mod:`skbio.tree.bp`)
================================================================

.. currentmodule:: skbio.tree.bp

A succinct balanced-parentheses representation of trees, ported from
improved-octo-waddle (https://github.com/biocore/improved-octo-waddle). The
public objects are re-exported from :mod:`skbio.tree`.

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._bp import BPTree
from ._bp_io import parse_newick, write_newick, parse_jplace
from ._bp_conv import to_skbio_treenode, from_skbio_treenode, to_skbio_treearray
from ._bp_insert import insert_fully_resolved

__all__ = [
    "BPTree",
    "parse_newick",
    "write_newick",
    "parse_jplace",
    "to_skbio_treenode",
    "from_skbio_treenode",
    "to_skbio_treearray",
    "insert_fully_resolved",
]
