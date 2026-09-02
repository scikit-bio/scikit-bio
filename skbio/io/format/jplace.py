r"""jplace format (:mod:`skbio.io.format.jplace`)
==============================================

.. currentmodule:: skbio.io.format.jplace

The jplace format (``jplace``) stores the results of a phylogenetic placement:
the assignment of query sequences (e.g., metagenomic reads) onto the edges of a
fixed reference tree [1]_. It is a JSON document produced by placement tools
such as pplacer, EPA-ng, and APPLES.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.tree.BPTree`                                       |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
A jplace file is a JSON object with the following top-level members [1]_:

- ``tree``: a Newick string in which each edge is annotated with an integer
  edge number, delimited by ``{}``. (``[]`` is treated as a standard Newick
  comment and ignored.)
- ``placements``: a list of placement records. Each record has a ``p`` member
  (a list of placement rows, each row a list matching ``fields``) and an ``n``
  member (a list of query/fragment names).
- ``fields``: the column names for the rows in ``p`` (must include
  ``edge_num``).
- ``version``: the jplace format version.
- ``metadata``: free-form metadata (e.g., the invocation command line).

Reading a jplace file yields the reference tree as a
:class:`~skbio.tree.BPTree`, the succinct balanced-parentheses representation
that carries the ``{}`` edge numbers and scales to the large reference trees
typical of placement studies. The placement records themselves are **dropped**;
use :func:`skbio.tree.bp.parse_jplace` to obtain the placements table alongside
the tree.

Writing a ``BPTree`` produces a jplace document containing the reference tree
(with its ``{}`` edge numbers) and an empty ``placements`` list -- a ``BPTree``
holds a reference tree but no placements. The tree parse and serialization run
in scikit-bio's compiled backend, so throughput is preserved for large files.

Caveats
~~~~~~~
- Multiplicities are not supported by :func:`skbio.tree.bp.parse_jplace`; any
  ``"nm"`` entry is ignored.
- Edge numbers use ``{}`` only. ``[]`` is treated as a standard Newick comment.
- A jplace document written from an edge-number-less ``BPTree`` repeats ``{0}``
  on every edge; the round trip is only meaningful when the tree carries edge
  numbers (as a tree read from jplace does).

Examples
--------
Read the reference tree of a jplace file into a :class:`~skbio.tree.BPTree`
(placements are dropped):

.. code-block:: python

   from skbio.tree import BPTree
   tree = BPTree.read("placement.jplace", format="jplace")

Write a :class:`~skbio.tree.BPTree` back out as a (placement-less) jplace file:

.. code-block:: python

   tree.write("reference.jplace", format="jplace")

References
----------
.. [1] Matsen FA, Hoffman NG, Gallagher A, Stamatakis A. (2012). A format for
   phylogenetic placements. PLoS ONE 7(2): e31009.


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import json

from skbio.io import create_format, JplaceFormatError
from skbio.tree import BPTree

jplace = create_format("jplace")

_REQUIRED_KEYS = ("tree", "placements", "fields", "version")


@jplace.sniffer()
def _jplace_sniffer(fh):
    # A jplace record is a JSON object. ``json.load`` fails fast on anything
    # that is not JSON (e.g., a Newick record, which begins with ``(``), so the
    # cost of sniffing a non-jplace file is a cheap parse error. We then require
    # all of the mandatory top-level members to be present, which is
    # unambiguous relative to every other scikit-bio format.
    try:
        data = json.load(fh)
    except (ValueError, UnicodeDecodeError):
        return False, {}
    if not isinstance(data, dict):
        return False, {}
    if all(key in data for key in _REQUIRED_KEYS):
        return True, {}
    return False, {}


def _read_tree(fh):
    """Parse a jplace document and return its reference tree (drop placements)."""
    from skbio.tree.bp._bp_io import parse_jplace

    # The compiled parser scans the whole document at once; materialize the
    # handle first. It raises ValueError/KeyError on malformed input; re-raise
    # as the format's error type per the registry convention.
    try:
        _, tree = parse_jplace(fh.read())
    except (ValueError, KeyError) as e:
        raise JplaceFormatError("Could not parse file as jplace: %s" % e) from e
    return tree


@jplace.reader(BPTree)
def _jplace_to_bp(fh, cls=None):
    return _read_tree(fh)


@jplace.writer(BPTree)
def _bp_to_jplace(obj, fh):
    from skbio.tree.bp._bp_io import write_jplace

    write_jplace(obj, fh)
