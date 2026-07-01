# encoding: utf-8
# cython: profile=False, boundscheck=False, wraparound=False
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Derived from improved-octo-waddle (https://github.com/biocore/improved-octo-waddle)
# originally authored by Daniel McDonald, distributed under the Modified BSD License.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._bp cimport BPTree
import numpy as np
import pandas as pd
import json
cimport numpy as cnp
cimport cython

cnp.import_array()


cdef inline cnp.double_t length_from_edge(unicode token):
    cdef:
        Py_ssize_t split_idx

    # 0.12345{0123} -> 0.12345
    # OR 0.12345[0123] -> 0.12345
    split_idx_curly = token.find('{')
    split_idx_square = token.find('[')
    split_idx = max(split_idx_curly, split_idx_square)
    if split_idx == -1:
        return np.double(token)
    else:
        return np.double(token[:split_idx])


cdef inline cnp.int32_t number_from_edge(unicode token):
    cdef:
        Py_ssize_t split_idx
        Py_ssize_t end

    # 0.12345{0123} -> 0123
    # OR 0.12345[0123] -> 0123
    split_idx_curly = token.find('{')
    split_idx_square = token.find('[')
    split_idx = max(split_idx_curly, split_idx_square)
    if split_idx == -1:
        return 0
    else:
        end = len(token)
        return np.int32(token[split_idx + 1:end - 1])


cdef void _set_node_metadata(cnp.uint32_t ptr, unicode token,
                             cnp.ndarray[object, ndim=1] names,
                             cnp.ndarray[cnp.double_t, ndim=1] lengths,
                             cnp.ndarray[cnp.int32_t, ndim=1] edges):
    """Inplace update of names and lengths given token details"""
    cdef:
        cnp.double_t length
        cnp.int32_t edge
        Py_ssize_t split_idx, i, end
        unicode name, token_parsed

    name = None
    length = 0.0
    edge = 0

    if token[0] == u':':
        token_parsed = token[1:]
        length = length_from_edge(token_parsed)
        edge = number_from_edge(token_parsed)
    elif u':' in token:
        split_idx = token.rfind(':')
        name = token[:split_idx]
        token_parsed = token[split_idx + 1:]
        length = length_from_edge(token_parsed)
        edge = number_from_edge(token_parsed)
        name = name.strip("'").strip()
    elif u'{' in token or u'[' in token:
        # strip as " {123}" is valid?
        token = token.strip()
        end = len(token)
        edge = np.int32(token.strip()[1:end - 1])
    else:
        name = token.replace("'", "").replace('"', "").strip()

    names[ptr] = name
    lengths[ptr] = length
    edges[ptr] = edge


def write_newick(BPTree tree, object output, bint include_edge):
    """Write a BPTree as a Newick string.

    Parameters
    ----------
    tree : BPTree
        The tree to serialize.
    output : file-like object
        An open, writable handle (anything with a ``write`` method) that the
        Newick string is written to.
    include_edge : bool
        If ``True``, write each node's edge number alongside its branch length
        using the ``:length{edge}`` convention. If ``False``, write the branch
        length only.

    Notes
    -----
    The string is written in place to ``output``; the function returns nothing.
    Node names containing any of the characters ``;``, ``,``, ``(``, ``)``,
    ``:``, or ``_`` are wrapped in single quotes.

    """
    cdef:
        list name_stack
        list edge_stack
        list length_stack
        list open_paren_stack
        object name
        cnp.float64_t length
        Py_ssize_t idx
        cnp.uint8_t v
        Py_ssize_t root_close

    length_stack = []
    name_stack = []
    edge_stack = []
    open_paren_stack = []
    root_close = tree.close(0)

    for idx, v in enumerate(tree.data):
        if v:
            if not tree.isleaf(idx):
                output.write('(')
            name_stack.append(tree.name(idx))
            length_stack.append(tree.length(idx))
            edge_stack.append(tree.edge(idx))
            open_paren_stack.append(idx)
        else:
            name = name_stack.pop()
            length = length_stack.pop()
            edge = edge_stack.pop()

            if name is not None:
                # if we have magical characters, make sure we quote
                if set(name) & {';', ',', '(', ')', ':', '_'}:
                    output.write("'%s'" % name)
                else:
                    output.write(name)

            if include_edge:
                output.write(':%f{%d}' % (length, edge))
            else:
                output.write(':%f' % length)

            if tree.nsibling(open_paren_stack.pop()) == 0:
                if idx != root_close:
                    output.write(')')
            else:
                output.write(',')

    output.write(';')


cpdef parse_newick(unicode data):
    """Parse a Newick string into a BPTree.

    Parameters
    ----------
    data : str
        A Newick-formatted string, terminated with a semicolon (``;``). Branch
        lengths and edge numbers (delimited by ``{}`` or ``[]``) are parsed if
        present.

    Returns
    -------
    BPTree
        The tree topology, with names, branch lengths, and edge numbers
        attached.

    Raises
    ------
    ValueError
        If ``data`` is not terminated with a semicolon, or if it describes a
        tree with only a single node.

    """
    cdef:
        cnp.uint32_t ptr, open_ptr
        Py_ssize_t token_ptr, tmp, lag, datalen
        BPTree topology
        unicode token, last_token
        cnp.ndarray[object, ndim=1] names
        cnp.ndarray[cnp.double_t, ndim=1] lengths
        cnp.ndarray[cnp.int32_t, ndim=1] edges

    data = data.strip()
    if not data.endswith(';'):
        raise ValueError("Newick does not appear terminated with a semicolon")

    datalen = len(data)
    topology = _newick_to_bp(data)

    if len(topology.data) <= 2:
        raise ValueError("Only trees with more than 1 node supported")

    names = np.full(len(topology.data), None, dtype=object)
    lengths = np.zeros(len(topology.data), dtype=np.double)
    edges = np.full(len(topology.data), 0, dtype=np.int32)

    ptr = 0
    token_ptr = _ctoken(data, datalen, 0)
    token = data[0:token_ptr]
    last_token = None

    # lag reflects the scenario where ((x))y, where the label y gets may end
    # up being associated with an earlier unnamed vertex. lag represents the
    # offset between the topology pointer and the token pointer effectively.
    lag = 0
    while token != ';':
        if token == '(':
            # an open parenthesis never has metadata associated with it
            ptr += 1

        if (token == ')' or token == ',') and last_token == ')':
            # determine if there are unnamed/unlengthed nodes
            lag += 1

        elif token not in '(),:;':
            ptr += lag
            lag = 0

            open_ptr = topology.open(ptr)
            _set_node_metadata(open_ptr, token, names, lengths, edges)

            if topology.isleaf(ptr):
                ptr += 2
            else:
                ptr += 1

        last_token = token
        tmp = _ctoken(data, datalen, token_ptr)
        token = data[token_ptr:tmp]
        token_ptr = tmp

    topology.set_names(names)
    topology.set_lengths(lengths)
    topology.set_edges(edges)

    return topology


cdef object _newick_to_bp(unicode data):
    """Convert newick to balanced parentheses"""
    cdef:
        Py_ssize_t i, topology_ptr, single_descendent
        Py_UCS4 c, last_c
        cnp.ndarray[cnp.uint8_t, ndim=1] topology

    potential_single_descendant = False

    topology = np.empty(len(data), dtype=np.uint8)
    topology_ptr = 0
    last_c = u'x'
    in_quote = False

    for i in range(len(data)):
        c = data[i]
        if c == u"'":
            in_quote = not in_quote
        else:
            if in_quote:
                continue
            elif c == u'(':
                # opening of a node
                topology[topology_ptr] = 1
                topology_ptr += 1
                last_c = c
                potential_single_descendant = True
            elif c == u')':
                # closing of a node
                if potential_single_descendant or last_c == u',':
                    # we have a single descendant or a last child (i.e., ",)")
                    topology[topology_ptr] = 1
                    topology[topology_ptr + 1] = 0
                    topology[topology_ptr + 2] = 0
                    topology_ptr += 3
                    potential_single_descendant = False
                else:
                    topology[topology_ptr] = 0
                    topology_ptr += 1
                last_c = c
            elif c == u',':
                if last_c != u')':
                    # we have a new tip
                    topology[topology_ptr] = 1
                    topology[topology_ptr + 1] = 0
                    topology_ptr += 2
                potential_single_descendant = False
                last_c = c
            else:
                # ignore non-structure
                pass

    return BPTree(topology[:topology_ptr])


cdef inline int _ccheck(Py_UCS4 c):
    """structure check"""
    cdef:
        Py_ssize_t i

    if c == u'(':
        return 1
    elif c == u')':
        return 1
    elif c == u',':
        return 1
    elif c == u';':
        return 1
    else:
        return 0


cdef inline int _is_quote(Py_UCS4 c):
    if c == u'"':
        return 1
    elif c == u"'":
        return 1
    else:
        return 0


cdef inline Py_ssize_t _ctoken(unicode data, Py_ssize_t datalen, Py_ssize_t start):
    cdef:
        Py_ssize_t idx, in_quote = 0
        Py_UCS4 c

    if start == (datalen - 1):
        return start + 1

    for idx in range(start, datalen):
        c = data[idx]

        if in_quote:
            if _is_quote(c):
                in_quote = 0
            continue
        else:
            if _is_quote(c):
                in_quote = 1
                continue

        if _ccheck(c):
            if idx == start:
                return idx + 1
            else:
                return idx

    return idx + 1


def parse_jplace(object data):
    """Parse a jplace string into a placement table and a tree.

    Parameters
    ----------
    data : str
        A jplace-formatted JSON string, as produced by phylogenetic placement
        tools.

    Returns
    -------
    pandas.DataFrame
        The placements, with one row per (fragment, placement) pair. Columns
        are ``fragment`` followed by the ``fields`` declared in the jplace
        document. Rows are restricted to edges present in the tree.
    BPTree
        The reference tree parsed from the jplace ``tree`` entry.

    Notes
    -----
    Implementation-specific caveats:

    1. Multiplicities are not supported. Placements are required to have an
       ``"n"`` entry, and any ``"nm"`` entry is ignored.
    2. Matsen et al. [1]_ define ``[]`` for edge labels and ``{}`` for edge
       numbers. Either ``[]`` or ``{}`` is supported, but not edges with both.
       Edge labels, if specified, are required to be integers.

    References
    ----------
    .. [1] Matsen FA, Hoffman NG, Gallagher A, Stamatakis A. (2012). A format
       for phylogenetic placements. PLoS ONE 7(2): e31009.

    """
    cdef:
        dict as_json
        list fields, placements, fragments, p, placement_data,
        list placement_inner_data, pquery, entry
        unicode frag, newick
        Py_ssize_t placement_idx, placement_inner_idx, fragment_idx,
        Py_ssize_t n_fragments
        BPTree tree
        object df
        set edges

    as_json = json.loads(data)
    newick = as_json['tree']
    placement_data = as_json['placements']

    fields = as_json['fields']
    fields = ['fragment', ] + fields

    placements = []
    for placement_idx in range(len(placement_data)):
        placement = placement_data[placement_idx]

        placement_inner_data = placement['p']

        if 'n' not in placement:
            raise KeyError("jplace parsing limited to entries with 'n' keys")

        fragments = placement['n']
        n_fragments = len(fragments)

        for placement_inner_idx in range(len(placement_inner_data)):
            pquery = placement_inner_data[placement_inner_idx]

            for fragment_idx in range(n_fragments):
                frag = fragments[fragment_idx]
                entry = [frag, ] + pquery
                placements.append(entry)

    tree = parse_newick(newick)
    edges = {tree.edge(i) for i, v in enumerate(tree.data) if v}
    df = pd.DataFrame(placements, columns=fields)
    df = df[df['edge_num'].isin(edges)]

    return df, tree
