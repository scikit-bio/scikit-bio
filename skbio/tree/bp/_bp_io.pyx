# encoding: utf-8
# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
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
    split_idx = token.find('{')
    if split_idx == -1:
        return np.double(token)
    else:
        return np.double(token[:split_idx])


cdef inline cnp.int32_t number_from_edge(unicode token):
    cdef:
        Py_ssize_t split_idx
        Py_ssize_t end

    # 0.12345{0123} -> 0123
    split_idx = token.find('{')
    if split_idx == -1:
        return 0
    end = token.find('}')
    if end == -1:
        end = len(token)
    return np.int32(token[split_idx + 1:end])


cdef inline unicode _unquote_name(unicode token, bint convert_underscores):
    """Recover a node name from a raw token.

    A single-quoted label is de-quoted and its doubled apostrophes (``''``) are
    collapsed to one (``'``), matching the Newick spec and the TreeNode reader.
    An unquoted label optionally has underscores converted to spaces (the
    classic Newick convention); quoted content is always taken literally.
    """
    token = token.strip()
    if len(token) >= 2 and token[0] == u"'" and token[len(token) - 1] == u"'":
        return token[1:len(token) - 1].replace(u"''", u"'")
    if convert_underscores:
        return token.replace(u'_', u' ')
    return token


cdef inline Py_ssize_t _find_length_colon(unicode token):
    """Index of the branch-length ``:`` (the first ``:`` outside single quotes)."""
    cdef:
        Py_ssize_t i, n = len(token)
        bint in_quote = False
        Py_UCS4 c

    for i in range(n):
        c = token[i]
        if c == u"'":
            in_quote = not in_quote
        elif c == u':' and not in_quote:
            return i
    return -1


cdef unicode _strip_comments(unicode data):
    """Remove Newick ``[comment]`` regions.

    Comments may nest and are skipped entirely (their content is discarded, as
    in the TreeNode reader). A ``[`` inside a single-quoted label is literal.
    """
    cdef:
        Py_ssize_t i, n = len(data), seg_start = 0
        int comment_depth = 0
        bint in_quote = False
        Py_UCS4 c
        list out = []

    for i in range(n):
        c = data[i]
        if in_quote:
            if c == u"'":
                in_quote = False
            continue
        if comment_depth > 0:
            if c == u'[':
                comment_depth += 1
            elif c == u']':
                comment_depth -= 1
                if comment_depth == 0:
                    seg_start = i + 1
            continue
        if c == u"'":
            in_quote = True
        elif c == u'[':
            out.append(data[seg_start:i])
            comment_depth += 1

    out.append(data[seg_start:n])
    return u''.join(out)


cdef void _set_node_metadata(cnp.uint32_t ptr, unicode token,
                             cnp.ndarray[object, ndim=1] names,
                             cnp.ndarray[cnp.double_t, ndim=1] lengths,
                             cnp.ndarray[cnp.int32_t, ndim=1] edges,
                             bint convert_underscores):
    """Inplace update of names, lengths, and edges given token details."""
    cdef:
        cnp.double_t length
        cnp.int32_t edge
        Py_ssize_t split_idx, curly
        unicode name, token_parsed

    name = None
    length = 0.0
    edge = 0

    if len(token) > 0 and token[0] == u':':
        # length (and optional edge number) only, e.g. ":0.5" or ":0.5{3}"
        token_parsed = token[1:]
        length = length_from_edge(token_parsed)
        edge = number_from_edge(token_parsed)
    else:
        split_idx = _find_length_colon(token)
        if split_idx != -1:
            # name:length[{edge}]
            name = _unquote_name(token[:split_idx], convert_underscores)
            token_parsed = token[split_idx + 1:]
            length = length_from_edge(token_parsed)
            edge = number_from_edge(token_parsed)
        elif u'{' in token:
            # an edge number with no branch length, e.g. "{3}" or "name{3}"
            curly = token.find(u'{')
            if curly > 0:
                name = _unquote_name(token[:curly], convert_underscores)
            edge = number_from_edge(token[curly:])
        else:
            name = _unquote_name(token, convert_underscores)

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
    Node names containing whitespace or any of the characters
    ``' _ ; , : ( ) [ ]`` are wrapped in single quotes, and embedded apostrophes
    are doubled (``'`` -> ``''``), so the output re-parses losslessly.

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
            if not tree.is_tip(idx):
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
                # Quote names with structural characters or whitespace, and
                # double embedded apostrophes, so the output re-parses cleanly.
                chars = set(name)
                if (chars & {"'", "_", ";", ",", ":", "(", ")", "[", "]"}
                        or any(ch.isspace() for ch in chars)):
                    output.write("'%s'" % name.replace("'", "''"))
                else:
                    output.write(name)

            if include_edge:
                output.write(':%s{%d}' % (length, edge))
            else:
                output.write(':%s' % length)

            if tree.next_sibling(open_paren_stack.pop()) == 0:
                if idx != root_close:
                    output.write(')')
            else:
                output.write(',')

    output.write(';\n')


cpdef parse_newick(unicode data, bint convert_underscores=True):
    """Parse a Newick string into a BPTree.

    Parameters
    ----------
    data : str
        A Newick-formatted string, terminated with a semicolon (``;``). Branch
        lengths and ``{}``-delimited edge numbers are parsed if present, and
        ``[]`` comments are ignored.
    convert_underscores : bool, optional
        If ``True`` (default), underscores in unquoted labels are converted to
        spaces, matching the TreeNode Newick reader. Quoted labels are always
        taken literally.

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

    data = _strip_comments(data.strip()).strip()
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
            _set_node_metadata(open_ptr, token, names, lengths, edges,
                               convert_underscores)

            if topology.is_tip(ptr):
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
        The reference tree parsed from the jplace ``tree`` entry, carrying its
        ``{}`` edge numbers.

    Raises
    ------
    ValueError
        If a required member (``tree``, ``placements``, ``fields``) is missing,
        or if ``fields`` does not include ``edge_num``.
    KeyError
        If a placement record lacks the required ``n`` member.

    Notes
    -----
    Implementation-specific caveats:

    1. Multiplicities are not supported. Placements are required to have an
       ``"n"`` entry, and any ``"nm"`` entry is ignored.
    2. Matsen et al. [1]_ allow ``[]`` or ``{}`` to delimit edge numbers. Here
       edge numbers use ``{}`` only; ``[]`` is treated as a standard Newick
       comment and ignored, so legacy ``[]``-delimited edge numbers are not
       recognized.

    References
    ----------
    .. [1] Matsen FA, Hoffman NG, Gallagher A, Stamatakis A. (2012). A format
       for phylogenetic placements. PLoS ONE 7(2): e31009.

    """
    cdef:
        dict as_json
        list fields, columns, placements, fragments, placement_data
        list placement_inner_data, pquery, entry
        unicode frag, newick
        Py_ssize_t placement_idx, placement_inner_idx, fragment_idx
        Py_ssize_t n_fragments
        BPTree tree
        object df
        set edges

    as_json = json.loads(data)

    for key in ('tree', 'placements', 'fields'):
        if key not in as_json:
            raise ValueError(
                "jplace document is missing the required '%s' member" % key
            )

    newick = as_json['tree']
    placement_data = as_json['placements']
    fields = list(as_json['fields'])

    if 'edge_num' not in fields:
        raise ValueError(
            "jplace 'fields' must include 'edge_num' to map placements onto "
            "the reference tree"
        )

    columns = ['fragment', ] + fields

    placements = []
    for placement_idx in range(len(placement_data)):
        placement = placement_data[placement_idx]

        if 'n' not in placement:
            raise KeyError("jplace parsing limited to entries with 'n' keys")

        placement_inner_data = placement['p']
        fragments = placement['n']
        n_fragments = len(fragments)

        for placement_inner_idx in range(len(placement_inner_data)):
            pquery = placement_inner_data[placement_inner_idx]

            for fragment_idx in range(n_fragments):
                frag = fragments[fragment_idx]
                entry = [frag, ] + pquery
                placements.append(entry)

    # Preserve underscores in reference-tree labels: jplace taxa are typically
    # accessions where '_' is significant, not a stand-in for a space.
    tree = parse_newick(newick, convert_underscores=False)
    edges = {tree.edge(i) for i, v in enumerate(tree.data) if v}
    df = pd.DataFrame(placements, columns=columns)
    df = df[df['edge_num'].isin(edges)]

    return df, tree


def _json_default(object obj):
    """Fallback JSON encoder that unwraps NumPy scalars to Python types."""
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(
        "Object of type %s is not JSON serializable" % type(obj).__name__
    )


def write_jplace(BPTree tree, object output, object fields=None,
                 object version=None, object metadata=None):
    """Write a reference tree as a jplace document.

    A ``BPTree`` holds a reference tree but no placements, so the document is
    written with an empty ``placements`` list; ``fields``, ``version``, and
    ``metadata`` fill the remaining required jplace members. This mirrors the
    ``jplace`` read path, which returns only the reference tree.

    Parameters
    ----------
    tree : BPTree
        The reference tree, serialized as Newick with ``{}`` edge numbers via
        :func:`write_newick`. Existing edge numbers are written as-is; a tree
        without edge numbers emits ``{0}`` on every edge.
    output : file-like object
        An open, writable handle that the JSON document is written to.
    fields : list of str, optional
        The placement field names to record. Defaults to ``["edge_num"]``.
    version : int, optional
        The jplace format version to record. Defaults to ``3``.
    metadata : dict, optional
        Free-form metadata to record. Defaults to ``{}``.

    """
    cdef:
        unicode tree_str
        object buf, document

    from io import StringIO

    buf = StringIO()
    write_newick(tree, buf, True)
    tree_str = buf.getvalue().strip()

    document = {
        "tree": tree_str,
        "placements": [],
        "fields": list(fields) if fields is not None else ["edge_num"],
        "version": version if version is not None else 3,
        "metadata": metadata if metadata is not None else {},
    }

    json.dump(document, output, default=_json_default)
