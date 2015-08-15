from __future__ import division
import numpy as np
from ._unifrac import _validate


def unifrac(branch_lengths, i, j):
    """Calculates unifrac(i,j) from branch lengths and cols i and j of m.

    This is the original, unweighted UniFrac metric.

    Parameters
    ----------
    branch_lengths: np.array
         branch_lengths should be row vector, same length as # nodes in tree.
    i : np.array
        a slice of states from m, same length as # nodes in tree.
        Slicing m (e.g. m[:,i]) returns a vector in the right format; note
        that it should be a row vector (the default), not a column vector.
    j : np.array
        a slice of states from m, same length as # nodes in tree.
        Slicing m (e.g. m[:,j]) returns a vector in the right format; note
        that it should be a row vector (the default), not a column vector.

    Returns
    -------
    float
        Unweighted unifrac metric

    Notes
    -----
    This is cogent.maths.unifrac.fast_tree.unifrac, but there are
    other metrics that can (should?) be ported like:
     - unnormalized_unifrac
     - G
     - unnormalized_G
    """

    _or, _and = np.logical_or(i, j), np.logical_and(i, j)
    return 1 - ((branch_lengths * _and).sum() / (branch_lengths * _or).sum())


def w_unifrac(branch_lengths, i, j):
    """Calculates weighted unifrac(i,j) from branch lengths and cols i,j of m.
    """
    return (branch_lengths * abs((i / i.sum()) - (j / j.sum()))).sum()


def _skbio_counts_to_envs(otu_ids, *counts):
    """This is a place holder method to help get the API hooked up

    This method should be unnecessary once the implementation is optimized to
    the interface.
    """
    envs = {}
    n_counts = len(counts)

    for packed in zip(otu_ids, *counts):
        # NOTE: this is ducttape to fit the API
        otu_id = packed[0]

        counts = {}
        for env in range(1, n_counts + 1):
            if packed[env]:
                counts[env] = packed[env]

        if counts:
            envs[otu_id] = counts

    return envs


def unweighted_unifrac(u_counts, v_counts, otu_ids, tree,
                       suppress_validation=False):
    """fit to unifrac API"""
    if not suppress_validation:
        _validate(u_counts, v_counts, otu_ids, tree)

    u_sum = sum(u_counts)
    v_sum = sum(v_counts)

    if not u_sum or not v_sum:
        if u_sum + v_sum:
            # u or v counts are all zeros
            return 1.0
        else:
            # u and v are zero
            return 0.0

    envs = _skbio_counts_to_envs(otu_ids, u_counts, v_counts)

    return fast_unifrac(tree, envs)


def weighted_unifrac(u_counts, v_counts, otu_ids, tree, normalized=False,
                     suppress_validation=False):
    u_sum = sum(u_counts)
    v_sum = sum(v_counts)

    if not u_sum or not v_sum:
        if u_sum + v_sum:
            # u or v counts are all zeros
            if normalized:
                return 1.0
        else:
            # u and v are zero
            return 0.0
    if not suppress_validation:
        _validate(u_counts, v_counts, otu_ids, tree)
    envs = _skbio_counts_to_envs(otu_ids, u_counts, v_counts)

    return fast_unifrac(tree, envs, weighted=True)


def fast_unifrac(t, envs, weighted=False, metric=unifrac):
    # weighted_unifrac_f=_weighted_unifrac,make_subtree=True):
    """
    Run fast unifrac.

    Parameters
    ----------
    t: skbio.TreeNode
        phylogenetic tree relating the sequences.
    envs: dict
        dict of {sequence:{env:count}} showing environmental abundance.
    weighted: bool
        if True, performs the weighted UniFrac procedure.
    metric: function
        distance metric to use.  currently you must use unifrac only
        if weighted=True.
        see fast_tree.py for metrics (e.g.: G, unnormalized_G, unifrac, etc.)

    Returns
    -------
    u : np.array
        Unifrac distance matrix
    """
    (envs, count_array,
     unique_envs, env_to_index,
     node_to_index, env_names,
     branch_lengths, nodes, t) = _fast_unifrac_setup(t, envs)

    bound_indices = bind_to_array(nodes, count_array)
    # weighted unifrac
    if weighted:
        tip_indices = [n._leaf_index for n in t.tips()]
        sum_descendants(bound_indices)
        tip_ds = branch_lengths.copy()[:, np.newaxis]
        bindings = bind_to_parent_array(t, tip_ds)
        tip_distances(tip_ds, bindings, tip_indices)
        u = w_unifrac(branch_lengths, count_array[:, 0], count_array[:, 1])
    # unweighted unifrac
    else:
        bool_descendants(bound_indices)
        u = unifrac(branch_lengths, count_array[:, 0], count_array[:, 1])

    return u


def _fast_unifrac_setup(t, envs, make_subtree=True):
    """Setup shared by fast_unifrac and by significance tests."""
    if make_subtree:
        t2 = t.copy()
        wanted = set(envs.keys())

        def delete_test(node):
            if node.is_tip() and node.name not in wanted:
                return True
            return False
        t2.remove_deleted(delete_test)
        t2.prune()
        t = t2

    # index tree
    node_index, nodes = index_tree(t)
    # get good nodes, defined as those that are in the env file.
    good_nodes = dict([(i.name, envs[i.name])
                       for i in t.tips() if i.name in envs])
    envs = good_nodes
    (count_array, unique_envs,
     env_to_index, node_to_index) = index_envs(envs, node_index)
    env_names = sorted(unique_envs)
    # Note: envs get sorted at the step above
    branch_lengths = get_branch_lengths(node_index)
    if not envs:
        raise ValueError, ''.join(["No valid samples/environments found.",
                                   "Check whether tree tips match ",
                                   "otus/taxa present ",
                                   "in samples/environments"])
    return (envs, count_array, unique_envs, env_to_index,
            node_to_index, env_names, branch_lengths, nodes, t)


def index_tree(t):
    """

    Indexes nodes in-place as n._leaf_index.

    Algorithm is as follows:
    for each node in post-order traversal over tree:
        if the node has children:
            set an index on each child
            for each child with children:
                add the child and its start and end tips to the result

    Parameters
    ----------
    t : skbio.TreeNode
        Phylogenetic tree

    Returns
    -------
    tuple
        Contains {node_id:node}, [node_id,first_child,last_child]
    """
    id_index = {}    # needs to be dict, not list, b/c adding out of order
    child_index = []
    curr_index = 0
    for n in t.traverse(self_before=False, self_after=True):
        for c in n.children:
            c._leaf_index = curr_index
            id_index[curr_index] = c
            curr_index += 1
            if c:    # c has children itself, so need to add to result
                child_index.append((c._leaf_index, c.children[0]._leaf_index,
                                    c.children[-1]._leaf_index))
    # handle root, which should be t itself
    t._leaf_index = curr_index
    id_index[curr_index] = t
    # only want to add to the child_index if t has children...
    if t.children:
        child_index.append((t._leaf_index, t.children[0]._leaf_index,
                            t.children[-1]._leaf_index))
    return id_index, child_index


def index_envs(env_counts, tree_index, array_constructor=int):
    """
    Calculates the taxon counts in each env (what is an env?)

    Parameters
    ----------
    env_counts : np.array
       Should be the output of count_envs(lines).
    tree_index : np.array
       Should be the id_index of index_tree(t).
    array_constructor: function
        int by default (may need to change to float later
        to handle microarray data).

    Returns
    -------
    np.array
        Array of taxon x env with counts of the taxon in each env.
    """
    num_nodes = len(tree_index)
    unique_envs, num_envs = get_unique_envs(env_counts)
    env_to_index = dict([(e, i) for i, e in enumerate(unique_envs)])
    result = np.zeros((num_nodes, num_envs), array_constructor)
    # figure out taxon label to index map
    node_to_index = {}
    for i, node in tree_index.items():
        if node.name is not None:
            node_to_index[node.name] = i
    # walk over env_counts, adding correct slots in array
    for name in env_counts:
        curr_row_index = node_to_index[name]
        for env, count in env_counts[name].items():
            result[curr_row_index, env_to_index[env]] = count
    # return all the data structures we created; will be useful for other tasks
    return result, unique_envs, env_to_index, node_to_index


def get_branch_lengths(tree_index):
    """
    Parameters
    ----------
    tree_index: dict


    Returns
    -------
    np.array
        Array of branch lengths, in tree index order.
    """
    result = np.zeros(len(tree_index), float)
    for i, node in tree_index.items():
        try:
            if node.length is not None:
                result[i] = node.length
        except AttributeError:
            pass
    return result


def _branch_correct(tip_distances, i, j):
    """Calculates weighted unifrac branch length correction.

    tip_distances  must be 0 except for tips.
    """
    result = tip_distances.ravel()*((i / i.sum())+(j / j.sum()))
    return result.sum()


def get_unique_envs(envs):
    """extract all unique envs from envs dict"""
    result = set()
    for v in envs.values():
        result.update(v.keys())
    # sort envs for convenience in testing and display
    return sorted(result), len(result)


def bind_to_array(tree_index, a):
    """Binds tree_index to array a, returning result in list.

    Takes as input list of (node, first_child, last_child)
    returns list of (node_row, child_rows) such that node_row points to the
    row of a that corresponds to the current node, and child_rows points to the
    row or rows of a that correspond to the direct children of the current
    node.

    Order is assumed to be traversal order, i.e. for the typical case of
    postorder traversal iterating over the items in the result and
    consolidating each time should give the same result as postorder
    traversal of the original tree. Should also be able to modify for
    preorder traversal.
    """
    # note: range ends with end+1, not end, b/c end is included
    return [(a[node], a[start:end+1]) for node, start, end in tree_index]


def traverse_reduce(bound_indices, f):
    """Applies a[i] = f(a[j:k]) over list of [(a[i], a[j:k])].

    If list is in traversal order, has same effect as consolidating the
    function over the tree, only much faster.

    Note that f(a[j:k]) must return an object that can be broadcast to the
    same shape as a[i], e.g. summing a 2D array to get a vector.
    """
    for i, s in bound_indices:
        i[:] = f(s, 0)


def bool_descendants(bound_indices):
    """For each internal node, sets col to True if any descendant is True."""
    traverse_reduce(bound_indices, np.logical_or.reduce)


# specific to weighted
def sum_descendants(bound_indices):
    """For each internal node, sets col to sum of values in descendants."""
    traverse_reduce(bound_indices, sum)


def bind_to_parent_array(t, a):
    """Binds tree to array a, returning result in list.

    Takes as input tree t with _leaf_index set.

    Returns list of (node_row, parent_row such that node_row points to the
    row of a that corresponds to the current row, and parent_row points to
    the row of the parent.

    Order will be preorder traversal, i.e. for propagating attributes from
    the root to the tip.

    Typical usage of this function is to set up an array structure for many
    preorder traversals on the same tree, especially where you plan to change
    the data between traversals.
    """
    result = []
    for n in t.traverse(self_before=True, self_after=False):
        if n is not t:
            result.append([a[n._leaf_index], a[n.parent._leaf_index]])
    return result


def unifrac_tasks_from_matrix(u, env_names):
    """Returns the UniFrac matrix, PCoA, and/or cluster from the matrix."""
    result = {}
    result['distance_matrix'] = (u, env_names)
    return result


def tip_distances(a, bound_indices, tip_indices):
    """Sets each tip to its distance from the root.
    Note: This will need its own unittest"""
    for i, s in bound_indices:
        i += s
    mask = np.zeros(len(a))
    np.put(mask, tip_indices, 1)
    a *= mask[:, np.newaxis]


def PD_whole_tree(t, envs):
    """Run PD on t and envs for each env.

    Note: this is specific for PD per se, use PD_generic_whole_tree if you
    want to calculate a related metric.
    """
    envs, count_array, unique_envs, env_to_index, node_to_index, env_names, \
        branch_lengths, nodes, t = _fast_unifrac_setup(t, envs)
    count_array = count_array.astype(bool)
    bound_indices = bind_to_array(nodes, count_array)
    # initialize result
    bool_descendants(bound_indices)
    result = (branch_lengths * count_array.T).sum(1)
    return unique_envs, result


def faith_pd(counts, otu_ids, tree):
    """skbio api"""
    envs = _skbio_counts_to_envs(otu_ids, u_counts, v_counts)

    return PD_whole_tree(tree, envs)[1]  # just return the result
