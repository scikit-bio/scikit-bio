# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True

from cython cimport floating
from cython.parallel cimport parallel, prange
from openmp cimport omp_get_max_threads
from libc.string cimport memmove
from libc.math cimport ldexp, ldexpf
from heapq import heappush


# This script hosts components of the greedy algorithms for phylogenetic reconstruction
# using the minimum evolution (ME) principle. Specifically, it supports GME and BME for
# de novo tree building, and FastNNI and BNNI for tree arrangement.

# Most functions don't create any intermediate arrays, but repeatedly use pre-allocated
# arrays for all operations. This improves computational efficiency.


# ------------------------------------------------------------------------------------
# NOTE on parallelization: The algorithms implemented here were designed to facilitate
# parallelization. Specifically, the preorder and postorder are stored, with which the
# code can identify the entire clade under any given node, without needing to redo the
# tree traversal. Operations that are topology-independent can be parallelized across
# nodes within the clade.
# 
# A challenge to parallelization is scheduling. Clades vary greatly in size, therefore
# iterations also have greatly varying compute loads. It is tricky to allocate chunks
# of iterations to individual threads such that all threads have roughly the same load.
# In the experimental code `_bal_avgdist_insert_p`, a "dynamic" policy with a chunk
# size inversely proportional to the clade size is adopted, thus fewer larger clades
# will be processed in each thread.
# ------------------------------------------------------------------------------------


def _get_num_threads():
    """Determine the number of threads that will be used by OpenMP."""
    return omp_get_max_threads()


cdef inline floating xldexp(floating x, int exp) noexcept nogil:
    """Calculate power of 2.

    This function is not used because tests found it is slower than pre-calculating
    all possible powers of 2 (`npots`) and lookup.

    Additionally, there is a bit hacking method to make this calculation even faster.

    """
    if floating is float:
        return ldexpf(x, exp)
    else:
        return ldexp(x, exp)


cdef (int, bint) config_prange(
    int ops,
    int chunksize,
    int minclade,
    bint adaptive = False,
) noexcept nogil:
    """Determine prange schedule.

    Returns
    -------
    chunksize : int
    use_threads : bint

    """
    cdef int chunk
    cdef bint use_threads

    # determine chunk size
    if adaptive:
        chunk = max(1, chunksize // ops)
    else:
        chunk = chunksize

    # determine whether to use threads
    use_threads = ops > max(chunk, minclade)

    # evenly distribute chunks
    # chunk = -(-ops // num_threads)
    return chunk, use_threads


def _preorder(
    Py_ssize_t[::1] order,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] stack,
    Py_ssize_t start = 0,
):
    r"""Perform preorder traversal.

    This function and :func:`_postorder` use stacks to avoid recursion. The stack array
    is pre-allocated. The output (ordered nodes) is also written into a pre-allocated
    array.

    This function and :func:`_postorder` are not actually used in the greedy algorithms
    in this module, which incrementally grow the tree as well as the orders. The two
    functions are implemented for reference and test purpose.

    """
    cdef Py_ssize_t curr, left
    cdef Py_ssize_t order_i = 0  # next index of order
    cdef Py_ssize_t stack_i = 1  # next index of stack
    stack[0] = start
    while stack_i:
        # pop a node from stack into order
        stack_i -= 1
        curr = stack[stack_i]
        order[order_i] = curr
        # index[curr] = order_i  # preorder index
        order_i += 1
        # append children to stack, right first such that left is processed first
        left = tree[curr, 0]
        if left:
            stack[stack_i] = tree[curr, 1]
            stack[stack_i + 1] = left
            stack_i += 2


def _postorder(
    Py_ssize_t[::1] order,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] stack,
    Py_ssize_t start = 0,
):
    """Perform postorder traversal.

    See also :func:`_preorder`.

    """
    cdef Py_ssize_t curr, last, left, right
    cdef Py_ssize_t order_i = 0
    cdef Py_ssize_t stack_i = 1
    cdef Py_ssize_t prev = 0  # last visited node
    stack[0] = start
    curr = tree[start, 0]
    while stack_i:
        if curr:
            stack[stack_i] = curr
            stack_i += 1
            curr = tree[curr, 0]  # go to left child
        else:
            last = stack[stack_i - 1]
            left, right = tree[last, 0], tree[last, 1]
            if left and prev != right:
                curr = right  # go to right child
            else:
                order[order_i] = last
                # index[last] = order_i  # postorder index
                order_i += 1
                prev = last
                stack_i -= 1


def _avgdist_matrix(
    floating[:, ::1] adm,
    floating[:, :] dm,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] postodr,
):
    r"""Calculate a matrix of average distances between all pairs of subtrees.

    This function will update adm, a float array of (n, n) representing pairwise
    distances between all nodes (tips, internal nodes and the root) in the tree.

    Implemented according to Appendix 4 of Desper and Gascuel (2002). Basically, this
    algorithm traverses the tree and calculates the average distance between each pair
    of subtrees. Here, a "subtree" is identified by a node, but it can be one of the
    two scenarios:

    0. Lower subtree (L): the subtree descending from a node (including the node).
    1. Upper subtree (U): the subtree branching from the node upward. The root of this
       subtree is the node's parent (NOT the node itself). Its immediate children are
       the parent's parent and the node's sibling.

    Then it iteratively applies Eq. 2 to calculate the average distance between two
    subtrees based on the average distances from the children (1 and 2) of one of the
    subtrees (B) to the other (A):

        d(A, B) = (|B_1| * d(A, B_1) + |B_2| * d(A, B_2)) / |B|

    TODO: It might be possible to only fill half of the matrix (upper or lower
    triangle). Same for other functions, especially those in BME.

    """
    cdef Py_ssize_t i, j, k
    cdef Py_ssize_t curr, parent, sibling, p_size, s_size
    cdef Py_ssize_t a, a_size, a_taxon
    cdef Py_ssize_t a1 = 0, a2 = 0, a1_size = 0, a2_size = 0
    cdef Py_ssize_t b, b_size, b1, b2

    # pointers to matrix rows to facilitate lookup of cells
    cdef floating* dm_0 = &dm[0, 0]
    cdef floating* dm_a
    cdef floating* adm_0 = &adm[0, 0]
    cdef floating* adm_a

    # total numbers of taxa and nodes in the tree
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3

    # Calculate the average distance between each pair of subtrees defined by nodes
    # a and b (A4.1 of the paper).

    # Step 1: Calculate non-nested subtree-to-subtree distances (i.e., one is not an
    # ancestor of another. Therefore each subtree is the lower (descending) tree of
    # the node. (A4.1 (a))
    # The paper suggests outer and inner postorder traverals for all nodes. Since
    # distances are symmetric, one can save half of the calculations, as done by the
    # following code.

    # Loop over nodes in postorder.
    # Skip the root, which is always the last node in the postorder traversal.
    for i in range(n - 1):
        a = postodr[i]        
        a_size = tree[a, 4]

        # check if a is a tip
        if a_size == 1:
            a_taxon = tree[a, 1]
        else:
            a_taxon = 0  # a can never be root (taxon 0), therefore 0 means none
            a1, a2 = tree[a, 0], tree[a, 1]
            a1_size, a2_size = tree[a1, 4], tree[a2, 4]

        dm_a = &dm[a_taxon, 0]
        adm_a = &adm[a, 0]

        # Iterate over the ancestors of a, and find the other subtree on the right.
        # We can skip the left subtree because it must have been calculated already.
        # If the right subtree contains a's ancestry, just skip.
        curr = a
        while curr:
            parent = tree[curr, 2]
            if tree[parent, 0] != curr:
                curr = parent
                continue
            sibling = tree[curr, 3]

            # Loop over nodes within the right subtree in postorder.
            # This postorder doesn't need to be re-calculated. Because all nodes within
            # a clade are continuous in postorder, one can take a slice of the full
            # postorder that represent the descending nodes of the current node. The
            # size of the slice is 2 x taxon count - 2.
            k = tree[sibling, 7]
            for j in range(k - tree[sibling, 4] * 2 + 2, k + 1):
                b = postodr[j]
                b_size = tree[b, 4]

                # If both a and b are tips, take the original taxon-to-taxon distance
                # (A4.1 (a) i).
                if b_size == 1 and a_taxon != 0:
                    adm_a[b] = adm[b, a] = dm_a[tree[b, 1]]

                # If a is an internal node, and b is either (a tip or an internal node),
                # calculate the average distance based on the two child subtrees of a
                # (A4.1 (a) ii).
                elif a_taxon == 0:
                    adm_a[b] = adm[b, a] = (
                        a1_size * adm[a1, b] + a2_size * adm[a2, b]
                    ) / a_size

                # If a is a tip, and b is an internal node, calculate the average
                # distance based on the two child subtrees of b (A4.1 (a) iii).
                else:
                    b1, b2 = tree[b, 0], tree[b, 1]
                    adm_a[b] = adm[b, a] = (
                        tree[b1, 4] * adm_a[b1] + tree[b2, 4] * adm_a[b2]
                    ) / b_size

            curr = parent

    # Step 2: Calculate subtree to root (taxon 0) distances (A4.1 (b)).

    # This is done through a postorder traversal.
    for i in range(n - 1):
        a = postodr[i]
        a_size = tree[a, 4]
        if a_size == 1:
            adm[a, 0] = adm_0[a] = dm_0[tree[a, 1]]
        else:
            a1, a2 = tree[a, 0], tree[a, 1]
            adm[a, 0] = adm_0[a] = (
                tree[a1, 4] * adm_0[a1] + tree[a2, 4] * adm_0[a2]
            ) / a_size

    # Step 3: Calculate nested subtree to subtree distances, in which the first node
    # (a) is a descendant of the second node (b), therefore the first subtree (A) is
    # the lower (descending) tree of node a, whereas the second subtree (B) is the
    # upper (ancestral) tree of node b (A4.1 (c)).

    # This is done through a preorder traversal.
    for i in range(1, n):
        a = preodr[i]
        parent, sibling = tree[a, 2], tree[a, 3]

        # The size of (upper) subtree b is the complement of its descendants. Same for
        # the parent subtree.
        a_size = m - tree[a, 4]
        p_size = m - tree[parent, 4]
        s_size = tree[sibling, 4]

        # Iterate over all subtrees below b.
        # The paper says this traversal can be done in any manner. Here, we use the
        # postorder. See * above.
        adm_a = &adm[a, 0]
        k = tree[a, 7]
        for j in range(k - tree[a, 4] * 2 + 2, k):
            b = postodr[j]
            adm_a[b] = adm[b, a] = (
                s_size * adm[b, sibling] + p_size * adm[b, parent]
            ) / a_size


def _bal_avgdist_matrix(
    floating[:, ::1] adm,
    floating[:, :] dm,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] postodr,
):
    r"""Calculate a matrix of balanced average distances between all pairs of subtrees.

    This function resembles :func:`_avgdist_matrix`, but it weighs subtrees equally
    regardless of their sizes. Specifically, it replaces Eq. 2 with Eq. 6. of Desper
    and Gascuel (2002):

        d(A, B) = (d(A, B_1) + d(A, B_2)) / 2

    Same for all functions starting with `bal_`.

    """
    cdef Py_ssize_t i, j, k
    cdef Py_ssize_t curr, parent, sibling
    cdef Py_ssize_t a, b, a_taxon
    cdef Py_ssize_t a1 = 0, a2 = 0

    cdef floating* dm_0 = &dm[0, 0]
    cdef floating* dm_a
    cdef floating* adm_0 = &adm[0, 0]
    cdef floating* adm_a

    cdef Py_ssize_t n = 2 * tree[0, 4] - 1

    # Step 1: Calculate non-nested subtree to subtree distances.
    for i in range(n - 1):
        a = postodr[i]
        if tree[a, 0] == 0:
            a_taxon = tree[a, 1]
        else:
            a_taxon = 0
            a1, a2 = tree[a, 0], tree[a, 1]
        dm_a = &dm[a_taxon, 0]
        adm_a = &adm[a, 0]
        curr = a
        while curr:
            parent = tree[curr, 2]
            if tree[parent, 0] != curr:
                curr = parent
                continue
            sibling = tree[curr, 3]
            k = tree[sibling, 7]
            for j in range(k - tree[sibling, 4] * 2 + 2, k + 1):
                b = postodr[j]
                if tree[b, 0] == 0 and a_taxon != 0:
                    adm_a[b] = adm[b, a] = dm_a[tree[b, 1]]
                elif a_taxon == 0:
                    adm_a[b] = adm[b, a] = 0.5 * (adm[a1, b] + adm[a2, b])
                else:
                    adm_a[b] = adm[b, a] = 0.5 * (
                        adm_a[tree[b, 0]] + adm_a[tree[b, 1]]
                    )
            curr = parent

    # Step 2: Calculate subtree to root distances.
    for i in range(n - 1):
        a = postodr[i]
        if tree[a, 0] == 0:
            adm[a, 0] = adm_0[a] = dm_0[tree[a, 1]]
        else:
            adm[a, 0] = adm_0[a] = 0.5 * (adm_0[tree[a, 0]] + adm_0[tree[a, 1]])

    # Step 3: Calculate nested subtree to subtree distances.
    for i in range(1, n):
        a = preodr[i]
        adm_a = &adm[a, 0]
        parent, sibling, k = tree[a, 2], tree[a, 3], tree[a, 7]
        for j in range(k - tree[a, 4] * 2 + 2, k):
            b = postodr[j]
            adm_a[b] = adm[b, a] = 0.5 * (adm[b, sibling] + adm[b, parent])


def _avgdist_taxon(
    floating[:, ::1] adk,
    Py_ssize_t taxon,
    floating[:, :] dm,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] postodr,
):
    """Calculate average distances between a new taxon (k) and existing subtrees.

    This function will update adk, a float array of (2, n) in which the two rows
    represent the average distances from the taxon to the lower and upper subtrees of
    each existing node, respectively.

    Implemented according to Appendix 3 of Desper and Gascuel (2002). Basically, this
    algorithm calculates all lower subtree distances via a postorder traversal, then
    calculates all upper subtree distances via a preorder traversal.

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, left, right, parent, sibling

    # total numbers of taxa and nodes in the tree
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3

    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]

    # Calculate the distance between taxon and the lower subtree of each node.
    for i in range(n):
        node = postodr[i]
        left, right = tree[node, 0], tree[node, 1]
        if left == 0:
            adkl[node] = dm[taxon, right]
        else:
            adkl[node] = (
                tree[left, 4] * adkl[left] + tree[right, 4] * adkl[right]
            ) / tree[node, 4]

    # Assign upper distance of root.
    adku[0] = dm[taxon, 0]

    # Calculate the distance between taxon k and the upper subtree of each node.
    for i in range(1, n):
        node = preodr[i]
        parent, sibling = tree[node, 2], tree[node, 3]
        adku[node] = (
            (m - tree[parent, 4]) * adku[parent] + tree[sibling, 4] * adkl[sibling]
        ) / (m - tree[node, 4])


def _bal_avgdist_taxon(
    Py_ssize_t n,
    Py_ssize_t k,
    floating[:, :] dm,
    floating[::1] adkl,
    floating[::1] adku,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] order,
):
    r"""Calculate balanced average distances between a new taxon and existing subtrees.

    This function resembles :func:`_avgdist_taxon` but uses the balanced framework.

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, left, right

    # Calculate distances to lower subtrees in reversed preorder (bottom-up).
    # This modifies the original algorithm, which is postorder. Both postorder and
    # reversed preorder guarantee that a node is visited after both of its children
    # have been visited.
    for i in range(n - 1, -1, -1):
        node = order[i]
        left, right = tree[node, 0], tree[node, 1]
        if left == 0:
            adkl[node] = dm[k, right]  # tip: adopt from `dm`
        else:
            adkl[node] = 0.5 * (adkl[left] + adkl[right])  # internal: mean of children

    # distance to root (upper subtree of taxon 0)
    adku[0] = dm[k, 0]

    # Calculate distances to upper subtrees in preorder (top-down).
    for i in range(1, n):
        node = order[i]
        adku[node] = 0.5 * (adku[tree[node, 2]] + adkl[tree[node, 3]])


def _bal_avgdist_taxon_cc(
    Py_ssize_t n,
    Py_ssize_t k,
    floating[:, ::1] dm,
    floating[::1] adkl,
    floating[::1] adku,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] order,
):
    r"""Calculate balanced average distances between a new taxon and existing subtrees.

    This function is identical to `_bal_avgdist_taxon` except that it assumes `dm` is
    C-contiguous, which makes lookup more efficient.

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, left, right
    cdef floating* dm_k = &dm[k, 0]
    for i in range(n - 1, -1, -1):
        node = order[i]
        left, right = tree[node, 0], tree[node, 1]
        if left == 0:
            adkl[node] = dm_k[right]
        else:
            adkl[node] = 0.5 * (adkl[left] + adkl[right])
    adku[0] = dm_k[0]
    for i in range(1, n):
        node = order[i]
        adku[node] = 0.5 * (adku[tree[node, 2]] + adkl[tree[node, 3]])


def _ols_lengths(
    floating[::1] lens,
    floating[:, ::1] adm,
    Py_ssize_t[:, ::1] tree,
):
    r"""Calculate branch lengths of a tree based on the OLS framework.

    Using an average distance matrix between all pairs of subtrees.

    Implemented according to Eqs. 3 & 4 of Desper and Gascuel (2002).

    TODO: Can be parallelized. Although this isn't a bottlenecking step.

    """
    cdef Py_ssize_t node, left, right, parent, sibling
    cdef Py_ssize_t l_size, r_size, p_size, s_size
    cdef floating lambda_
    cdef Py_ssize_t m = tree[0, 4] + 1

    for node in range(1, 2 * m - 3):
        left = tree[node, 0]
        parent = tree[node, 2]
        sibling = tree[node, 3]

        # External (terminal) branch: based on the triplet (iAB) of self (i), parent
        # (A), and sibling (B) (Eq. 4).
        if left == 0:
            lens[node] = 0.5 * (
                adm[parent, node] + adm[sibling, node] - adm[parent, sibling]
            )

        # Internal branch: based on the quartet (AB|CD) of parent (A, upper), sibling
        # (B), and children (C and D) (Eq. 3).
        else:
            right = tree[node, 1]
            l_size = tree[left, 4]
            r_size = tree[right, 4]
            p_size = m - tree[parent, 4]
            s_size = tree[sibling, 4]
            lambda_ = <floating>(p_size * r_size + s_size * l_size) / (
                (p_size + s_size) * (l_size + r_size)
            )
            lens[node] = 0.5 * (
                lambda_ * (adm[parent, left] + adm[sibling, right])
                + (1 - lambda_) * (adm[parent, right] + adm[sibling, left])
                - (adm[parent, sibling] + adm[left, right])
            )

    # root branch
    left, right = tree[0, 0], tree[0, 1]
    lens[0] = 0.5 * (adm[left, 0] + adm[right, 0] - adm[left, right])


def _ols_lengths_d2(
    floating[::1] lens,
    floating[:, ::1] ad2,
    Py_ssize_t[:, ::1] tree,
):
    r"""Calculate branch lengths of a tree based on an OLS framework.

    Using only average distances between pairs of distant-2 subtrees.

    This function produces the same result as `_ols_lengths`. The latter relies on
    Eq. 3 of Desper and Gascuel (2002), which involves distances between distant-3
    subtrees. Therefore, we need to modify the equation such that it only takes
    distances between distant-2 subtrees as input.

    Specifically, with the following local structure:

                |
              parent
              /   \
           node  sibling
           /  \
        left  right

    We will need to calculate the distances between left / right (lower) and parent
    (upper) / sibling (lower). This can be achieved by (according to Fig. 2a):

        d(left(L), parent(U)) = d(left(L), node(U)) + d(node(L), parent(U))
            - d(node(L), node(U))

    Do the same for d(right(L), parent(U)), d(left(L), sibling(L)), and d(right(L),
    sibling(L)). Plug results into Eq. 3, we will get:

        l(node) = 0.5 * (
            d(left(L), parent(U)) + d(right(L), parent(U)) + d(node(L), parent(U))
            + d(node(L), sibling(U)) - d(sibling(L), parent(U)) - d(left(L), right(L))
        ) - d(node(L), node(U))

    Note that this equation is free of the lambda factor.

    Here, d(node(L), node(U)) is the distance between the lower and upper subtrees of
    the same node. This can be calculated using the equation in A4.1(c):

        d(node(L), node(U)) = (|sibling(L)| * d(node(L), sibling(L)) + |parent(U)|
            * d(node(L), parent(U))) / |node(U)|

    Therefore, we will get l(node).

    """
    cdef Py_ssize_t node, left, sibling
    cdef Py_ssize_t m = tree[0, 4] + 1

    cdef floating* ad2l = &ad2[0, 0]
    cdef floating* ad2u = &ad2[1, 0]

    for node in range(1, 2 * m - 3):
        left = tree[node, 0]

        # External (terminal) branch: based on Eq. 4 of the paper.
        if left == 0:
            lens[node] = 0.5 * (ad2l[node] + ad2u[node] - ad2u[tree[node, 3]])

        # Internal branch: based on the equation discussed above.
        else:
            sibling = tree[node, 3]
            lens[node] = 0.5 * (
                ad2u[left]
                + ad2u[tree[node, 1]]
                + ad2u[node]
                + ad2l[node]
                - ad2u[sibling]
                - ad2l[left]
            ) - (
                tree[sibling, 4]
                * ad2l[node]
                + (m - tree[tree[node, 2], 4])
                * ad2u[node]
            ) / (m - tree[node, 4])

    # root branch
    left = tree[0, 0]
    lens[0] = 0.5 * (ad2u[left] + ad2u[tree[0, 1]] - ad2l[left])


def _bal_lengths(
    floating[::1] lens,
    floating[:, ::1] adm,
    Py_ssize_t[:, ::1] tree,
):
    r"""Calculate branch lengths of a tree based on the balanced framework.

    Using a balanced average distance matrix between all pairs of subtrees.

    This function resembles :func:`_ols_lengths` but it uses the balanced framework.

    Implemented according to Eqs. 3 & 4 and the description on top of pg. 691 of Desper
    and Gascuel (2002).

    """
    cdef Py_ssize_t n = tree.shape[0]
    cdef Py_ssize_t node, left, right, parent, sibling

    for node in range(1, n):
        left, parent, sibling = tree[node, 0], tree[node, 2], tree[node, 3]

        # node is tip
        if left == 0:
            if tree[parent, 0] != node:
                lens[node] = 0.5 * (
                    adm[parent, node] + adm[sibling, node] - adm[parent, sibling]
                )
            else:
                lens[node] = 0.5 * (
                    adm[parent, node] + adm[node, sibling] - adm[parent, sibling]
                )

        # node is internal
        else:
            right = tree[node, 1]

            # node is left, sibling is right
            if tree[parent, 0] != node:
                right = tree[node, 1]
                lens[node] = 0.25 * (
                    adm[parent, left]
                    + adm[sibling, right]
                    + adm[parent, right]
                    + adm[sibling, left]
                ) - 0.5 * (adm[parent, sibling] + adm[left, right])

            # sibling is left, node is right
            else:
                lens[node] = 0.25 * (
                    adm[parent, left]
                    + adm[right, sibling]
                    + adm[parent, right]
                    + adm[left, sibling]
                ) - 0.5 * (adm[parent, sibling] + adm[left, right])

    left, right = tree[0, 0], tree[0, 1]
    lens[0] = 0.5 * (adm[0, left] + adm[0, right] - adm[left, right])


def _ols_min_branch_d2(
    floating[::1] lens,
    floating[:, ::1] ad2,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
):
    r"""Find the branch with the minimum length change after inserting a new taxon.

    It returns the node at the lower end of the branch.

    Implemented according to Eq. 7 of Desper and Gascuel (2002).

    IMPORTANT NOTE: Should there are ties (which are common), this function returns
    the first minimum branch seen during preorder traversal. This behavior resembles
    FastME. Alternatively, one can return the first minimum branch by the order of
    addition using `return lens[:n].argmin()`, which could produce different results
    at the presence of ties. It must be noted that all candidates in a tie are equally
    optimal in the current iteration of the greedy algorithm.

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, parent, sibling, size, p_size, s_size
    cdef floating L, numerator, lambda_0, lambda_1

    cdef Py_ssize_t min_node = 0
    cdef floating min_len = 0
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3

    cdef floating* ad2l = &ad2[0, 0]
    cdef floating* ad2u = &ad2[1, 0]
    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]

    lens[min_node] = min_len

    # Traverse tree in preorder and calculate the length change of each branch from its
    # previous branch.
    for i in range(1, n):
        node = preodr[i]
        parent = tree[node, 2]
        sibling = tree[node, 3]
        size = tree[node, 4]
        p_size = m - tree[parent, 4]
        s_size = tree[sibling, 4]

        numerator = s_size + size * p_size
        lambda_0 = numerator / ((s_size + size) * (p_size + 1))
        lambda_1 = numerator / ((s_size + p_size) * (size + 1))

        # factor 0.5 is omitted
        lens[node] = L = lens[parent] + (
            (lambda_0 - lambda_1) * (adkl[sibling] + ad2u[node])
            + (lambda_1 - 1) * (ad2l[sibling] + adku[parent])
            + (1 - lambda_0) * (ad2u[sibling] + adkl[node])
        )

        if L < min_len:
            min_len, min_node = L, node

    return min_node


def _bal_min_branch(
    Py_ssize_t n,
    floating[::1] lens,
    floating[:, ::1] adm,
    floating[::1] adkl,
    floating[::1] adku,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] order,
):
    """Find the branch with the minimum length change after inserting a new taxon.

    This function resembles :func:`_ols_min_branch_d2` but it 1) uses the balanced
    framework and 2) calculates based on the entire matrix. See also the note of the
    latter.

    Implemented according to Eq. 10 of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, parent, sibling
    cdef floating length

    cdef Py_ssize_t min_i = 0
    cdef floating min_len = 0

    cdef floating cell

    # Identify the smallest branch length change through a preorder traversal.
    # NOTE: Length change of root (`lens[0]`) has been set to zero.
    for i in range(1, n):
        node = order[i]
        parent = tree[node, 2]
        sibling = tree[node, 3]
        if tree[parent, 0] == node:
            cell = adm[node, sibling]  # left child
        else:
            cell = adm[sibling, node]  # right child
        length = lens[parent] + (
            adm[parent, sibling] + adkl[node] - cell - adku[parent]
        )
        lens[node] = length
        if length < min_len:
            min_len, min_i = length, i

    return min_i


def _bal_min_branch2(
    Py_ssize_t n,
    floating[::1] lens,
    floating[:, ::1] adm,
    floating[::1] adkl,
    floating[::1] adku,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] order,
    Py_ssize_t[::1] sizes,
):
    """Find the branch with the minimum length change after inserting a new taxon.

    This is an alternative version that is faster, but it doesn't maintain strict
    preorder. Result may be different.

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, left, right
    cdef floating L_intm, L_left, L_right

    cdef Py_ssize_t min_i = 0
    cdef floating min_len = 0

    # Identify the smallest branch length change through a preorder traversal.
    # NOTE: Length change of root (`lens[0]`) has been set to zero.
    for i in range(n):
        node = order[i]
        if left := tree[node, 0]:  # internal node
            right = tree[node, 1]
            L_intm = lens[node] - adm[left, right] - adku[node]
            lens[left] = L_left = L_intm + adm[node, right] + adkl[left]
            lens[right] = L_right = L_intm + adm[node, left] + adkl[right]
            if L_left <= L_right:
                if L_left < min_len:
                    min_len, min_i = L_left, i + 1
            elif L_right < min_len:
                min_len, min_i = L_right, i + sizes[left] + 1

    return min_i


def _avgdist_d2_insert(
    floating[:, ::1] ad2,
    Py_ssize_t target,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
):
    r"""Update average distances between distant-2 subtrees after taxon insertion.

    This function will update `ad2`, a float array of (2, n) representing pairwise
    distances between all distant-2 subtrees in the tree. Here, `distant-2 subtrees`
    refer to subtrees that are two branches away from each other. Specifically, there
    are two scenarios:

    - Row 0: Distance between the lower subtree of the current node and the lower
      subtree of its parent.
    - Row 1: Distance between the lower subtree of the current node and the upper
      subtree of its sibling.

    This function assumes that the taxon will be inserted into the branch connecting
    the target node and its parent. After insertion, the taxon will become the sibling
    of the target.

                                    parent
             parent                  /  \
             /  \        =>       link  sibling
        target  sibling           /  \
                             target  taxon (k)

    This function should be executed *before* calling :func:`_insert_taxon`, which will
    mutate the tree.

    Implemented according to Eq. 8 of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t i, ii
    cdef Py_ssize_t node, left, right, parent, sibling, size, curr, p_size, size_1

    # dimensions and positions
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    cdef floating* ad2l = &ad2[0, 0]
    cdef floating* ad2u = &ad2[1, 0]
    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]

    ###### Special case: insert into the root branch. ######

    if target == 0:
        # k (lower) to root (parent, upper): pre-calculated
        ad2u[tip] = adku[0]

        # k (lower) to link (sibling, lower): equals to k to root (lower).
        ad2l[link] = ad2l[tip] = adkl[0]

        # Link (lower) to root (parent, upper): de novo calculation according to the
        # equation in A4.1(b). It is basically the distance between the upper and lower
        # subtrees of the root itself.
        left, right = tree[0, 0], tree[0, 1]
        ad2u[link] = (
            tree[left, 4] * ad2u[left] + tree[right, 4] * ad2u[right]
        ) / tree[0, 4]

        # Calculate all node (lower) to parent (upper, containing k) distances. These
        # parents include the new link.
        for node in range(1, n):
            p_size = m - tree[tree[node, 2], 4]
            ad2u[node] = (adkl[node] + ad2u[node] * p_size) / (p_size + 1)

        return

    ###### Regular case: insert into any other branch. ######

    parent, sibling, size = tree[target, 2], tree[target, 3], tree[target, 4]
    size_1 = size + 1

    # Temporarily copy the distances of target to link (will edit later).
    ad2u[link] = ad2u[target]
    ad2l[link] = ad2l[target]

    # Distance between k (lower) and link (parent, upper) equals to that between k and
    # the upper subtree of target.
    ad2u[tip] = adku[target]

    # Distance between target (lower) and link (parent, upper) needs to be calculated
    # using the equation in A4.1(c). Basically, it is the distance between the lower
    # and upper subtrees of the same target.
    ad2u[target] = (
        tree[sibling, 4] * ad2l[target] + (m - tree[parent, 4]) * ad2u[target]
    ) / (m - size)

    # Transfer the pre-calculated distance between target (lower) and k (sibling,
    # lower).
    ad2l[target] = ad2l[tip] = adkl[target]

    # Within the clade below target, calculate the distance between each node (lower)
    # and its parent (upper, containing k).
    ii = tree[target, 6]
    for i in range(ii + 1, ii + size * 2 - 1):
        node = preodr[i]
        p_size = m - tree[tree[node, 2], 4]
        ad2u[node] = (adkl[node] + ad2u[node] * p_size) / (p_size + 1)

    # Iterate over the ancestors of target, starting from link and ending at root.
    curr = link
    while curr:
        # Calculate the distance between each pair of lower (containing k) and upper
        # ancestors.
        ad2u[curr] = (adku[parent] + ad2u[curr] * size) / size_1

        # Calculate the distance between each ancestor (lower, containing k) and its
        # sibling (lower).
        ad2l[curr] = ad2l[sibling] = (adkl[sibling] + ad2l[curr] * size) / size_1

        # Within the clade below each sibling, calculate the distance between each node
        # (lower) and its parent (upper, containing k).
        ii = tree[sibling, 6]
        for i in range(ii + 1, ii + tree[sibling, 4] * 2 - 1):
            node = preodr[i]
            p_size = m - tree[tree[node, 2], 4]
            ad2u[node] = (adkl[node] + ad2u[node] * p_size) / (p_size + 1)

        curr = parent
        parent, sibling, size = tree[curr, 2], tree[curr, 3], tree[curr, 4]
        size_1 = size + 1


### Old serial version (double assignment) ###

# def _bal_avgdist_insert(
#     floating[:, ::1] adm,
#     Py_ssize_t target,
#     floating[:, ::1] adk,
#     Py_ssize_t[:, ::1] tree,
#     Py_ssize_t[::1] postodr,
#     floating[::1] npots,
#     Py_ssize_t[::1] stack,
#     Py_ssize_t[::1] paths,
#     Py_ssize_t[::1] lvls,
# ):
#     r"""Update balanced average distance matrix after taxon insertion.

#     This function resembles :func:`_avgdist_d2_insert` but it 1) uses the balanced
#     framework and 2) updates the entire matrix. The latter makes it the dominant term
#     of the entire algorithm.

#     Two additional parameters are provided: `npots` is a pre-calculated array of
#     2^(-l) powers (l is the depth difference between two nodes). `stack` is an
#     integer array to store ancestral nodes of target.

#     """
#     cdef Py_ssize_t i, j, ii, jj, anc_i
#     cdef Py_ssize_t parent, sibling, depth
#     cdef Py_ssize_t curr, anc, cousin, depoff, depth_diff
#     cdef Py_ssize_t a, b
#     cdef floating npot, diff

#     # dimensions and positions
#     cdef Py_ssize_t n = 2 * tree[0, 4] - 1
#     cdef Py_ssize_t link = n
#     cdef Py_ssize_t tip = n + 1

#     cdef floating* adkl = &adk[0, 0]
#     cdef floating* adku = &adk[1, 0]

#     cdef floating* npots_2 = &npots[2]

#     ###### Special case: insert into the root branch. ######

#     if target == 0:
#         # Transfer distance between k and root (upper).
#         adm[0, tip] = adm[tip, 0] = adku[0]

#         # k to link: equals to k to root (lower).
#         adm[link, tip] = adm[tip, link] = adkl[0]

#         # Root to link: de novo calculation according to the equation in A4.1(b). It is
#         # basically the distance between the upper and lower subtrees of the root.
#         a1, a2 = tree[0, 0], tree[0, 1]
#         adm[0, link] = adm[link, 0] = 0.5 * (adm[a1, 0] + adm[a2, 0])

#         # Iterate over all nodes but the root.
#         for a in range(1, n):

#             # Transfer distances between the node (lower) and k.
#             adm[a, tip] = adm[tip, a] = adkl[a]
            
#             # Calculate the distance between the node (lower) and link (upper, with two
#             # taxa (0 and k) added) using Eq. 8.
#             adm[a, link] = adm[link, a] = 0.5 * (adkl[a] + adm[a, 0])

#             # Calculate the distances between the node (upper, containing k) and each
#             # of its descendant (lower).
#             ii = tree[a, 7]
#             npot = npots[tree[a, 5] + 1]
#             for i in range(ii - tree[a, 4] * 2 + 2, ii):
#                 b = postodr[i]
#                 adm[a, b] = adm[b, a] = adm[a, b] + npot * (adkl[b] - adm[0, b])

#         return

#     ###### Regular case: insert into any other branch. ######

#     parent, sibling, depth = tree[target, 2], tree[target, 3], tree[target, 5]

#     ### Step 1: Distances within the clade below target. ###

#     # Locate the clade below target (including target)
#     ii = tree[target, 7]
#     for i in range(ii - tree[target, 4] * 2 + 2, ii + 1):
#         a = postodr[i]

#         # Transfer pre-calculated distance between k (lower) and any node within the clade
#         # (including target, lower).
#         adm[a, tip] = adm[tip, a] = adkl[a]

#         # Distance from any descendant (lower) to link (upper) equals to that to target.
#         # (The last assignment: target to link is incorrect. It will be fixed below.)
#         adm[a, link] = adm[link, a] = adm[a, target]

#         # Within the clade, find all ancestor (a) - descendant (b) pairs, and calculate the
#         # distance between the upper subtree of a (containing k) and the lower subtree of b.
#         jj = tree[a, 7]
#         npot = npots[tree[a, 5] - depth + 1]
#         for j in range(jj - tree[a, 4] * 2 + 2, jj):
#             b = postodr[j]
#             adm[a, b] = adm[b, a] = adm[a, b] + npot * (adkl[b] - adm[target, b])

#     ### Step 2: Distances around the insertion point. ###

#     # Distance between k (lower) and link (upper) equals to that between k and the
#     # upper subtree of target.
#     adm[tip, link] = adm[link, tip] = adku[target]

#     # Distance between target (lower) and link (upper) needs to be calculated using the
#     # equation in A4.1(c). Basically, it is the distance between the lower and upper
#     # subtrees of the same target.
#     adm[target, link] = adm[link, target] = 0.5 * (
#         adm[target, sibling] + adm[target, parent]
#     )

#     ### Step 3: Distances among nodes outside the clade. ###

#     # Iterate over ancestors of target in ascending order.
#     anc_i = 0
#     curr = target
#     while curr:
#         stack[anc_i] = anc = tree[curr, 2]
#         depth_diff = depth - 2 * tree[anc, 5]

#         # Transfer the pre-calculated distance between k and the ancestor (upper).
#         adm[anc, tip] = adm[tip, anc] = adku[anc]

#         # Calculate the distance between link (lower, containing k) and the ancestor
#         # (upper).
#         adm[anc, link] = adm[link, anc] = 0.5 * (adku[anc] + adm[anc, target])

#         # Calculate the distance between each previous ancestor (lower, containing k)
#         # and the current ancestor (upper).
#         diff = adku[anc] - adm[target, anc]
#         for i in range(anc_i):
#             a = stack[i]
#             adm[anc, a] = adm[a, anc] = adm[anc, a] + npots[i + 2] * diff

#         # Identify the cousin clade descending from the ancestor.
#         cousin = tree[curr, 3]
#         ii = tree[cousin, 7]
#         for i in range(ii - tree[cousin, 4] * 2 + 2, ii + 1):
#             a = postodr[i]

#             # Transfer the pre-calculated distances between k and each descendant
#             # (lower).
#             adm[a, tip] = adm[tip, a] = adkl[a]

#             # Calculate the distance between link (lower, containing k) and each
#             # descendant (lower).
#             adm[a, link] = adm[link, a] = 0.5 * (adkl[a] + adm[a, target])

#             # Calculate the distance between each previous ancestor (lower, containing
#             # k) and each descendant (lower).
#             diff = adkl[a] - adm[a, target]
#             for j in range(anc_i):
#                 b = stack[j]
#                 adm[a, b] = adm[b, a] = adm[a, b] + npots[j + 2] * diff

#             # Iterate over descendants of each member of the clade, and calculate the
#             # distance between the former (upper, containing k) and the latter (lower).
#             jj = tree[a, 7]
#             npot = npots[depth_diff + tree[a, 5]]
#             for j in range(jj - tree[a, 4] * 2 + 2, jj):
#                 b = postodr[j]
#                 adm[a, b] = adm[b, a] = adm[a, b] + npot * (adkl[b] - adm[b, target])

#         curr = anc
#         anc_i += 1


def _bal_avgdist_insert(
    Py_ssize_t n,
    Py_ssize_t tag_i,
    floating[:, ::1] adm,
    floating[::1] adkl,
    floating[::1] adku,
    floating[::1] diffs,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] order,
    Py_ssize_t[::1] sizes,
    Py_ssize_t[::1] depths,
    floating[::1] npots,
    Py_ssize_t[::1] ancs,
    Py_ssize_t[::1] ancx,
):
    r"""Update balanced average distance matrix after taxon insertion.

    This function resembles :func:`_avgdist_d2_insert` but it 1) uses the balanced
    framework and, 2) updates the entire matrix (`adm`), which makes it the dominant
    term of the entire algorithm.

    This function is the serial version. A parallel version can be found in
    :func:`_bal_avgdist_insert_p`.

    ***

    This function assumes that the taxon will be inserted into the branch connecting
    the target node and its parent. After insertion, the taxon will become the right
    sibling of the target.

                                    parent
             parent                  /  \
             /  \        =>       link  sibling
        target  sibling           /  \
                             target  taxon (k)

    This function only updates distance [i, j] but not [j, i], saving half of the
    assignment operations, which are expensive given the size of the matrix. i and j
    suffice i < j in a preorder traversal (parent - left - right), where 1) ancestors
    are always earlier than their descendants, and 2) left children (and their
    descendants) are always earlier than right children (and their descendants).
    Therefore, the order of i, j can be easily determined during tree traversal.

                 anc3 (root)
                 /   \
              anc2    x
             /   \
            +    anc1
           /     /   \
          +   link    x
              /  \
         target  taxon (k)
           /
          o
         / \
        o   o

    The tree can be partitioned into two parts: 1) the clade descending from target
    (o) and 2) the clades that are cousins to the target, which in turn are divided
    into left cousins (+) and right cousins (x). These clades can be identified by
    navigating all ancestors of target up to root.

    0. For a tree with N taxa, the number of ancestors (i.e., depth) of any node is
       between O(logN) (balanced tree, best case) and O(N) (skewed tree, worst case).

    The distances to be updated can be divided into two compute-intensive categories:

    1. Within each clade (o, + and x), the distance between the upper subtree of each
       ancestor and the lower subtree of each of its descendants needs to be updated.
       The total number of ancestor-descendant pairs (i.e., sum of depths) is between
       O(nlogn) (balanced tree, best case) and O(n^2) (skewed tree, worst case).

    2. Within each cousin clade (+ and x), the distance between the lower subtree of
       each node and the lower subtree of each ancestor of target needs to be updated.

    ***

    Several additional arrays are involved:

    - `sizes` stores the number of nodes (including self) within the clade descending
      from each node. It facilitates the determination of range of node constituting
      the clade.

    - `depths` stores the number of branches from each node to root (i.e., depth). It
      facilitates the calculation of path distance between two arbitrary nodes.

    - `npots` is a pre-calculated array of 2^(-l) powers (l is the depth difference
      between two nodes).

    - `ancs` is an integer array storing ancestral nodes of target.

    - `ancx` is an integer array storing (index * stride) of ancestral nodes of
      target. It facilitates locating individual cells in `adm`.

    """
    # node indices as in tree
    cdef Py_ssize_t a, b
    cdef Py_ssize_t tag = order[tag_i]   # target node
    cdef Py_ssize_t lnk = n              # link node
    cdef Py_ssize_t kay = n + 1          # new taxon (k)
    cdef Py_ssize_t par = tree[tag, 2]   # parent
    cdef Py_ssize_t sib = tree[tag, 3]   # sibling
    cdef Py_ssize_t cur                  # current
    cdef Py_ssize_t anc                  # ancestor
    cdef Py_ssize_t cuz                  # cousin (not used)
    cdef Py_ssize_t lft                  # left child

    # node indices as in order
    cdef Py_ssize_t i, j
    cdef Py_ssize_t cur_i, anc_i, cuz_i

    # sizes (number of nodes in clade)
    cdef Py_ssize_t size = sizes[tag]
    cdef Py_ssize_t cur_s, cuz_s

    # depths (number of branches to root)
    cdef Py_ssize_t depth = depths[tag]

    # intermediate variables
    cdef floating cell    # specific value in `adm` (a matrix)
    cdef floating diff    # difference of distance; usually from `adk`
    cdef floating npot    # negative power of 2; from `npots`

    # pointers to specific rows in `adm`
    cdef floating* adm_0 = &adm[0, 0]
    cdef floating* adm_t = &adm[tag, 0]
    cdef floating* adm_l = &adm[lnk, 0]
    cdef floating* adm_k = &adm[kay, 0]
    cdef floating* adm_r  # ancestor
    cdef floating* adm_a  # arbitrary node

    cdef Py_ssize_t deg  # number of branches between two nodes
    cdef Py_ssize_t level  # number of branches above target
    cdef Py_ssize_t stride = adm.shape[1]  # width of `adm`
    cdef floating* npots_2 = &npots[2]  # power array offset

    sizes[kay] = 1  # k is a tip with size = 1

    ##### Special case: insert into the root branch. #####

    #    root (target)      root
    #    /  \               /  \
    # left  right   =>   link  taxon (k)
    #                    /  \
    #                 left  right

    if tag_i == 0:
        # Update sizes and depths (see illustration above).
        sizes[0] = n + 2
        sizes[lnk] = n
        depths[lnk] = depths[kay] = 1

        # Transfer distance between k and root (upper).
        adm_t[kay]= adku[0]

        # k to link: equals to k to root (lower).
        adm_l[kay] = adkl[0]

        # Root to link: de novo calculation according to the equation in A4.1(b). It is
        # basically the distance between the upper and lower subtrees of the root.
        adm_t[lnk] = 0.5 * (adm_t[tree[0, 0]] + adm_t[tree[0, 1]])

        # Iterate over all nodes but root (first in preorder). It is important to do
        # this in reverse preorder, such that the intermediates (see below) of all
        # descendants of a node are already calculated when they are used.
        for i in range(n - 1, 0, -1):
            a = order[i]
            adm_a = &adm[a, 0]

            # Transfer distances between the node (lower) and k.
            adm_a[kay] = diff = adkl[a]
            cell = adm_t[a]

            # Calculate the distance between the node (lower) and link (upper, with two
            # taxa (0 and k) added) using Eq. 8.
            adm_l[a] = 0.5 * (diff + cell)

            # Calculate this intermediate and store in `adkl`, such that it won't need
            # to be repeatedly calculated in the loop below. `adkl` will be completely
            # re-calculated by `_bal_avgdist_taxon`, so it is fine to use it to store
            # intermediates.
            diffs[i] = diff - cell

            # Update depth (entire tree +1 deeper).
            # NOTE: Although depth is unnecessary for tips, having an `if` check here
            # may not benefit performance because interleaved internal/tips are bad for
            # branch prediction.
            depths[a] = deg = depths[a] + 1

            # Calculate the distances between the node (upper, containing k) and each
            # of its descendants (lower). We can now reuse the intermediate calculated
            # above. The order of iteration is not important.
            npot = npots[deg]
            for j in range(i + 1, i + sizes[a]):
                adm_a[order[j]] += npot * diffs[j]

            # NOTE: `npots[deg] * adkl[b]` can be replaced with `xldexp(adkl[b], -deg)`
            # but tests found the current method is faster.

        return

    ##### Regular case: insert into any other branch. #####

    #     parent               parent
    #     /    \               /    \
    # target  sibling   =>   link  sibling
    #                        /  \
    #                   target  taxon (k)

    ### Step 1: Distances around the insertion point. ###

    # Update sizes and depths.
    sizes[lnk] = size + 2
    depths[lnk] = depth = depths[tag]
    depths[kay] = depths[tag] = depth + 1

    # Distance between k (lower) and link (upper) equals to that between k and the
    # upper subtree of target.
    adm_l[kay] = adku[tag]

    # Distance between target (lower) and link (upper) needs to be calculated using the
    # equation in A4.1(c). Basically, it is the distance between the lower and upper
    # subtrees of the same target.
    par, sib = tree[tag, 2], tree[tag, 3]
    cell = adm[tag, sib] if tree[par, 0] == tag else adm[sib, tag]
    adm_l[tag] = 0.5 * (cell + adm[par, tag])

    # Transfer pre-calculated distance between target (lower) and k (lower).
    adm_t[kay] = adkl[tag]

    ### Step 2: Distances within the clade below target. ###

    # Locate the clade below target (excluding target).
    for i in range(tag_i + size - 1, tag_i, -1):
        a = order[i]
        adm_a = &adm[a, 0]

        # Transfer pre-calculated distance between k (lower) and any node within
        # the clade (lower).
        adm_a[kay] = diff = adkl[a]

        # Distance between any descendant (lower) and link (upper) equals to that to
        # target.
        adm_l[a] = cell = adm_t[a]

        # Calculate the distance between each descendant (lower) and target (upper).
        adm_t[a] = 0.5 * (diff + cell)

        diffs[i] = diff - cell  # intermediate
        depths[a] = deg = depths[a] + 1  # depth +1

        # Within the clade, find all ancestor (a) - descendant (b) pairs, and calculate
        # the distance between a (upper, with k) and b (lower).
        npot = npots[deg - depth]
        for j in range(i + 1, i + sizes[a]):
            adm_a[order[j]] += npot * diffs[j]

    ### Step 3: Distances among nodes outside the clade. ###

    # Iterate over ancestors of target in ascending order till root.
    level = 0
    deg = 2 - depth
    cur_i = tag_i
    cur_s = size
    while cur_i:
        cur = order[cur_i]

        # Determine ancestor node and its pointer in `adm`.
        ancs[level] = anc = tree[cur, 2]
        ancx[level] = stride * anc
        adm_r = &adm[anc, 0]

        # Transfer the pre-calculated distance between k and ancestor (upper).
        adm_r[kay] = diff = adku[anc]

        # Calculate the distance between link (lower, with k) and ancestor (upper).
        cell = adm_r[tag]
        adm_r[lnk] = 0.5 * (diff + cell)

        # Calculate the distance between each previous ancestor (lower, with k)
        # and the current ancestor (upper).
        diff -= cell
        for i in range(level):
            adm_r[ancs[i]] += npots_2[i] * diff

        # Determine whether cousin is the right or left child of the shared ancestor.
        # This helps to determine the order of coordinates of each distance to be
        # updated. Basically, if the left child of parent is current, then current is
        # left and cousin is right. In this case, coordinates are [current, cousin].
        if (lft := tree[anc, 0]) == cur:

            # Determine the indices of ancestor and cousin in preorder. The current
            # solution is awkward as it is separate from the determination process in
            # tree, due to the lack of mapping from tree indices to preorder indices.

            # Since current is left, parent is -1 from current, and right cousin is
            # current + its size.
            anc_i = cur_i - 1
            cuz_i = cur_i + cur_s
            cuz_s = sizes[order[cuz_i]]

            # Locate the cousin clade descending from the ancestor (including cousin).
            for i in range(cuz_i + cuz_s - 1, cuz_i - 1, -1):
                a = order[i]
                adm_a = &adm[a, 0]

                # Transfer the pre-calculated distances between k and each descendant
                # (lower).
                adm_k[a] = diff = adkl[a]

                # Calculate the distance between link (lower, with k) and each
                # descendant (lower).
                cell = adm_t[a]
                adm_l[a] = 0.5 * (diff + cell)

                diffs[i] = diff = diff - cell  # intermediate

                # Calculate the distance between each previous ancestor (lower, with k
                # and each descendant (lower).
                for j in range(level):
                    # This is equivalent to: adm[ancs[j], a] += ... but faster (maybe)
                    adm_0[ancx[j] + a] += npots_2[j] * diff

                # Calculate the distance between a (upper, with k) and each of its
                # descendants (lower).
                npot = npots[depths[a] + deg]
                for j in range(i + 1, i + sizes[a]):
                    adm_a[order[j]] += npot * diffs[j]

        # Cousin is left, and coordinates are [cousin, current]. Otherwise the same.
        else:

            # ancestor index = current index - cousin size - 1
            cuz_s = sizes[lft]
            anc_i = cur_i - cuz_s - 1

            # This is equivalent to range(cuz_i + cuz_s - 1, cuz_i - 1, -1), in which
            # cuz_i = anc_i + 1.
            for i in range(anc_i + cuz_s, anc_i, -1):
                a = order[i]
                adm_a = &adm[a, 0]
                diff, cell = adkl[a], adm_a[tag]
                adm_a[kay] = diff
                adm_a[lnk] = 0.5 * (diff + cell)
                diffs[i] = diff = diff - cell
                for j in range(level):
                    adm_a[ancs[j]] += npots_2[j] * diff
                npot = npots[depths[a] + deg]
                for j in range(i + 1, i + sizes[a]):
                    adm_a[order[j]] += npot * diffs[j]

        # level up
        cur_i = anc_i
        cur_s = sizes[anc]
        sizes[anc] += 2  # every ancestor size +2 (link and k added)
        level += 1
        deg += 2


def _count_pairs(
    Py_ssize_t n,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] order,
    Py_ssize_t[::1] sizes,
    Py_ssize_t[::1] pairs,
):
    """Calculate the number of ancestor-descendant pairs within each clade.

    This function is executed after completing the serial phase and before entering the
    parallel phase. The `pairs` array is useful in the distribution of workloads across
    threads. It can be updated accumulatively as the tree grows. This function merely
    creates the initial status. Additionally, this function is useful for test purpose.

    The total number of ancestor-descendant pairs equals to the sum of depths of nodes
    within a clade.

    """
    cdef Py_ssize_t i, a, lft, rgt
    for i in range(n - 1, -1, -1):  # traverse tree in reverse preorder
        a = order[i]
        if lft := tree[a, 0]:
            rgt = tree[a, 1]
            pairs[a] = pairs[lft] + pairs[rgt] + sizes[lft] + sizes[rgt]
        else:
            pairs[a] = 0


def _bal_insert_plan(
    Py_ssize_t n,
    Py_ssize_t tag_i,
    floating[:, ::1] adm,
    floating[::1] adkl,
    floating[::1] adku,
    floating[::1] diffs,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] order,
    Py_ssize_t[::1] sizes,
    Py_ssize_t[::1] pairs,
    Py_ssize_t[::1] depths,
    Py_ssize_t[::1] ancs,
    Py_ssize_t[::1] ancx,
    Py_ssize_t[::1] segs,
    Py_ssize_t[::1] lvls,
    Py_ssize_t[::1] oops,
    bint flat=True,
):
    r"""Update balanced average distance matrix after taxon insertion.

    ancs : int[:]
        Ancestors of target (from parent to root, in ascending order).

    segs : int[:]
        Indices of nodes that mark the starts of segments.
    lvls : int[:]
        Number of generations from target to shared ancestor per segment.
        Target is -1, target's parent is 0, target's grandparent is 1,...

    It is easy to know that there can be no more than n segments (when every node is
    its own segment). The actual bounds are:

    - min. = depth - 1

        When all ancestors are right children of their parents, such that their left
        cousins do not need to be marked.

    - max. = depth * 2 - 1

        When all ancestors are left children, and their right cousins are marked.

    When the tree is extremely skewed and all ancestors are left children, there will
    be n - 2 segments, because target's parent always has 3 children and they are the
    same segment.

    Therefore, filling index n_segs to mark the end of the segment array is safe.

    ***

    description : root anc3 anc2 anc1 anc0 target rcuz1 rcuz3     n
          depth :    0    1    2    3    4      5     4     2
        segment :    0    5   12   19   31     42    65    73   100
          level :    4    3    2    1    0     -1     1     3
       ancestor :    4    3    2    1    0

    """
    cdef Py_ssize_t a  # nodes as in tree
    cdef Py_ssize_t tag = order[tag_i]   # target node
    cdef Py_ssize_t lnk = n              # link node
    cdef Py_ssize_t kay = n + 1          # new taxon (k)
    cdef Py_ssize_t par = tree[tag, 2]   # parent
    cdef Py_ssize_t sib = tree[tag, 3]   # sibling
    cdef Py_ssize_t cur                  # current
    cdef Py_ssize_t anc                  # ancestor
    cdef Py_ssize_t cuz                  # cousin
    cdef Py_ssize_t lft                  # left child

    cdef Py_ssize_t i  # nodes as in order
    cdef Py_ssize_t cur_i, anc_i, cuz_i

    cdef Py_ssize_t size = sizes[tag]
    cdef Py_ssize_t cuz_s
    cdef Py_ssize_t depth = depths[tag]

    cdef floating cell, diff  # intermediates

    cdef floating* adm_t = &adm[tag, 0]
    cdef floating* adm_l = &adm[lnk, 0]
    cdef floating* adm_k = &adm[kay, 0]
    cdef floating* adm_r
    cdef floating* adm_a

    cdef Py_ssize_t level
    cdef Py_ssize_t stride = adm.shape[1]

    # Accumulative number of horizontal operations to be assigned to each node in the
    # spine. This number += (cousin size + 1) * (level - 1) at each node.
    cdef Py_ssize_t achori

    # Accumulative number of vertical operations to be subtracted from each node in the
    # spine. This number += (size - 1) at each node.
    cdef Py_ssize_t acvert

    # Left and right indices in the segment array. Ancestors are appended to the left
    # side of the array. Right cousins are appended to the right side.
    cdef Py_ssize_t li, ri

    pairs[kay] = 0
    sizes[kay] = 1

    ###### Special case: insert into the root branch. ######

    if tag_i == 0:
        oops[0] = pairs[0] - sizes[0] + 1

        pairs[lnk] = pairs[0]
        pairs[0] += sizes[0] + 1

        sizes[0] = n + 2
        sizes[lnk] = n
        depths[lnk] = depths[kay] = 1

        # There is only one segment, which is the entire tree.
        segs[0] = 0
        segs[1] = n
        lvls[0] = -1

        adm_t[kay] = adku[0]
        adm_l[kay] = adkl[0]
        adm_t[lnk] = 0.5 * (adm_t[tree[0, 0]] + adm_t[tree[0, 1]])

        # Unlike the serial version, order is not important here.
        if flat:
            for i in range(1, n):
                a = order[i]
                diff, cell = adkl[a], adm_t[a]
                adm[a, kay] = diff
                adm_l[a] = 0.5 * (diff + cell)
                diffs[i] = diff - cell
                depths[a] += 1
                # this +1 here is important; it cannot be postponed to `_insert_taxon`
                # TODO: double-check the maths

        return

    ###### Regular case: insert into any other branch. ######

    sizes[lnk] = size + 2
    pairs[lnk] = pairs[tag] + size + 1
    depths[lnk] = depth
    depths[kay] = depths[tag] = depth + 1

    # accumulative horizontal and vertical operations
    achori = 0
    acvert = size - 1

    # Target clade workload = number of pairs - (number of nodes - 1). This is because
    # target node itself is not involved in calculation.
    oops[depth] = pairs[tag] - acvert

    # Depth equals to the number of ancestors of target. Therefore, target is placed at
    # index = depth in the list of segment.
    segs[depth] = tag_i
    lvls[depth] = -1
    li = depth - 1
    ri = depth + 1

    ### Step 1: Distances around the insertion point. ###

    adm_l[kay] = adku[tag]
    cell = adm_t[sib] if tree[par, 0] == tag else adm[sib, tag]
    adm_l[tag] = 0.5 * (cell + adm[par, tag])
    adm_t[kay] = adkl[tag]

    ### Step 2: Distances within the clade below target. ###

    if flat:
        for i in range(tag_i + 1, tag_i + size):
            a = order[i]
            diff, cell = adkl[a], adm_t[a]
            adm[a, kay] = diff
            adm_l[a] = cell
            adm_t[a] = 0.5 * (diff + cell)
            diffs[i] = diff - cell
            depths[a] += 1

    ### Step 3: Distances among nodes outside the clade. ###

    level = 0  # number of generations up from target
    cur_i = tag_i
    while cur_i:
        cur = order[cur_i]
        ancs[level] = anc = tree[cur, 2]
        ancx[level] = stride * anc

        # current is left, cousin is right
        if (lft := tree[anc, 0]) == cur:
            anc_i = cur_i - 1
            cuz_i = cur_i + sizes[cur]
            cuz = order[cuz_i]
            cuz_s = sizes[cuz]

            segs[li] = anc_i
            lvls[li] = level

            # Add right cousin to segments.
            segs[ri] = cuz_i
            lvls[ri] = level
            # TODO: oops is essentially the same as regular nodes. Can it be skipped?
            oops[ri] = cuz_s * level + pairs[cuz]
            ri += 1

            if flat:
                for i in range(cuz_i, cuz_i + cuz_s):
                    a = order[i]
                    adm_a = &adm[a, 0]
                    diff, cell = adkl[a], adm_t[a]
                    adm_k[a] = diff
                    adm_l[a] = 0.5 * (diff + cell)
                    diffs[i] = diff - cell

        # cousin is left, current is right
        else:
            cuz_i = cur_i - sizes[lft]
            anc_i = cuz_i - 1
            cuz = order[cuz_i]
            cuz_s = sizes[cuz]

            segs[li] = anc_i
            lvls[li] = level

            if flat:
                for i in range(cuz_i, cuz_i + cuz_s):
                    a = order[i]
                    adm_a = &adm[a, 0]
                    diff, cell = adkl[a], adm_a[tag]
                    adm_a[kay] = diff
                    adm_a[lnk] = 0.5 * (diff + cell)
                    diffs[i] = diff - cell

        if flat:
            adm_r = &adm[anc, 0]
            diff, cell = adku[anc], adm_r[tag]
            adm_r[kay] = diff
            adm_r[lnk] = 0.5 * (diff + cell)
            diffs[anc_i] = diff - cell

        # Transfer upper distance to lower distance, which will be utilized later.
        else:
            adkl[anc] = adku[anc]

        cur_i = anc_i

        # Update accumulative size and pair.
        # +1 is because of the ancestor itself is added to the clade.
        achori += (cuz_s + 1) * level
        # -1 is because only descendants but not the clade root are counted.
        acvert += sizes[anc] - 1
        oops[li] = achori + pairs[anc] - acvert

        li -= 1
        level += 1

    segs[ri] = n  # seal the segments
    # Otherwise, since segs is reused across iterations, a later iteration may have
    # a shorter segs than a previous one, and the remnant indices in the previous
    # seg could confuse. n is a value that indices will never reach, therefore safe.

    ### update sizes of spine (let's postpone this)
    # for i in range(level):
    #     anc = ancs[i]
    #     sizes[anc] += 2
    #     pairs[anc] += size + 2 * i + 3


# def _bal_insert_mark(
#     Py_ssize_t n,
#     Py_ssize_t tag_i,
#     floating[:, ::1] adm,
#     floating[::1] adkl,
#     floating[::1] adku,
#     Py_ssize_t[:, ::1] tree,
#     Py_ssize_t[::1] order,
#     Py_ssize_t[::1] sizes,
#     Py_ssize_t[::1] pairs,
#     Py_ssize_t[::1] depths,
#     Py_ssize_t[::1] ancs,
#     Py_ssize_t[::1] ancx,
#     Py_ssize_t[::1] segs,
#     Py_ssize_t[::1] lvls,
#     Py_ssize_t[::1] oops,
# ):
#     r"""Update balanced average distance matrix after taxon insertion."""
#     # nodes as in tree
#     cdef Py_ssize_t tag = order[tag_i]   # target node
#     cdef Py_ssize_t lnk = n              # link node
#     cdef Py_ssize_t kay = n + 1          # new taxon (k)
#     cdef Py_ssize_t par = tree[tag, 2]   # parent
#     cdef Py_ssize_t sib = tree[tag, 3]   # sibling
#     cdef Py_ssize_t cur                  # current
#     cdef Py_ssize_t anc                  # ancestor
#     cdef Py_ssize_t cuz                  # cousin
#     cdef Py_ssize_t lft                  # left child

#     # nodes as in order
#     cdef Py_ssize_t cur_i, anc_i, cuz_i

#     cdef Py_ssize_t size = sizes[tag]
#     cdef Py_ssize_t cuz_s
#     cdef Py_ssize_t depth = depths[tag]

#     cdef floating* adm_t = &adm[tag, 0]
#     cdef floating* adm_l = &adm[lnk, 0]

#     cdef Py_ssize_t level
#     cdef Py_ssize_t stride = adm.shape[1]

#     # Accumulative number of horizontal operations to be assigned to each node in the
#     # spine. This number += (cousin size + 1) * (level - 1) at each node.
#     cdef Py_ssize_t achori

#     # Accumulative number of vertical operations to be subtracted from each node in the
#     # spine. This number += (size - 1) at each node.
#     cdef Py_ssize_t acvert

#     # Left and right indices in the segment array. Ancestors are appended to the left
#     # side of the array. Right cousins are appended to the right side.
#     cdef Py_ssize_t li, ri

#     pairs[kay] = 0
#     sizes[kay] = 1

#     ###### Special case: insert into the root branch. ######

#     if tag_i == 0:
#         oops[0] = pairs[0] - sizes[0] + 1

#         pairs[lnk] = pairs[0]
#         pairs[0] += sizes[0] + 1

#         sizes[0] = n + 2
#         sizes[lnk] = n
#         depths[lnk] = depths[kay] = 1

#         # There is only one segment, which is the entire tree.
#         segs[0] = 0
#         segs[1] = n
#         lvls[0] = -1

#         adm_t[kay] = adku[0]
#         adm_l[kay] = adkl[0]
#         adm_t[lnk] = 0.5 * (adm_t[tree[0, 0]] + adm_t[tree[0, 1]])

#         return

#     ###### Regular case: insert into any other branch. ######

#     sizes[lnk] = size + 2
#     pairs[lnk] = pairs[tag] + size + 1
#     depths[lnk] = depth
#     depths[kay] = depths[tag] = depth + 1

#     # accumulative horizontal and vertical operations
#     achori = 0
#     acvert = size - 1

#     # Target clade workload = number of pairs - (number of nodes - 1). This is because
#     # target node itself is not involved in calculation.
#     oops[depth] = pairs[tag] - acvert

#     # Depth equals to the number of ancestors of target. Therefore, target is placed at
#     # index = depth in the list of segment.
#     segs[depth] = tag_i
#     lvls[depth] = -1
#     li = depth - 1
#     ri = depth + 1

#     # Update distances around the insertion point.

#     adm_l[kay] = adku[tag]
#     adm_l[tag] = 0.5 * (adm[par, tag] + (
#         adm_t[sib] if tree[par, 0] == tag else adm[sib, tag]
#     ))
#     adm_t[kay] = adkl[tag]

#     # Mark "spine" of the tree.

#     level = 0  # number of generations up from target
#     cur_i = tag_i
#     while cur_i:
#         cur = order[cur_i]
#         ancs[level] = anc = tree[cur, 2]
#         ancx[level] = stride * anc

#         # Transfer upper distance to lower distance, which will be utilized later.
#         adkl[anc] = adku[anc]

#         # current is left, cousin is right
#         if (lft := tree[anc, 0]) == cur:
#             anc_i = cur_i - 1
#             cuz_i = cur_i + sizes[cur]
#             cuz = order[cuz_i]
#             cuz_s = sizes[cuz]

#             segs[li] = anc_i
#             lvls[li] = level

#             # Add right cousin to segments.
#             segs[ri] = cuz_i
#             lvls[ri] = level
#             # TODO: oops is essentially the same as regular nodes. Can it be skipped?
#             oops[ri] = cuz_s * level + pairs[cuz]
#             ri += 1

#         # cousin is left, current is right
#         else:
#             cuz_i = cur_i - sizes[lft]
#             anc_i = cuz_i - 1
#             cuz = order[cuz_i]
#             cuz_s = sizes[cuz]

#             segs[li] = anc_i
#             lvls[li] = level

#         cur_i = anc_i

#         # Update accumulative size and pair.
#         # +1 is because of the ancestor itself is added to the clade.
#         achori += (cuz_s + 1) * level
#         # -1 is because only descendants but not the clade root are counted.
#         acvert += sizes[anc] - 1
#         oops[li] = achori + pairs[anc] - acvert

#         li -= 1
#         level += 1

#     segs[ri] = n


def _update_size_pair(
    Py_ssize_t tag_i,
    Py_ssize_t depth,
    Py_ssize_t[::1] sizes,
    Py_ssize_t[::1] pairs,
    Py_ssize_t[::1] ancs,
):
    """Update size and pair of nodes along the spine."""
    cdef Py_ssize_t i
    cdef Py_ssize_t size = sizes[tag_i]

    for i in range(depth):
        anc = ancs[i]
        sizes[anc] += 2
        pairs[anc] += size + 2 * i + 3


def _chunk_nodes(
    Py_ssize_t n,
    Py_ssize_t goal,
    Py_ssize_t tag_i,
    Py_ssize_t[::1] order,
    Py_ssize_t[::1] sizes,
    Py_ssize_t[::1] pairs,
    Py_ssize_t[::1] chunks,
    Py_ssize_t[::1] segs,
    Py_ssize_t[::1] lvls,
    Py_ssize_t[::1] oops,
    Py_ssize_t[::1] chusegs,
):
    """Partition tree into chunks of nodes with roughly even workloads.

    This function traverses the tree in preorder and iteratively puts single nodes or
    the entire clade under each node into chunks with a fixed upper limit of workload
    (number of operations). Specifically, it puts a clade into a chunk when there is
    sufficient space left to hold it. Otherwise, it decomposes the clade into smaller
    clades or nodes and tries to put them in subsequent iterations.

    There are two types of workloads:

    1. "Vertical": Ancestor (upper) - descendant (lower) pairs. Applies to all nodes
       except for target and it ancestors.

       For each node, the workload is the number of nodes descending from it, i.e.,
       n - 1.

       For each clade, the workload is the total number of ancestor-descendant pairs,
       which is equal to the sum of depths of all nodes. This number varies between
       O(nlogn) (balanced tree) and O(n^2) (skewed tree). An incrementally updated
       array, `pairs`, stores these numbers.

    2. "Horizontal": direct relative (lower) - cousin (lower) pairs. Applies to all
       nodes outside the clade descending from target's parent.

       For each node, the workload is the number of more recent ancestors descending
       from the ancestor shared with target. This value is stored in `lvls`.

       For each clade, the workload is the number of nodes times `gen`.

    Workloads of segment clades are special and cannot be calculated using the above
    rules. The array `oops` stores the correct values.

    ***

    Because the workload of a clade can be calculated without traversing it, the code
    can "leap" to the next clade when it determines that the current clade can fit into
    a chunk without further decomposing it. Therefore, this function is between O(1)
    (entire tree is a chunk) and O(n) (each node is a chunk).

    ***

    Although each chunk has a fixed capacity. The actual workload per chunk can vary
    below or above this capacity. Specifically:

    1. When the current chunk still has room but the current node cannot fit into it,
       the chunk will be closed and the remaining room won't be utilized. This is the
       most common scenario.

    2. When a single node has a large workload that exceeds the chunk capacity, it will
       be forced to occupy one chunk. This is because each node is the minimum unit of
       parallelization, and cannot further break down. This scenario is rare but could
       happen in theory.

    """
    cdef Py_ssize_t a
    cdef Py_ssize_t s  # number of nodes within each clade
    cdef Py_ssize_t new_i

    # workload (number of operations) per node and per clade (subtree)
    cdef Py_ssize_t node_ops, tree_ops

    # Total workload of the tree was pre-calculated and stored at root.
    cdef Py_ssize_t total_ops = oops[0]

    # chunk capacity
    # TODO: replace `1` with a minimum number
    cdef Py_ssize_t capacity = max(1, total_ops // goal)

    ### Initialization: Create the first chunk, and put the first node (root) in it.
    # chusegs[0] = 0               # next segment start current chunk (done already)
    cdef Py_ssize_t n_chunks = 1   # current number of chunks
    cdef Py_ssize_t lvl = lvls[0]  # current segment's level
    cdef Py_ssize_t ii = 1         # next segment index
    cdef Py_ssize_t seg = segs[1]  # next segment (node)
    cdef Py_ssize_t i = 1          # current node index
    cdef Py_ssize_t space = capacity - lvl + 1  # remaining space of current chunk

    while i < n:
        a = order[i]
        s = sizes[a]
        # if i == seg and i < index:
        #     s -= 2  # size +2 for ancestors already !!!!!!!!!!!!!!!!

        ### Calculate workload per node & clade. ###
        # entering a new segment
        if i == seg:
            lvl = lvls[ii]

            # when segment start is target or an ancestor of it
            if i <= tag_i:
                node_ops = lvl
                tree_ops = oops[ii]

            # when segment start is a right cousin
            else:
                node_ops = s + lvl - 1
                tree_ops = pairs[a] + s * lvl

            # Locate next segment start
            # ii += 1
            # seg = segs[ii]

        # regular scenario (same as right cousin)
        else:
            node_ops = s + lvl - 1
            tree_ops = pairs[a] + s * lvl

        ### Put workload into a chunk. ###
        # put whole clade (and leap to next clade)
        if tree_ops <= space:
            space -= tree_ops
            i += s

        # put node only, and decompose
        elif node_ops <= space:
            space -= node_ops
            i += 1

        # not enough space; has to create a new chunk
        else:

            # put whole clade to next chunk
            if tree_ops <= capacity:
                space = capacity - tree_ops
                new_i = i + s

            # put node to next chunk and decompose
            else:
                space = capacity - node_ops
                new_i = i + 1

            # close current chunk and move to next chunk
            chunks[n_chunks] = i
            chusegs[n_chunks] = ii
            n_chunks += 1

            i = new_i

        # segment +1
        # if is_seg:
        #     lvl = lvls[ii]
        #     ii += 1
        #     seg = segs[ii]

        ### Locate next segment start. ###
        # This loop is necessary only when: 1) current node is an ancestor of target,
        # otherwise the entire clade should not cross segments. 2) the whole clade was
        # leaped over, otherwise node index only +1, and one doesn't need a loop to
        # find the next segment. However, having another `if` switch here may not
        # increase performance.
        # This task can be more efficiently achieved with a binary tree search, through
        # the `bisect` module or code from scratch. However, it is anticipated that the
        # next segment is relatively close to the current node. So this optimization
        # may not be necessary.
        # TODO: Revisit later.
        while seg < i:
            lvl = lvls[ii]
            ii += 1
            seg = segs[ii]

        # this might be faster:
        # while segs[ii] < i:
        #     ii += 1
        # lvl = lvls[ii - 1]
        # seg = segs[ii]

    chunks[n_chunks] = n  # upper bound of last chunk
    return n_chunks


def _bal_insert_pass(
    Py_ssize_t n,
    Py_ssize_t tag_i,
    Py_ssize_t aft_i,
    Py_ssize_t[::1] order,
    floating[:, ::1] adm,
    floating[::1] adkl,
    Py_ssize_t[::1] depths,
):
    cdef Py_ssize_t tag = order[tag_i]
    cdef Py_ssize_t lnk = n
    cdef Py_ssize_t kay = n + 1
    cdef Py_ssize_t i, a
    cdef floating diff, cell

    cdef floating* adm_t = &adm[tag, 0]
    cdef floating* adm_l = &adm[lnk, 0]
    cdef floating* adm_k = &adm[kay, 0]
    cdef floating* adm_a

    if tag_i == 0:
        for i in prange(1, n, nogil=True):
            a = order[i]
            diff = adkl[a]
            cell = adm_t[a]
            adm[a, kay] = diff
            adm_l[a] = 0.5 * (diff + cell)
            adkl[a] = diff - cell
            depths[a] += 1
        return

    for i in prange(n, nogil=True):
        a = order[i]

        # before target clade (ancestors and left cousins)
        if i < tag_i:
            adm_a = &adm[a, 0]
            diff = adkl[a]
            cell = adm_a[tag]
            adm_a[kay] = diff
            adm_a[lnk] = 0.5 * (diff + cell)

        # target (skip)
        elif i == tag_i:
            continue
        else:
            diff = adkl[a]
            cell = adm_t[a]

            # within target clade
            if i < aft_i:
                adm[a, kay] = diff
                adm_l[a] = cell
                adm_t[a] = 0.5 * (diff + cell)
                depths[a] += 1

            # after target clade (right cousins)
            else:
                adm_k[a] = diff
                adm_l[a] = 0.5 * (diff + cell)

        adkl[a] = diff - cell


# def _bal_avgdist_fill(
#     Py_ssize_t n,
#     Py_ssize_t tag_i,
#     Py_ssize_t depth,
#     floating[:, ::1] adm,
#     floating[::1] adkl,
#     Py_ssize_t[::1] order,
#     Py_ssize_t[::1] sizes,
#     Py_ssize_t[::1] depths,
#     floating[::1] npots,
#     Py_ssize_t[::1] ancs,
#     Py_ssize_t[::1] segs,
#     Py_ssize_t[::1] lvls,
# ):
#     r"""Update balanced average distance matrix after taxon insertion.

#     Fill direct (L) - cousin (L) pairs.

#     """
#     cdef Py_ssize_t a, b  # nodes as in tree
#     cdef Py_ssize_t i, j  # nodes as in order
#     cdef floating diff
#     cdef floating* adm_a
#     cdef floating* npots_2 = &npots[2]
#     cdef Py_ssize_t size
#     cdef Py_ssize_t depth_2 = depth - 2

#     # Initial status
#     cdef Py_ssize_t ii = 0
#     cdef Py_ssize_t seg = segs[ii]
#     cdef Py_ssize_t level = -1
#     cdef Py_ssize_t deg = 2 * level - depth_2

#     for i in range(n):
#         a = order[i]
#         adm_a = &adm[a, 0]
#         diff = adkl[a]

#         if i <= tag_i:
#             if i == seg:
#                 level = lvls[ii]
#                 ii += 1
#                 seg = segs[ii]
#                 deg = 2 * level - depth_2

#             # if not in spine, trigger vertical fill
#             else:
#                 if (size := sizes[a]) > 1:
#                     npot = npots[depths[a] + deg]
#                     for j in range(i + 1, i + size):
#                         b = order[j]
#                         adm_a[b] += npot * adkl[b]

#             # horizontal fill is guaranteed
#             for j in range(level):
#                 adm_a[ancs[j]] += npots_2[j] * diff

#         else:
#             if i == seg:
#                 level = lvls[ii]
#                 ii += 1
#                 seg = segs[ii]
#                 deg = 2 * level - depth_2

#             # both fills are guaranteed
#             if (size := sizes[a]) > 1:
#                 npot = npots[depths[a] + deg]
#                 for j in range(i + 1, i + size):
#                     b = order[j]
#                     adm_a[b] += npot * adkl[b]

#             for j in range(level):
#                 adm[ancs[j], a] += npots_2[j] * diff


def _bal_avgdist_fill(
    Py_ssize_t n,
    Py_ssize_t tag_i,
    Py_ssize_t aft_i,
    Py_ssize_t depth,
    floating[:, ::1] adm,
    floating[::1] adkl,
    floating[::1] diffs,
    Py_ssize_t[::1] order,
    Py_ssize_t[::1] sizes,
    Py_ssize_t[::1] depths,
    floating[::1] npots,
    Py_ssize_t[::1] ancs,
    Py_ssize_t[::1] ancx,
    Py_ssize_t[::1] segs,
    Py_ssize_t[::1] lvls,
    Py_ssize_t n_chunks,
    Py_ssize_t[::1] chunks,
    Py_ssize_t[::1] chusegs,
    bint flat=False,
):
    r"""Update balanced average distance matrix after taxon insertion."""
    cdef Py_ssize_t a, b  # nodes as in tree
    cdef Py_ssize_t i, j  # nodes as in order

    cdef Py_ssize_t tag = order[tag_i]
    cdef Py_ssize_t lnk = n
    cdef Py_ssize_t kay = n + 1

    cdef floating diff, cell

    cdef floating* adm_0 = &adm[0, 0]
    cdef floating* adm_t = &adm[tag, 0]
    cdef floating* adm_l = &adm[lnk, 0]
    cdef floating* adm_k = &adm[kay, 0]
    cdef floating* adm_a

    cdef floating* npots_2 = &npots[2]
    cdef Py_ssize_t size

    cdef Py_ssize_t depth_2 = depth - 2  # intermediate

    cdef Py_ssize_t c    # chunk index
    cdef Py_ssize_t ii   # segment index
    cdef Py_ssize_t seg  # segment bound

    # number of branches from insertion point to shared ancestor
    cdef Py_ssize_t level

    # number of branches from insertion point to any node
    cdef Py_ssize_t deg

    with nogil, parallel():

        ### Flat phase - O(n) ###
        if flat:
            if tag_i == 0:
                for i in prange(1, n):
                    a = order[i]
                    diff = adkl[a]
                    cell = adm_t[a]
                    adm[a, kay] = diff
                    adm_l[a] = 0.5 * (diff + cell)
                    diffs[i] = diff - cell
                    depths[a] += 1

            else:
                for i in prange(n):
                    a = order[i]
                    diff = adkl[a]

                    # before target clade (ancestors and left cousins)
                    if i < tag_i:
                        adm_a = &adm[a, 0]
                        cell = adm_a[tag]
                        adm_a[kay] = diff
                        adm_a[lnk] = 0.5 * (diff + cell)
                        diffs[i] = diff - cell

                    # skip target
                    elif i > tag_i:
                        cell = adm_t[a]

                        # within target clade
                        if i < aft_i:
                            adm[a, kay] = diff
                            adm_l[a] = cell
                            adm_t[a] = 0.5 * (diff + cell)
                            depths[a] += 1

                        # after target clade (right cousins)
                        else:
                            adm_k[a] = diff
                            adm_l[a] = 0.5 * (diff + cell)

                        diffs[i] = diff - cell

        ### Nested phase - O(nlogn) to O(n^2) ###

        for c in prange(n_chunks, schedule="dynamic"):
            ii = chusegs[c]  # target segment
            seg = segs[ii]
            level = lvls[ii - 1] if ii > 0 else -1
            deg = 2 * level - depth_2

            for i in range(chunks[c], chunks[c + 1]):
                a = order[i]
                adm_a = &adm[a, 0]
                diff = diffs[i]

                # Current node is at or before target in preorder. This means two things:
                # 1) It is before any ancestor lower than the ancestor shared between it
                #    and target. Therefore, the horizontally filled cell coordinates must
                #    be [a, anc].
                # 2) Every change point it encounters is target or an ancestor, which are
                #    skipped during vertical filling.
                # This `if` would cost little compute due to branch prediction.
                if i <= tag_i:
                    if i == seg:
                        level = lvls[ii]

                        # NOTE: This cannot be written as `ii += 1`, otherwise Cython will
                        # complain: "Cannot read reduction variable in loop body."
                        ii = ii + 1
                        seg = segs[ii]
                        deg = 2 * level - depth_2

                    # if not in spine, trigger vertical fill
                    else:
                        # NOTE: There is a dilemma here: the current code has an `if` which
                        # isn't good for branch prediction as tips (size = 1) are frequent.
                        # But removing this `if` will add two array reads to this memory-
                        # bound algorithm. Test suggests keeping the current code is fine.
                        if (size := sizes[a]) > 1:
                            npot = npots[depths[a] + deg]
                            for j in range(i + 1, i + size):
                                adm_a[order[j]] += npot * diffs[j]

                    # horizontal fill is guaranteed
                    for j in range(level):
                        adm_a[ancs[j]] += npots_2[j] * diff
                        # adm_a[ancs[j]] += ldexp(diff, -2 - j)

                # Current node is after target. This means: 1) Horizontal Cell coordinates
                # are [anc, a]. 2) Change points are right cousins.
                else:
                    if i == seg:
                        level = lvls[ii]
                        ii = ii + 1
                        seg = segs[ii]
                        deg = 2 * level - depth_2

                    # both fills are guaranteed
                    if (size := sizes[a]) > 1:
                        npot = npots[depths[a] + deg]
                        for j in range(i + 1, i + size):
                            adm_a[order[j]] += npot * diffs[j]

                    for j in range(level):
                        adm_0[ancx[j] + a] += npots_2[j] * diff
                        # adm_0[ancx[j] + a] += ldexp(diff, -2 - j)


def _insert_taxon(
    Py_ssize_t taxon,
    Py_ssize_t tag_i,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] sizes,
):
    r"""Insert a taxon between a target node and its parent.

    For example, with the following local structure of the original tree:

          A
         / \
        B   C

    With target=B, this function inserts a taxon into the branch A-B. The structure
    becomes:

            A
           / \
        link  C
         / \
        B  taxon

    The link and taxon will be appended to the end of the tree array, but the pre- and
    postorders need to be muted such that new nodes can be inserted. Specifically:

        Preorder:  A - B - C => A - link - B - taxon - C
        Postorder: B - C - A => B - taxon - link - C - A

    A special case is that the taxon is inserted into the root branch (node=0). The
    tree becomes:

            A
           / \
        link taxon
         / \
        B   C

    The inserted taxon always becomes the right child.

    """
    # This function can be simplified by Python and NumPy APIs. Although I hoped that
    # NumPy vectorization can accelerate the code, especially the pre- and postorder
    # parts, the reality according to my tests is that Cython code is significantly
    # faster than NumPy, and greatly reduces the overall runtime.

    cdef Py_ssize_t left, right, parent, sibling, size
    cdef Py_ssize_t side

    cdef Py_ssize_t after  # preorder index of node immediately after clade

    # determine tree dimensions
    # typically n = 2 * taxon - 3, but this function doesn't enforce this
    cdef Py_ssize_t m = taxon - 1
    cdef Py_ssize_t n = m * 2 - 1
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    cdef Py_ssize_t target = preodr[tag_i]

    cdef Py_ssize_t* node

    cdef size_t intsize = sizeof(Py_ssize_t)

    # Special case (root branch): taxon k becomes the sibling of all existing taxa
    # except for the root (taxon 0).
    if target == 0:
        node = &tree[0, 0]

        # children
        left, right = node[0], node[1]
        tree[left, 2] = tree[right, 2] = link

        # root
        node[0] = link
        node[1] = tip
        # sizes[0] = n + 2

        # link
        node = &tree[link, 0]
        node[0] = left
        node[1] = right
        node[2] = 0
        node[3] = tip
        # sizes[link] = n

        # tip
        node = &tree[tip, 0]
        node[0] = 0
        node[1] = taxon
        node[2] = 0
        node[3] = link
        # sizes[tip] = 1

        # preorder
        # for i in range(1, n):
        #     preodr[i + 1] = preodr[i]
        memmove(&preodr[2], &preodr[1], <size_t>((n - 1) * intsize))

        preodr[1] = link
        preodr[n + 1] = tip

        # postorder
        # postodr[n - 1] = link
        # postodr[n] = tip
        # postodr[n + 1] = 0

        # entire tree depth + 1
        # if use_depth:
        #     depths[link] = 1
        #     depths[tip] = 1
        #     for i in range(1, n):
        #         depths[i] += 1

    # Regular case (any other branch): The link becomes the parent of the target node,
    # and child of its original parent. Taxon k becomes the sibling
    else:
        node = &tree[target, 0]
        left = node[0]
        right = node[1]
        parent = node[2]
        sibling = node[3]
        size = sizes[target]

        side = int(tree[parent, 0] != target)
        tree[parent, side] = link
        tree[sibling, 3] = link
        node[2] = link
        node[3] = tip

        # preorder index of node after clade
        after = tag_i + size

        # link
        node = &tree[link, 0]
        node[0] = target
        node[1] = tip
        node[2] = parent
        node[3] = sibling
        # sizes[link] = size + 2

        # tip
        node = &tree[tip, 0]
        node[0] = 0
        node[1] = taxon
        node[2] = link
        node[3] = target
        # sizes[tip] = 1

        # clade depth +1
        # if use_depth:
        #     depth = depths[target]
        #     depths[link] = depth
        #     depths[tip] = depth + 1
        #     for i in range(index, after):
        #         depths[preodr[i]] += 1

        # preorder shift: nodes after clade +2, tip inserted after clade, nodes within
        # clade +1, link inserted before clade
        # for i in range(after, n):
        #     preodr[i + 2] = preodr[i]
        memmove(&preodr[after + 2], &preodr[after], <size_t>(
            (n - after) * intsize
        ))
        preodr[after + 1] = tip
        # for i in range(index, after):
        #     preodr[i + 1] = preodr[i]
        memmove(&preodr[tag_i + 1], &preodr[tag_i], <size_t>(
            (after - tag_i) * intsize
        ))
        preodr[tag_i] = link

        # postorder shift: all nodes after clade +2, tip and link inserted after clade
        # for i in range(n - 1, post_i, -1):
        #     k = postodr[i]
        #     tree[k, 7] += 2
        #     postodr[i + 2] = k
        # postodr[post_i + 2] = link
        # postodr[post_i + 1] = tip

        # size +2 from parent to root
        # curr = link
        # while curr:
        #     parent = tree[curr, 2]
        #     sizes[parent] += 2
        #     curr = parent


def _avgdist_swap(
    floating[:, ::1] adm,
    Py_ssize_t target,
    Py_ssize_t side,
    Py_ssize_t[:, ::1] tree,
):
    r"""Update average distance matrix after branch swapping.

    This function will update adm, a float array of (n, n) representing pairwise
    distances between all subtrees in the tree.

    It assumes that one child of the target node on a certain side (left: 0, right: 1)
    exchanges position with the sibling of the target:

                 |                         |
               parent                    parent
               /   \                     /   \
           target  sibling    =>     target  child
            /  \                      /  \
        other  child              other  sibling

    Two categories of average distances will be updated:

        1. From target (upper) to any node within other and sibling.
        2. From target (lower) to any node within parent and child.

    Only target to itself doesn't need to be updated. All other nodes in the tree do.

    This function should be executed *after* calling :func:`_swap_branches`.

    Implemented according to A4.3(b) of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t node

    # total numbers of taxa and nodes in the tree
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3

    # the two branches that were swapped
    cdef Py_ssize_t child = tree[target, 3]  # former child, now sibling
    cdef Py_ssize_t sibling = tree[target, side]  # former sibling, now child

    # the two branches that stay still
    cdef Py_ssize_t other = tree[target, 1 - side]  # the other child
    cdef Py_ssize_t parent = tree[target, 2]  # the parent

    cdef Py_ssize_t c_size = tree[child, 4]
    cdef Py_ssize_t o_size = tree[other, 4]
    cdef Py_ssize_t s_size = tree[sibling, 4]
    cdef Py_ssize_t p_size = m - tree[parent, 4]

    # Loop over all nodes except for target. These nodes can be divided into ones
    # within the clade below target (former sibling and other), and ones that are
    # outside the clade (former child and parent). Their distances to the target will
    # be updated separately.
    cdef Py_ssize_t start = tree[target, 6]
    cdef Py_ssize_t end = start + tree[target, 4] * 2 - 1

    # # 1) subtrees within (parent, child) vs sibling (lower) + other (lower)
    # for node in range(start):
    #     adm[node, target] = adm[target, node] = (
    #         s_size * adm[node, sibling] + o_size * adm[node, other]
    #     ) / (s_size + o_size)

    # # 2) # subtrees within (sibling, other) vs parent (upper) + child (lower)
    # for node in range(start + 1, end):
    #     adm[node, target] = adm[target, node] = (
    #         p_size * adm[node, parent] + c_size * adm[node, child]
    #     ) / (c_size + p_size)

    # # 3) same as 1)
    # for node in range(end, n):
    #     adm[node, target] = adm[target, node] = (
    #         s_size * adm[node, sibling] + o_size * adm[node, other]
    #     ) / (s_size + o_size)

    cdef floating temp_val = adm[target, target]

    for node in range(n):
        # subtrees within (sibling, other) vs parent (upper) + child (lower)
        if start < tree[node, 6] < end:
            adm[node, target] = adm[target, node] = (
                p_size * adm[node, parent] + c_size * adm[node, child]
            ) / (c_size + p_size)

        # subtrees within (parent, child) vs sibling (lower) + other (lower)
        else:
            adm[node, target] = adm[target, node] = (
                s_size * adm[node, sibling] + o_size * adm[node, other]
            ) / (s_size + o_size)

    # this cell doesn't need to be filled and let's reset it
    adm[target, target] = temp_val


def _bal_avgdist_swap(
    floating[:, ::1] adm,
    Py_ssize_t target,
    Py_ssize_t side,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    floating[::1] npots,
    Py_ssize_t[::1] stack,
):
    r"""Update balanced average distance matrix after branch swapping.

    This function resembles :func:`_avgdist_swap`, but it uses a balanced framework.
    Specifically, it follows Eq. 18 and Appendix 5.3 of Desper and Gascuel (2002).

    This function is the dominant term of the entire BNNI algorithm.

    """
    cdef Py_ssize_t node, curr, before, after, cousin, a, b
    cdef Py_ssize_t i, j, ii, jj, anc_i, start, end
    cdef Py_ssize_t depth, depth_2, depth_diff
    cdef floating npot, diff

    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3
    cdef Py_ssize_t child = tree[target, 3]
    cdef Py_ssize_t sibling = tree[target, side]
    cdef Py_ssize_t other = tree[target, 1 - side]
    cdef Py_ssize_t parent = tree[target, 2]

    cdef floating* admx

    # Step 1: Update distances between subtrees within each of the four clades (Eq. 18
    # and A5.3(a)).

    # 1.1: Work on the three downward clades rooted at sibling, child, and other. For
    # each of them, identify their neighbor (i.e., the other clade on the same side of
    # the target branch) before and after swapping.
    cdef Py_ssize_t[9] arr
    arr[0], arr[1], arr[2] = sibling, parent, other
    arr[3], arr[4], arr[5] = child, other, parent
    arr[6], arr[7], arr[8] = other, child, sibling
    for i in range(0, 9, 3):
        curr, before, after = arr[i], arr[i + 1], arr[i + 2]
        depth = tree[curr, 5]
        depth_2 = depth - 2
        ii = tree[curr, 6]
        for i in range(ii, ii + tree[curr, 4] * 2 - 1):
            a = preodr[i]
            npot = npots[tree[a, 5] - depth_2]
            jj = tree[a, 6]
            admx = &adm[a, 0]
            for j in range(jj + 1, jj + tree[a, 4] * 2 - 1):
                b = preodr[j]
                admx[b] = adm[b, a] = admx[b] - npot * (
                    adm[b, before] - adm[b, after]
                )

    # 1.2: Work on the upward clade rooted at the parent node. The process is similar
    # to that in `_bal_avgdist_insert`. Basically, we will consider three categories of
    # distances:
    # - previous ancestor vs. current ancestor
    # - previous & current ancestor vs. cousin or its descendant
    # - cousin or its descendant vs. descendant of the former
    before, after = sibling, child
    depth = tree[parent, 5]
    depth_2 = depth + 2

    # iterate over ancestors of the target node, starting from its parent
    anc_i = 0
    curr = parent
    while True:
        stack[anc_i] = curr
        depth_diff = depth_2 - 2 * tree[curr, 5] + 1

        # each previous ancestor vs. current ancestor
        diff = adm[curr, before] - adm[curr, after]
        for i in range(anc_i):
            a = stack[i]
            adm[curr, a] = adm[a, curr] = adm[curr, a] - npots[i + 2] * diff

        if not curr:
            break
        anc_i += 1

        # iterate over nodes within each cousin clade (cousin and its descendants)
        cousin = tree[curr, 3]
        ii = tree[cousin, 6]
        for i in range(ii, ii + tree[cousin, 4] * 2 - 1):
            a = preodr[i]
            admx = &adm[a, 0]

            # each previous & current ancestor vs. current node
            diff = adm[a, before] - adm[a, after]
            for j in range(anc_i):
                b = stack[j]
                admx[b] = adm[b, a] = admx[b] - npots[j + 2] * diff

            # current node vs. each descendant
            jj = tree[a, 6]
            npot = npots[depth_diff + tree[a, 5]]  #### trick
            for j in range(jj + 1, jj + tree[a, 4] * 2 - 1):
                b = preodr[j]
                admx[b] = adm[b, a] = admx[b] - npot * (
                    adm[b, before] - adm[b, after]
                )

        curr = tree[curr, 2]

    # Step 2: Update distances between subtrees within each side of the target branch
    # and the other end of the branch (A5.3(b)).
    start = tree[target, 6]
    end = start + tree[target, 4] * 2 - 1
    for i in range(start):
        node = preodr[i]
        adm[node, target] = adm[target, node] = 0.5 * (
            adm[node, sibling] + adm[node, other]
        )
    for i in range(start + 1, end):
        node = preodr[i]
        adm[node, target] = adm[target, node] = 0.5 * (
            adm[node, parent] + adm[node, child]
        )
    for i in range(end, n):
        node = preodr[i]
        adm[node, target] = adm[target, node] = 0.5 * (
            adm[node, sibling] + adm[node, other]
        )

    # Step 3: Update the distance between the lower and upper subtrees of the target
    # node itself (A5.3(c)).

    # This step is omitted as we will not fill the diagonal of the matrix (i.e., the
    # distance between the lower and upper subtree of the same node). See `_avgdist_
    # matrix`.

    # adm[target, target] = 0.25 * (
    #     adm[parent, sibling]
    #     + adm[parent, other]
    #     + adm[child, sibling]
    #     + adm[child, other]
    # )


cdef void _ols_swap(
    floating* L,
    Py_ssize_t* side,
    Py_ssize_t p_size,
    Py_ssize_t s_size,
    Py_ssize_t l_size,
    Py_ssize_t r_size,
    floating ad_before,
    floating ad_left,
    floating ad_right,
) noexcept nogil:
    r"""Calculate the change in overall tree length after a given swap.

    In a given quartet AB|CD, the swap takes place between subtrees B and C. The
    decrease of overall branch length, delta L (larger is better), can be calculated
    according to Eq. 9 of Desper and Gascuel (2002):

        delta L = ((lambda - 1)(d(A, C) + d(B, D)) - (lambda' - 1)(d(A, B) + d(C, D))
            - (lambda - lambda')(d(A, D) + d(B, C))) / 2

    Where:

        lambda  = (|A||D| + |B||C|) / ((|A| + |B|)(|C| + |D|))
        lambda' = (|A||D| + |B||C|) / ((|A| + |C|)(|B| + |D|))

    This function calculates delta L for two scenarios: 1) swap between left child and
    sibling, and 2) swap between right child and sibling.

    These two calculations share some steps, which are exploited by the following code.

    Factor 0.5 (/2) is omitted. The output value is negated, in order to match Python's
    minheap data structure (smallest popped first).

    Then the smaller change if negative is retained, according to A4.2.

    """
    cdef Py_ssize_t size_0 = p_size * s_size + l_size * r_size
    cdef Py_ssize_t size_1 = p_size * r_size + s_size * l_size
    cdef Py_ssize_t size_2 = p_size * l_size + s_size * r_size

    cdef floating ad_bl = ad_before - ad_left
    cdef floating ad_br = ad_before - ad_right

    cdef floating ad_lr_size = (ad_left - ad_right) / (size_1 + size_2)

    cdef floating L1 = size_1 * (ad_br / (size_0 + size_1) - ad_lr_size) - ad_bl
    cdef floating L2 = size_2 * (ad_bl / (size_0 + size_2) + ad_lr_size) - ad_br

    if L1 >= 0 and L2 >= 0:
        L[0] = 0
    elif L1 <= L2:
        L[0] = L1
        side[0] = 0
    else:
        L[0] = L2
        side[0] = 1


def _ols_all_swaps(
    floating[::1] lens,
    Py_ssize_t[:, ::1] tree,
    floating[:, ::1] adm,
):
    r"""Evaluate all possible swaps at all internal branches of a tree.

    Using an OLS framework.

    The results (length change, smaller is better) are saved to `lens`. The child side
    of each node that achieves this length change is saved to column 7 of the tree.

    Implemented according to A4.2(a) of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t node, left, right, parent, sibling
    
    # root is zero
    lens[0] = 0

    # calculate on internal branches (nodes with children)
    for node in range(1, tree.shape[0]):
        if tree[node, 0]:
            left = tree[node, 0]
            right = tree[node, 1]
            parent = tree[node, 2]
            sibling = tree[node, 3]

            # calculate length change
            _ols_swap(
                &lens[node],
                &tree[node, 7],
                m - tree[parent, 4],
                tree[sibling, 4],
                tree[left, 4],
                tree[right, 4],
                adm[parent, sibling] + adm[left, right],
                adm[parent, left] + adm[sibling, right],
                adm[parent, right] + adm[sibling, left],
            )
        
        # tips are zero
        else:
            lens[node] = 0


def _ols_corner_swaps(
    Py_ssize_t target,
    list heap,
    floating[::1] lens,
    Py_ssize_t[:, ::1] tree,
    floating[:, ::1] adm,
):
    r"""Update swaps of the four corner branches of a swapped branch.

    Using an OLS framework.

    Specifically, the four corner branches concern all nodes illustrated below.

                   |
                   x
                  / \
             parent  x
               /   \
           target  sibling
            /  \      / \
        left  right  x   x
         / \    / \
        x   x  x   x

    Implemented according to A4.2(c) of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t i, left, right, parent, sibling

    # update four corner branches if they are internal
    # 0: left, 1: right, 2: parent, 3: sibling
    # for 0, 1, 3, check left child; for 2, check parent
    for i in range(4):
        node = tree[target, i]
        if (node if i == 2 else tree[node, 0]):
            left = tree[node, 0]
            right = tree[node, 1]
            parent = tree[node, 2]
            sibling = tree[node, 3]

            # calculate length change
            _ols_swap(
                &lens[node],
                &tree[node, 7],
                m - tree[parent, 4],
                tree[sibling, 4],
                tree[left, 4],
                tree[right, 4],
                adm[parent, sibling] + adm[left, right],
                adm[parent, left] + adm[sibling, right],
                adm[parent, right] + adm[sibling, left],
            )

            # if length is reduced, push the swap into the heap
            if lens[node]:
                heappush(heap, (lens[node], node, tree[node, 7]))


def _bal_all_swaps(
    floating[::1] gains,
    Py_ssize_t[::1] sides,
    Py_ssize_t[::1] nodes,
    floating[:, ::1] adm,
    Py_ssize_t[:, ::1] tree,
):
    r"""Evaluate possible swaps at all internal branches of a tree.

    Using a balanced framework.

    Implemented according to Eq. 12 and A4.2(a) of Desper and Gascuel (2002).

        \delta L = ((d(A, B) + d(C, D)) - (d(A, C) + d(B, D))) / 4

    Larger and positive \delta L is favored.

    Factor 0.25 (/4) is omitted in the calculation.

    `lens` should be initiated with zeros, because root and tips are assigned zeros
    automatically and are skipped during the iteration.

    """
    cdef Py_ssize_t branch, node, left, right, parent, sibling
    cdef floating L1, L2, Lcomm
    for branch in range(nodes.shape[0]):
        node = nodes[branch]
        left = tree[node, 0]
        if not left:
            continue
        right, parent, sibling = tree[node, 1], tree[node, 2], tree[node, 3]
        Lcomm = adm[parent, sibling] + adm[left, right]
        L1 = adm[parent, left] + adm[sibling, right]
        L2 = adm[parent, right] + adm[sibling, left]
        if L1 >= Lcomm and L2 >= Lcomm:
            gains[branch] = 0
        elif L1 <= L2:
            gains[branch], sides[branch] = Lcomm - L1, 0
        else:
            gains[branch], sides[branch] = Lcomm - L2, 1
