# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True

from cython cimport floating
from cython.parallel cimport prange
from libc.string cimport memmove
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
    floating[:, ::1] adk,
    Py_ssize_t taxon,
    floating[:, :] dm,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] postodr,
):
    r"""Calculate balanced average distances between a new taxon and existing subtrees.

    This function resembles :func:`_avgdist_taxon` but uses the balanced framework.

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, left, right
    cdef Py_ssize_t n = 2 * tree[0, 4] - 1
    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]
    for i in range(n):
        node = postodr[i]
        left, right = tree[node, 0], tree[node, 1]
        if left == 0:
            adkl[node] = dm[taxon, right]
        else:
            adkl[node] = 0.5 * (adkl[left] + adkl[right])
    adku[0] = dm[taxon, 0]
    for i in range(1, n):
        node = preodr[i]
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


# def _bal_lengths(
#     floating[::1] lens,
#     floating[:, ::1] adm,
#     Py_ssize_t[:, ::1] tree,
# ):
#     r"""Calculate branch lengths of a tree based on the balanced framework.

#     Using a balanced average distance matrix between all pairs of subtrees.

#     This function resembles :func:`_ols_lengths` but it uses the balanced framework.

#     Implemented according to Eqs. 3 & 4 and the description on top of pg. 691 of Desper
#     and Gascuel (2002).

#     """
#     cdef Py_ssize_t node, left, right, parent, sibling

#     for node in range(1, 2 * tree[0, 4] - 1):
#         left = tree[node, 0]
#         parent = tree[node, 2]
#         sibling = tree[node, 3]
#         if left == 0:
#             lens[node] = 0.5 * (
#                 adm[parent, node] + adm[sibling, node] - adm[parent, sibling]
#             )
#         else:
#             right = tree[node, 1]
#             lens[node] = 0.25 * (
#                 adm[parent, left]
#                 + adm[sibling, right]
#                 + adm[parent, right]
#                 + adm[sibling, left]
#             ) - 0.5 * (adm[parent, sibling] + adm[left, right])
#     left, right = tree[0, 0], tree[0, 1]
#     # lens[0] = 0.5 * (adm[left, 0] + adm[right, 0] - adm[left, right])
#     lens[0] = 0.5 * (adm[0, left] + adm[0, right] - adm[left, right])


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
    cdef Py_ssize_t node, left, right, parent, sibling

    for node in range(1, 2 * tree[0, 4] - 1):
        left, parent, sibling = tree[node, 0], tree[node, 2], tree[node, 3]

        # node is tip
        if left == 0:
            if tree[parent, 0] == node:
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
            if tree[parent, 0] == node:
                right = tree[node, 1]
                lens[node] = 0.25 * (
                    adm[parent, left]
                    + adm[sibling, right]
                    + adm[parent, right]
                    + adm[sibling, left]
                ) - 0.5 * (adm[parent, sibling] + adm[right, left])

            # sibling is left, node is right
            else:
                lens[node] = 0.25 * (
                    adm[parent, left]
                    + adm[right, sibling]
                    + adm[parent, right]
                    + adm[left, sibling]
                ) - 0.5 * (adm[parent, sibling] + adm[right, left])

    left, right = tree[0, 0], tree[0, 1]
    lens[0] = 0.5 * (adm[0, left] + adm[0, right] - adm[right, left])


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
    floating[::1] lens,
    floating[:, ::1] adm,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
):
    """Find the branch with the minimum length change after inserting a new taxon.

    This function resembles :func:`_ols_min_branch_d2` but it 1) uses the balanced
    framework and 2) calculates based on the entire matrix. See also the note of the
    latter.

    Implemented according to Eq. 10 of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, parent, sibling
    cdef floating L

    cdef Py_ssize_t min_node = 0
    cdef floating min_len = 0

    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]

    #####
    cdef floating cell

    lens[min_node] = min_len
    for i in range(1, 2 * tree[0, 4] - 1):
        node = preodr[i]
        parent = tree[node, 2]
        sibling = tree[node, 3]

        # factor 0.25 is omitted
        # lens[node] = L = lens[parent] + (
        #     adm[sibling, parent] + adkl[node] - adm[sibling, node] - adku[parent]
        # )
        # lens[node] = L = lens[parent] + (
        #     adm[parent, sibling] + adkl[node] - adm[sibling, node] - adku[parent]
        # )

        #####
        if tree[parent, 0] == node:
            cell = adm[sibling, node]  # node is left, sibling is right
        else:
            cell = adm[node, sibling]
        lens[node] = L = lens[parent] + (
            adm[parent, sibling] + adkl[node] - cell - adku[parent]
        )

        if L < min_len:
            min_len, min_node = L, node
    return min_node


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
#     floating[::1] powers,
#     Py_ssize_t[::1] stack,
#     Py_ssize_t[::1] paths,
#     Py_ssize_t[::1] gens,
# ):
#     r"""Update balanced average distance matrix after taxon insertion.

#     This function resembles :func:`_avgdist_d2_insert` but it 1) uses the balanced
#     framework and 2) updates the entire matrix. The latter makes it the dominant term
#     of the entire algorithm.

#     Two additional parameters are provided: `powers` is a pre-calculated array of
#     2^(-l) powers (l is the depth difference between two nodes). `stack` is an
#     integer array to store ancestral nodes of target.

#     """
#     cdef Py_ssize_t i, j, ii, jj, anc_i
#     cdef Py_ssize_t parent, sibling, depth
#     cdef Py_ssize_t curr, anc, cousin, depoff, depth_diff
#     cdef Py_ssize_t a, b
#     cdef floating power, diff

#     # dimensions and positions
#     cdef Py_ssize_t n = 2 * tree[0, 4] - 1
#     cdef Py_ssize_t link = n
#     cdef Py_ssize_t tip = n + 1

#     cdef floating* adkl = &adk[0, 0]
#     cdef floating* adku = &adk[1, 0]

#     cdef floating* powers_2 = &powers[2]

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
#             power = powers[tree[a, 5] + 1]
#             for i in range(ii - tree[a, 4] * 2 + 2, ii):
#                 b = postodr[i]
#                 adm[a, b] = adm[b, a] = adm[a, b] + power * (adkl[b] - adm[0, b])

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
#         power = powers[tree[a, 5] - depth + 1]
#         for j in range(jj - tree[a, 4] * 2 + 2, jj):
#             b = postodr[j]
#             adm[a, b] = adm[b, a] = adm[a, b] + power * (adkl[b] - adm[target, b])

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
#             adm[anc, a] = adm[a, anc] = adm[anc, a] + powers[i + 2] * diff

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
#                 adm[a, b] = adm[b, a] = adm[a, b] + powers[j + 2] * diff

#             # Iterate over descendants of each member of the clade, and calculate the
#             # distance between the former (upper, containing k) and the latter (lower).
#             jj = tree[a, 7]
#             power = powers[depth_diff + tree[a, 5]]
#             for j in range(jj - tree[a, 4] * 2 + 2, jj):
#                 b = postodr[j]
#                 adm[a, b] = adm[b, a] = adm[a, b] + power * (adkl[b] - adm[b, target])

#         curr = anc
#         anc_i += 1


def _bal_avgdist_insert(
    floating[:, ::1] adm,
    Py_ssize_t target,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] postodr,
    floating[::1] powers,
    Py_ssize_t[::1] stack,
    Py_ssize_t[::1] paths,
    Py_ssize_t[::1] gens,
):
    r"""Update balanced average distance matrix after taxon insertion.

    This function resembles :func:`_avgdist_d2_insert` but it 1) uses the balanced
    framework and, 2) updates the entire matrix (`adm`), which makes it the dominant
    term of the entire algorithm.

    This function is the serial version. A parallel version can be found in
    :func:`_bal_avgdist_insert_p`.

    This function assumes that the taxon will be inserted into the branch connecting
    the target node and its parent. After insertion, the taxon will become the right
    sibling of the target.

                                    parent
             parent                  /  \
             /  \        =>       link  sibling
        target  sibling           /  \
                             target  taxon (k)

    This function only updates distance [i, j] but not [j, i], saving half of the
    assignment operations, which are expensive given the size of the matrix. Here,
    i and j suffice i < j in a postorder traversal (left - right - parent), in which
    ancestors are always latter than their descendants, right children (and their
    descendants) are always latter than left children (and their descendants).
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
    navigating all ancestors of target up until root.

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

    Two additional parameters are provided: `powers` is a pre-calculated array of
    2^(-l) powers (l is the depth difference between two nodes). `stack` is an
    integer array to store ancestral nodes of target.

    """
    cdef Py_ssize_t i, j, ii, jj, anc_i
    cdef Py_ssize_t parent, sibling, depth
    cdef Py_ssize_t curr, anc, cousin, depoff
    cdef Py_ssize_t a, b

    # Intermediate variables
    cdef floating cell  # specific value in `adm` (a matrix)
    cdef floating diff  # difference in distance; usually from `adk`
    cdef floating power  # negative power of 2; from `powers`

    # Dimensions and positions
    cdef Py_ssize_t n = 2 * tree[0, 4] - 1
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    # Pointers to lower and upper tree distances
    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]

    # Pointers to specific rows in `adm`
    cdef floating* adm_t = &adm[target, 0]
    cdef floating* adm_l = &adm[link, 0]
    cdef floating* adm_k = &adm[tip, 0]
    cdef floating* adm_c  # each ancestor of target
    cdef floating* adm_a  # each node within a clade

    cdef floating* powers_2 = &powers[2]

    ###### Special case: insert into the root branch. ######

    #    root (target)      root
    #    /  \               /  \
    # left  right   =>   link  taxon (k)
    #                    /  \
    #                 left  right

    if target == 0:
        # Transfer distance between k and root (upper).
        adm_t[tip]= adku[0]

        # k to link: equals to k to root (lower).
        adm_k[link] = adkl[0]

        # Root to link: de novo calculation according to the equation in A4.1(b). It is
        # basically the distance between the upper and lower subtrees of the root.
        adm_t[link] = 0.5 * (adm_t[tree[0, 0]] + adm_t[tree[0, 1]])

        # Iterate over all nodes but root (last in postorder).
        for i in range(n - 1):
            a = postodr[i]
            adm_a = &adm[a, 0]

            # Transfer distances between the node (lower) and k.
            adm_k[a] = diff = adkl[a]
            cell = adm_t[a]

            # Calculate the distance between the node (lower) and link (upper, with two
            # taxa (0 and k) added) using Eq. 8.
            adm_l[a] = 0.5 * (diff + cell)

            # Calculate this intermediate once, such that it won't need to be repeatedly
            # calculated in the subsequent loop. It is fine to store this value in `adk`,
            # which will be entirely re-calculated by `_bal_avgdist_taxon`.
            adkl[a] = diff - cell

            # Calculate the distances between the node (U, containing k) and each of its
            # descendants (L). We can now reuse the intermediate calculated above. It is
            # important to do this in postorder.
            ii = tree[a, 7]
            power = powers[tree[a, 5] + 1]
            for i in range(ii - tree[a, 4] * 2 + 2, ii):
                b = postodr[i]
                adm_a[b] += power * adkl[b]

        return

    ###### Regular case: insert into any other branch. ######

    parent, sibling, depth = tree[target, 2], tree[target, 3], tree[target, 5]

    ### Step 1: Distances around the insertion point. ###

    # Distance between k (lower) and link (upper) equals to that between k and the
    # upper subtree of target.
    adm_l[tip] = adku[target]

    # Distance between target (lower) and link (upper) needs to be calculated using the
    # equation in A4.1(c). Basically, it is the distance between the lower and upper
    # subtrees of the same target.
    if tree[parent, 0] == target:
        cell = adm[sibling, target]
    else:
        cell = adm[target, sibling]
    adm_l[target] = 0.5 * (cell + adm[parent, target])

    # Transfer pre-calculated distance between target (lower) and k (lower).
    adm_k[target] = adkl[target]

    ### Step 2: Distances within the clade below target. ###

    # Locate the clade below target (excluding target).
    depoff = 1 - depth
    ii = tree[target, 7]
    for i in range(ii - tree[target, 4] * 2 + 2, ii):
        a = postodr[i]
        adm_a = &adm[a, 0]

        # Transfer pre-calculated distance between k (lower) and any node within
        # the clade (lower).
        adm_k[a] = diff = adkl[a]

        # Distance from any descendant (lower) to link (upper) equals to that to
        # target.
        adm_l[a] = cell = adm_t[a]

        # Calculate the distance between each node within the clade (lower) and target
        # (upper).
        adm_t[a] = 0.5 * (diff + cell)

        # intermediate
        adkl[a] = diff - cell

        # Within the clade, find all ancestor (a) - descendant (b) pairs, and
        # calculate the distance between the upper subtree of a (with k) and
        # the lower subtree of b.
        power = powers[tree[a, 5] + depoff]
        jj = tree[a, 7]
        for j in range(jj - tree[a, 4] * 2 + 2, jj):
            b = postodr[j]
            adm_a[b] += power * adkl[b]

    ### Step 3: Distances among nodes outside the clade. ###

    # Iterate over ancestors of target in ascending order till root.
    anc_i = 0
    curr = target
    depoff = 2 - depth
    while curr:
        stack[anc_i] = anc = tree[curr, 2]
        adm_c = &adm[anc, 0]
        diff, cell = adku[anc], adm_c[target]

        # Transfer the pre-calculated distance between k and ancestor (upper).
        adm_c[tip] = diff

        # Calculate the distance between link (lower, with k) and ancestor (upper).
        adm_c[link] = 0.5 * (diff + cell)

        # Calculate the distance between each previous ancestor (lower, with k)
        # and the current ancestor (upper).
        diff -= cell
        for i in range(anc_i):
            adm_c[stack[i]] += powers_2[i] * diff

        # Locate the cousin clade descending from the ancestor (including cousin).
        cousin = tree[curr, 3]
        ii = tree[cousin, 7] + 1

        # Determine whether cousin is the right or left child of the shared ancestor.
        # This helps to determine the order of coordinates of each distance to be
        # updated. If cousin is right, then coordinates are [cousin, curr].
        if tree[anc, 0] == curr:
            for i in range(ii - tree[cousin, 4] * 2 + 1, ii):
                a = postodr[i]
                adm_a = &adm[a, 0]
                diff, cell = adkl[a], adm_a[target]

                # Transfer the pre-calculated distances between k and each descendant
                # (lower).
                adm_a[tip] = adkl[a]

                # Calculate the distance between link (lower, with k) and each
                # descendant (lower).
                adm_a[link] = 0.5 * (diff + cell)

                # intermediate
                adkl[a] = diff = diff - cell

                # Calculate the distance between each previous ancestor (lower, with k
                # and each descendant (lower).
                for j in range(anc_i):
                    adm_a[stack[j]] += powers_2[j] * diff

                # Iterate over descendants of each member of the clade, and calculate the
                # distance between the former (upper, with k) and the latter (lower).
                jj = tree[a, 7]
                power = powers[tree[a, 5] + depoff]
                for j in range(jj - tree[a, 4] * 2 + 2, jj):
                    b = postodr[j]
                    adm_a[b] += power * adkl[b]

        # If cousin is left, coordinates are [curr, cousin]. Everything else is the
        # same, but code might be slightly slower due to the coordinates.
        else:
            for i in range(ii - tree[cousin, 4] * 2 + 1, ii):
                a = postodr[i]
                adm_a = &adm[a, 0]
                diff, cell = adkl[a], adm_t[a]
                adm_k[a] = diff
                adm_l[a] = 0.5 * (diff + cell)
                adkl[a] = diff = diff - cell
                for j in range(anc_i):
                    adm[stack[j], a] += powers_2[j] * diff
                jj = tree[a, 7]
                power = powers[tree[a, 5] + depoff]
                for j in range(jj - tree[a, 4] * 2 + 2, jj):
                    b = postodr[j]
                    adm_a[b] += power * adkl[b]

        curr = anc
        anc_i += 1
        depoff += 2


### parallelize as it goes; no intermediate ###

def _bal_avgdist_insert_p(
    floating[:, ::1] adm,
    Py_ssize_t target,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] postodr,
    floating[::1] powers,
    Py_ssize_t[::1] stack,
    Py_ssize_t[::1] paths,
    Py_ssize_t[::1] gens,
    int chunksize = 10,
    int minclade = 100,
    bint adaptive = False,
):
    r"""Update balanced average distance matrix after taxon insertion."""
    cdef Py_ssize_t i, j, ii, jj, anc_i
    cdef Py_ssize_t parent, sibling, depth
    cdef Py_ssize_t curr, anc, cousin, depoff
    cdef Py_ssize_t a, b

    cdef floating cell, diff, power

    # Parallelization
    cdef int ops     # total operations
    cdef int chunk   # chunk size
    cdef bint worth  # whether use threads

    cdef Py_ssize_t n = 2 * tree[0, 4] - 1
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]

    cdef floating* adm_t = &adm[target, 0]
    cdef floating* adm_l = &adm[link, 0]
    cdef floating* adm_k = &adm[tip, 0]
    cdef floating* adm_c
    cdef floating* adm_a

    cdef floating* powers_2 = &powers[2]

    ###### Special case: insert into the root branch. ######

    if target == 0:
        adm_t[tip]= adku[0]
        adm_k[link] = adkl[0]
        adm_t[link] = 0.5 * (adm_t[tree[0, 0]] + adm_t[tree[0, 1]])

        # Iterate over all nodes but the root.
        # NOTE: Nodes are ordered by the order they were inserted into the tree. Thus,
        # the sizes of clades do not have an obvious ascending or descending pattern.
        # This eases workload distribution across threads.
        chunk, worth = config_prange(n - 1, chunksize, minclade, adaptive)
        for a in prange(
            1, n, nogil=True, schedule="dynamic", chunksize=chunk, use_threads_if=worth
        ):
            adm_a = &adm[a, 0]
            adm_k[a] = adkl[a]
            adm_l[a] = 0.5 * (adkl[a] + adm_t[a])
            ii = tree[a, 7]
            power = powers[tree[a, 5] + 1]
            for i in range(ii - tree[a, 4] * 2 + 2, ii):
                b = postodr[i]
                adm_a[b] += power * (adkl[b] - adm_t[b])

        return

    ###### Regular case: insert into any other branch. ######

    parent, sibling, depth = tree[target, 2], tree[target, 3], tree[target, 5]

    ### Step 1: Distances around the insertion point. ###

    adm_l[tip] = adku[target]
    cell = adm[sibling, target] if tree[parent, 0] == target else adm_t[sibling]
    adm_l[target] = 0.5 * (cell + adm[parent, target])
    adm_k[target] = adkl[target]

    ### Step 2: Distances within the clade below target. ###

    depoff = 1 - depth
    ops = tree[target, 4] * 2 - 2
    if ops > 0:
        ii = tree[target, 7]
        chunk, worth = config_prange(ops, chunksize, minclade, adaptive)
        for i in prange(
            ii - ops, ii, nogil=True, schedule="dynamic", chunksize=chunk,
            use_threads_if=worth
        ):
            a = postodr[i]
            adm_a = &adm[a, 0]
            adm[tip, a] = adkl[a]
            adm[link, a] = adm[target, a]
            jj = tree[a, 7]
            power = powers[tree[a, 5] + depoff]
            for j in range(jj - tree[a, 4] * 2 + 2, jj):
                b = postodr[j]
                adm_a[b] += power * (adkl[b] - adm_t[b])

        # NOTE: This loop cannot be merged into the above loop as the latter has been
        # parallelized (i.e., no longer in postorder).
        # TODO: Parallelize this loop
        for i in range(ii - ops, ii):
            a = postodr[i]
            adm_t[a] = 0.5 * (adm_t[a] + adkl[a])

    ### Step 3: Distances among nodes outside the clade. ###

    anc_i = 0
    curr = target
    depoff = 2 - depth
    while curr:
        stack[anc_i] = anc = tree[curr, 2]
        adm_c = &adm[anc, 0]
        adm_c[tip] = adku[anc]
        adm_c[link] = 0.5 * (adku[anc] + adm_c[target])
        diff = adku[anc] - adm_c[target]
        for i in range(anc_i):
            adm_c[stack[i]] += powers_2[i] * diff

        cousin = tree[curr, 3]
        ii = tree[cousin, 7] + 1
        ops = tree[cousin, 4] * 2 - 1
        chunk, worth = config_prange(ops, chunksize, minclade, adaptive)

        # Cousin is right
        if tree[anc, 0] == curr:
            for i in prange(
                ii - ops, ii, nogil=True, schedule="dynamic", chunksize=chunk,
                use_threads_if=worth
            ):
                a = postodr[i]
                adm_a = &adm[a, 0]
                adm_a[tip] = adkl[a]
                adm_a[link] = 0.5 * (adkl[a] + adm_a[target])
                diff = adkl[a] - adm_a[target]
                for j in range(anc_i):
                    adm_a[stack[j]] += powers_2[j] * diff
                jj = tree[a, 7]
                power = powers[tree[a, 5] + depoff]
                for j in range(jj - tree[a, 4] * 2 + 2, jj):
                    b = postodr[j]
                    adm_a[b] += power * (adkl[b] - adm[b, target])

        # Cousin is left
        else:
            for i in prange(
                ii - ops, ii, nogil=True, schedule="dynamic", chunksize=chunk,
                use_threads_if=worth
            ):
                a = postodr[i]
                adm_a = &adm[a, 0]
                adm_k[a] = adkl[a]
                adm_l[a] = 0.5 * (adkl[a] + adm_t[a])
                diff = adkl[a] - adm_t[a]
                for j in range(anc_i):
                    adm[stack[j], a] += powers_2[j] * diff
                jj = tree[a, 7]
                power = powers[tree[a, 5] + depoff]
                for j in range(jj - tree[a, 4] * 2 + 2, jj):
                    b = postodr[j]
                    adm_a[b] += power * (adkl[b] - adm_t[b])

        curr = anc
        anc_i += 1
        depoff += 2


### centralized parallelization : separate vertical fill ###

def _bal_avgdist_insert_p2(
    floating[:, ::1] adm,
    Py_ssize_t target,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] postodr,
    floating[::1] powers,
    Py_ssize_t[::1] stack,
    Py_ssize_t[::1] paths,
    Py_ssize_t[::1] gens,
    int chunksize = 10,
    int minclade = 100,
    bint adaptive = False,
):
    r"""Update balanced average distance matrix after taxon insertion."""
    cdef Py_ssize_t i, j, ii, jj, anc_i
    cdef Py_ssize_t parent, sibling, size, depth
    cdef Py_ssize_t curr, anc, cousin, depoff
    cdef Py_ssize_t a, b

    cdef floating cell, diff, power

    # Parallelization
    cdef int ops     # total operations
    cdef int chunk   # chunk size
    cdef bint worth  # whether use threads

    cdef int minops = 5000

    # dimensions and positions
    cdef Py_ssize_t n = 2 * tree[0, 4] - 1
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]

    cdef floating* adm_t = &adm[target, 0]
    cdef floating* adm_l = &adm[link, 0]
    cdef floating* adm_k = &adm[tip, 0]
    cdef floating* adm_c
    cdef floating* adm_a

    cdef floating* powers_2 = &powers[2]

    paths[target] = 0

    ###### Special case: insert into the root branch. ######

    if target == 0:
        adm_t[tip] = adku[0]
        adm_k[link] = adkl[0]
        adm_t[link] = 0.5 * (adm_t[tree[0, 0]] + adm_t[tree[0, 1]])
        for a in prange(1, n, nogil=True, use_threads_if=n >= minops):
            adm_k[a] = diff = adkl[a]
            cell = adm_t[a]
            adm_l[a] = 0.5 * (diff + cell)
            adkl[a] = diff - cell
            paths[a] = tree[a, 5] + 1

    ###### Regular case: insert into any other branch. ######

    else:
        parent, sibling, depth = tree[target, 2], tree[target, 3], tree[target, 5]

        ### Step 1: Distances around the insertion point. ###

        adm_l[tip] = adku[target]
        cell = adm[sibling, target] if tree[parent, 0] == target else adm_t[sibling]
        adm_l[target] = 0.5 * (cell + adm[parent, target])
        adm_k[target] = adkl[target]

        ### Step 2: Distances within the clade below target. ###

        # Locate the clade below target (excluding target). Skip if target is a tip.
        # ops = 12
        depoff = 1 - depth
        ops = tree[target, 4] * 2 - 2
        ii = tree[target, 7]
        for i in prange(ii - ops, ii, nogil=True, use_threads_if=ops > minops):
            a = postodr[i]
            adm_k[a] = diff = adkl[a]
            adm_l[a] = cell = adm_t[a]
            adm_t[a] = 0.5 * (diff + cell)
            adkl[a] = diff - cell
            paths[a] = tree[a, 5] + depoff

        ### Step 3: Distances among nodes outside the clade. ###

        # Iterate over ancestors of target in ascending order.
        anc_i = 0
        curr = target
        depoff = 2 - depth
        while curr:
            stack[anc_i] = anc = tree[curr, 2]
            paths[anc] = 0
            adm_c = &adm[anc, 0]
            adm_c[tip] = diff = adku[anc]
            cell = adm_c[target]
            adm_c[link] = 0.5 * (diff + cell)
            diff -= cell
            for i in range(anc_i):
                adm_c[stack[i]] += powers_2[i] * diff

            cousin = tree[curr, 3]
            ii = tree[cousin, 7]
            ops = tree[cousin, 4] * 2 - 1

            # Cousin is right
            # ops = 13 + 5n
            if tree[anc, 0] == curr:
                for i in prange(ii - ops + 1, ii + 1, nogil=True, use_threads_if=ops * anc_i > minops):
                    a = postodr[i]
                    adm_a = &adm[a, 0]
                    diff, cell = adkl[a], adm_a[target]
                    adm_a[tip] = diff
                    adm_a[link] = 0.5 * (diff + cell)
                    adkl[a] = diff = diff - cell
                    for j in range(anc_i):
                        adm_a[stack[j]] += powers_2[j] * diff
                    paths[a] = tree[a, 5] + depoff

            # Cousin is left
            else:
                for i in prange(ii - ops + 1, ii + 1, nogil=True, use_threads_if=ops * anc_i > minops):
                    a = postodr[i]
                    adm_a = &adm[a, 0]
                    diff, cell = adkl[a], adm_t[a]
                    adm_k[a] = diff
                    adm_l[a] = 0.5 * (diff + cell)
                    adkl[a] = diff = diff - cell
                    for j in range(anc_i):
                        adm[stack[j], a] += powers_2[j] * diff
                    paths[a] = tree[a, 5] + depoff

            curr = anc
            anc_i += 1
            depoff += 2

    ###### Fill all ancestor - descendant pairs. ######

    # chunk, worth = config_prange(n, chunksize, minclade, adaptive)
    # for a in prange(
    #     n, nogil=True, schedule="dynamic", chunksize=chunk, use_threads_if=worth
    # ):
    #     path = paths[a]
    #     size = tree[a, 4]
    #     if path > 0 and size > 0:
    #         power = powers[path]
    #         adm_a = &adm[a, 0]
    #         jj = tree[a, 7]
    #         for j in range(jj - size * 2 + 2, jj):
    #             b = postodr[j]
    #             adm_a[b] += power * adkl[b]


def _bal_avgdist_insert_p3(
    floating[:, ::1] adm,
    Py_ssize_t target,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] postodr,
    floating[::1] powers,
    Py_ssize_t[::1] stack,
    Py_ssize_t[::1] paths,
    Py_ssize_t[::1] gens,
    int chunksize = 10,
    int minclade = 100,
    bint adaptive = False,
):
    r"""Update balanced average distance matrix after taxon insertion.

    This function is the parallel version of :func:`_bal_avgdist_insert`.

    Within each subtree of n taxa, the total number of ancestor-descendant pairs is
    between O(nlogn) (balanced tree, best case) and O(n^2) (skewed tree, worst case).

    In postorder traversal (left - right - parent), ancestors are always latter than
    descendants, right children (and their descendants) are always latter than left
    children (and their descendants).

    """
    cdef Py_ssize_t i, j, ii, jj, anc_i
    cdef Py_ssize_t parent, sibling, size, depth
    cdef Py_ssize_t curr, anc, cousin, depoff
    cdef Py_ssize_t a, b

    # Intermediate variables:
    cdef floating cell  # specific value in `adm` (a matrix)
    cdef floating diff  # difference in distance; usually from `adk`
    cdef floating power  # negative power of 2; from `powers`

    # Parallelization
    cdef int ops     # total operations
    cdef int chunk   # chunk size
    cdef bint worth  # whether use threads

    cdef int minops = 5000

    # dimensions and positions
    cdef Py_ssize_t n = 2 * tree[0, 4] - 1
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    # Pointers to lower and upper tree distances
    cdef floating* adkl = &adk[0, 0]
    cdef floating* adku = &adk[1, 0]

    # Pointers to specific rows in `adm`
    cdef floating* adm_t = &adm[target, 0]
    cdef floating* adm_l = &adm[link, 0]
    cdef floating* adm_k = &adm[tip, 0]
    cdef floating* adm_c  # each ancestor of target
    cdef floating* adm_a  # each node within a clade

    cdef floating* powers_2 = &powers[2]

    paths[target] = 0
    gens[target] = 0

    ###### Special case: insert into the root branch. ######

    if target == 0:
        # Transfer distance between k and root (upper).
        adm_t[tip] = adku[0]

        # k to link: equals to k to root (lower).
        adm_k[link] = adkl[0]

        # Root to link: de novo calculation according to the equation in A4.1(b). It is
        # basically the distance between the upper and lower subtrees of the root.
        adm_t[link] = 0.5 * (adm_t[tree[0, 0]] + adm_t[tree[0, 1]])

        # Iterate over all nodes but the root.
        # NOTE: Nodes are ordered by the order they were inserted into the tree. Thus,
        # the sizes of clades do not have an obvious ascending or descending pattern.
        # This eases workload distribution across threads.
        # ops = 10
        # for a in prange(1, n, nogil=True, use_threads_if=n >= minops):
        for a in range(1, n):

            # Transfer distances between the node (lower) and k.
            adm_k[a] = diff = adkl[a]
            cell = adm_t[a]

            # Calculate the distance between the node (lower) and link (upper, with two
            # taxa (0 and k) added) using Eq. 8.
            adm_l[a] = 0.5 * (diff + cell)

            # Cache intermediate.
            adkl[a] = diff - cell

            # Calculate the path length between the node and k, which directly descends
            # from the root.
            paths[a] = tree[a, 5] + 1
            gens[a] = 0

    ###### Regular case: insert into any other branch. ######

    else:
        parent, sibling, depth = tree[target, 2], tree[target, 3], tree[target, 5]

        ### Step 1: Distances around the insertion point. ###

        # Distance between k (lower) and link (upper) equals to that between k and the
        # upper subtree of target.
        adm_l[tip] = adku[target]

        # Distance between target (lower) and link (upper) needs to be calculated using the
        # equation in A4.1(c). Basically, it is the distance between the lower and upper
        # subtrees of the same target.
        if tree[parent, 0] == target:
            cell = adm[sibling, target]
        else:
            cell = adm[target, sibling]
        adm_l[target] = 0.5 * (cell + adm[parent, target])

        # Transfer pre-calculated distance between target (lower) and k (lower).
        adm_k[target] = adkl[target]

        ### Step 2: Distances within the clade below target. ###

        # Locate the clade below target (excluding target). Skip if target is a tip.
        # ops = 12
        depoff = 1 - depth
        ops = tree[target, 4] * 2 - 2
        ii = tree[target, 7]
        for i in range(ii - ops, ii):
            a = postodr[i]

            # Transfer pre-calculated distance between k (lower) and any node within
            # the clade (lower).
            adm_k[a] = diff = adkl[a]

            # Distance from any descendant (lower) to link (upper) equals to that to
            # target.
            adm_l[a] = cell = adm_t[a]

            # Calculate the distance between node (lower) and target (upper).
            adm_t[a] = 0.5 * (diff + cell)

            adkl[a] = diff - cell

            # Path length
            paths[a] = tree[a, 5] + depoff
            gens[a] = 0

        ### Step 3: Distances among nodes outside the clade. ###

        # Iterate over ancestors of target in ascending order.
        anc_i = 0
        curr = target
        depoff = 2 - depth
        while curr:
            stack[anc_i] = anc = tree[curr, 2]
            paths[anc] = 0
            gens[anc] = 0
            adm_c = &adm[anc, 0]

            # Transfer the pre-calculated distance between k and the ancestor (upper).
            adm_c[tip] = diff = adku[anc]

            # Calculate the distance between link (lower, containing k) and the ancestor
            # (upper).
            cell = adm_c[target]
            adm_c[link] = 0.5 * (diff + cell)

            # Calculate the distance between each previous ancestor (lower, containing k)
            # and the current ancestor (upper).
            diff = diff - cell
            for i in range(anc_i):
                adm_c[stack[i]] += powers_2[i] * diff

            # Locate the cousin clade descending from the ancestor (including cousin).
            cousin = tree[curr, 3]
            ii = tree[cousin, 7]
            ops = tree[cousin, 4] * 2 - 1

            # Cousin is right, curr is left.
            # ops = 13 + 5n
            if tree[anc, 0] == curr:
                for i in range(ii - ops + 1, ii + 1):
                    a = postodr[i]
                    adm_a = &adm[a, 0]
                    diff, cell = adkl[a], adm_a[target]

                    # Transfer the pre-calculated distances between k and each descendant
                    adm_a[tip] = diff

                    # Calculate the distance between link (lower, containing k) and each
                    # descendant (lower).
                    adm_a[link] = 0.5 * (diff + cell)

                    adkl[a] = diff - cell
                    paths[a] = tree[a, 5] + depoff
                    gens[a] = anc_i

            # Cousin is left, curr is left.
            else:
                for i in range(ii - ops + 1, ii + 1):
                    a = postodr[i]
                    adm_a = &adm[a, 0]
                    diff, cell = adkl[a], adm_t[a]
                    adm_k[a] = diff
                    adm_l[a] = 0.5 * (diff + cell)
                    adkl[a] = diff - cell
                    paths[a] = tree[a, 5] + depoff
                    gens[a] = anc_i

            curr = anc
            anc_i += 1
            depoff += 2


def _fill_horizontal(
    Py_ssize_t target,
    floating[:, ::1] adm,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    floating[::1] powers,
    Py_ssize_t[::1] stack,
    Py_ssize_t[::1] gens,
    int chunksize = 10,
    int minclade = 100,
    bint adaptive = False,
):
    cdef Py_ssize_t a
    cdef int chunk   # chunk size
    cdef bint worth  # whether use threads
    cdef Py_ssize_t n = 2 * tree[0, 4] - 1
    cdef floating* adkl = &adk[0, 0]
    cdef floating* powers_2 = &powers[2]
    cdef Py_ssize_t post_t = tree[target, 7]
    cdef floating* adm_a
    cdef Py_ssize_t gen
    cdef floating diff
    chunk, worth = config_prange(n, chunksize, minclade, adaptive)
    for a in prange(
        n, nogil=True, schedule="dynamic", chunksize=chunk, use_threads_if=worth
    ):
        gen = gens[a]
        if gen > 0:
            diff = adkl[a]
            if tree[a, 7] > post_t:  # postorder
                adm_a = &adm[a, 0]
                for i in range(gen):
                    adm_a[stack[i]] += powers_2[i] * diff
            else:
                for i in range(gen):
                    adm[stack[i], a] += powers_2[i] * diff


def _fill_vertical(
    floating[:, ::1] adm,
    floating[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    floating[::1] powers,
    Py_ssize_t[::1] paths,
    int chunksize = 10,
    int minclade = 100,
    bint adaptive = False,
):
    cdef Py_ssize_t a, b, j, jj, path, size
    cdef floating power
    cdef int chunk   # chunk size
    cdef bint worth  # whether use threads
    cdef Py_ssize_t n = 2 * tree[0, 4] - 1
    cdef floating* adkl = &adk[0, 0]
    chunk, worth = config_prange(n, chunksize, minclade, adaptive)
    for a in prange(
        n, nogil=True, schedule="dynamic", chunksize=chunk, use_threads_if=worth
    ):
        path = paths[a]
        size = tree[a, 4]
        if path > 0 and size > 0:
            power = powers[path]
            adm_a = &adm[a, 0]

            # preorder
            jj = tree[a, 6]
            for j in range(jj + 1, jj + size * 2 - 1):
                b = preodr[j]
                adm_a[b] += power * adkl[b]

            # postorder
            # jj = tree[a, 7]
            # for j in range(jj - size * 2 + 2, jj):
            #     b = postodr[j]
            #     adm_a[b] += power * adkl[b]


def _insert_taxon(
    Py_ssize_t taxon,
    Py_ssize_t target,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] postodr,
    bint use_depth=True,
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
    # parts, the reality according to my tests is that cell-by-cell Cython code is
    # significantly faster than NumPy, and greatly reduces the overall runtime of the
    # entire algorithms. This effect is more obvious when the dataset is small, but
    # less so when it is large (but it is still there).
    #
    # The reason might be that when moving a block of elements within the same array,
    # NumPy needs to create a temporary array, but Cython can do the job in place.
    #
    # There might be a chance to re-consider NumPy (or even CuPy) API in the future.
    cdef Py_ssize_t left, right, parent, sibling, size, depth, pre_i, post_i
    cdef Py_ssize_t i, k, side, pre_i_after, curr

    # determine tree dimensions
    # typically n = 2 * taxon - 3, but this function doesn't enforce this
    cdef Py_ssize_t m = tree[0, 4]
    cdef Py_ssize_t n = m * 2 - 1
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    cdef Py_ssize_t* node

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
        node[4] = m + 1
        node[7] = n + 1

        # link
        node = &tree[link, 0]
        node[0] = left
        node[1] = right
        node[2] = 0
        node[3] = tip
        node[4] = m
        node[5] = 1
        node[6] = 1
        node[7] = n - 1

        # tip
        node = &tree[tip, 0]
        node[0] = 0
        node[1] = taxon
        node[2] = 0
        node[3] = link
        node[4] = 1
        node[5] = 1
        node[6] = n + 1
        node[7] = n

        # entire tree depth + 1
        if use_depth:
            for i in range(1, n):
                tree[i, 5] += 1

        # preorder
        for i in range(n - 1, 0, -1):
            tree[i, 6] += 1
            preodr[i + 1] = preodr[i]
        preodr[1] = link
        preodr[n + 1] = tip

        # postorder
        postodr[n - 1] = link
        postodr[n] = tip
        postodr[n + 1] = 0

    # Regular case (any other branch): The link becomes the parent of the target node,
    # and child of its original parent. Taxon k becomes the sibling
    else:
        node = &tree[target, 0]
        left = node[0]
        right = node[1]
        parent = node[2]
        sibling = node[3]
        size = node[4]
        depth = node[5]
        pre_i = node[6]
        post_i = node[7]

        side = int(tree[parent, 0] != target)
        tree[parent, side] = link
        tree[sibling, 3] = link
        node[2] = link
        node[3] = tip

        # preorder index of node after clade
        pre_i_after = pre_i + size * 2 - 1

        # link
        node = &tree[link, 0]
        node[0] = target
        node[1] = tip
        node[2] = parent
        node[3] = sibling
        node[4] = size + 1
        node[5] = depth
        node[6] = pre_i
        node[7] = post_i + 2

        # tip
        node = &tree[tip, 0]
        node[0] = 0
        node[1] = taxon
        node[2] = link
        node[3] = target
        node[4] = 1
        node[5] = depth + 1
        node[6] = pre_i_after + 1
        node[7] = post_i + 1

        # clade depth +1
        if use_depth:
            for i in range(pre_i, pre_i_after):
                tree[preodr[i], 5] += 1

        # preorder shift: nodes after clade +2, tip inserted after clade, nodes within
        # clade +1, link inserted before clade
        for i in range(n - 1, pre_i_after - 1, -1):
            tree[preodr[i], 6] += 2
            # k = preodr[i]
            # tree[k, 6] += 2
            # preodr[i + 2] = k

        memmove(&preodr[pre_i_after + 2], &preodr[pre_i_after], <size_t>(
            (n - pre_i_after) * sizeof(Py_ssize_t)
        ))

        preodr[pre_i_after + 1] = tip

        for i in range(pre_i_after - 1, pre_i - 1, -1):
            tree[preodr[i], 6] += 1
            # k = preodr[i]
            # tree[k, 6] += 1
            # preodr[i + 1] = k

        memmove(&preodr[pre_i + 1], &preodr[pre_i], <size_t>(
            (pre_i_after - pre_i) * sizeof(Py_ssize_t)
        ))

        preodr[pre_i] = link

        # postorder shift: all nodes after clade +2, tip and link inserted after clade
        for i in range(n - 1, post_i, -1):
            tree[postodr[i], 7] += 2
            # k = postodr[i]
            # tree[k, 7] += 2
            # postodr[i + 2] = k

        memmove(&postodr[post_i + 3], &postodr[post_i + 1], <size_t>(
            (n - post_i - 1) * sizeof(Py_ssize_t)
        ))

        postodr[post_i + 2] = link
        postodr[post_i + 1] = tip

        # size +1 from link to root
        curr = link
        while curr:
            parent = tree[curr, 2]
            tree[parent, 4] += 1
            curr = parent


def _insert_taxon_1(
    Py_ssize_t taxon,
    Py_ssize_t target,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] postodr,
    bint use_depth=True,
):
    r"""Insert a taxon between a target node and its parent.

    Move memory instead of +1 or +2.

    """
    cdef Py_ssize_t left, right, parent, sibling, size, depth, pre_i, post_i1
    cdef Py_ssize_t i, k, side, pre_i_after, curr

    # determine tree dimensions
    # typically n = 2 * taxon - 3, but this function doesn't enforce this
    cdef Py_ssize_t m = tree[0, 4]
    cdef Py_ssize_t n = m * 2 - 1
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

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
        node[4] = m + 1
        node[7] = n + 1

        # link
        node = &tree[link, 0]
        node[0] = left
        node[1] = right
        node[2] = 0
        node[3] = tip
        node[4] = m
        node[5] = 1
        node[6] = 1
        node[7] = n - 1

        # tip
        node = &tree[tip, 0]
        node[0] = 0
        node[1] = taxon
        node[2] = 0
        node[3] = link
        node[4] = 1
        node[5] = 1
        node[6] = n + 1
        node[7] = n

        # entire tree depth + 1
        if use_depth:
            for i in range(1, n):
                tree[i, 5] += 1
                tree[i, 6] += 1
        else:
            for i in range(1, n):
                tree[i, 6] += 1

        # preorder
        memmove(&preodr[2], &preodr[1], <size_t>((n - 1) * intsize))

        preodr[1] = link
        preodr[n + 1] = tip

        # postorder
        postodr[n - 1] = link
        postodr[n] = tip
        postodr[n + 1] = 0

    # Regular case (any other branch): The link becomes the parent of the target node,
    # and child of its original parent. Taxon k becomes the sibling
    else:
        node = &tree[target, 0]
        left = node[0]
        right = node[1]
        parent = node[2]
        sibling = node[3]
        size = node[4]
        depth = node[5]
        pre_i = node[6]
        post_i1 = node[7] + 1

        side = int(tree[parent, 0] != target)
        tree[parent, side] = link
        tree[sibling, 3] = link
        node[2] = link
        node[3] = tip

        # preorder index of node after clade
        pre_i_after = pre_i + size * 2 - 1

        # link
        node = &tree[link, 0]
        node[0] = target
        node[1] = tip
        node[2] = parent
        node[3] = sibling
        node[4] = size + 1
        node[5] = depth
        node[6] = pre_i
        node[7] = post_i1 + 1

        # tip
        node = &tree[tip, 0]
        node[0] = 0
        node[1] = taxon
        node[2] = link
        node[3] = target
        node[4] = 1
        node[5] = depth + 1
        node[6] = pre_i_after + 1
        node[7] = post_i1

        # preorder shift: nodes after clade +2, tip inserted after clade, nodes within
        # clade +1, link inserted before clade
        for i in range(pre_i_after, n):
            tree[preodr[i], 6] += 2

        memmove(&preodr[pre_i_after + 2], &preodr[pre_i_after], <size_t>(
            (n - pre_i_after) * intsize
        ))

        preodr[pre_i_after + 1] = tip

        # optional: clade depth +1
        if use_depth:
            for i in range(pre_i, pre_i_after):
                k = preodr[i]
                tree[k, 5] += 1
                tree[k, 6] += 1
        else:
            for i in range(pre_i, pre_i_after):
                tree[preodr[i], 6] += 1

        memmove(&preodr[pre_i + 1], &preodr[pre_i], <size_t>(
            (pre_i_after - pre_i) * intsize
        ))

        preodr[pre_i] = link

        # postorder shift: all nodes after clade +2, tip and link inserted after clade
        for i in range(post_i1, n):
            tree[postodr[i], 7] += 2

        memmove(&postodr[post_i1 + 2], &postodr[post_i1], <size_t>(
            (n - post_i1) * intsize
        ))

        postodr[post_i1 + 1] = link
        postodr[post_i1] = tip

        # size +1 from link to root
        curr = link
        while curr:
            parent = tree[curr, 2]
            tree[parent, 4] += 1
            curr = parent


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
    floating[::1] powers,
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
    cdef floating power, diff

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
            power = powers[tree[a, 5] - depth_2]
            jj = tree[a, 6]
            admx = &adm[a, 0]
            for j in range(jj + 1, jj + tree[a, 4] * 2 - 1):
                b = preodr[j]
                admx[b] = adm[b, a] = admx[b] - power * (
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
            adm[curr, a] = adm[a, curr] = adm[curr, a] - powers[i + 2] * diff

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
                admx[b] = adm[b, a] = admx[b] - powers[j + 2] * diff

            # current node vs. each descendant
            jj = tree[a, 6]
            power = powers[depth_diff + tree[a, 5]]  #### trick
            for j in range(jj + 1, jj + tree[a, 4] * 2 - 1):
                b = preodr[j]
                admx[b] = adm[b, a] = admx[b] - power * (
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
