# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

cimport cython
import numpy as np
cimport numpy as cnp
cnp.import_array()


def _preorder(
    Py_ssize_t[::1] order,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] stack,
    Py_ssize_t start = 0,
):
    r"""Perform preorder traversal.

    This function and :func:`_postorder` use stacks to avoid recursion. The stack
    array is pre-allocated. The output (ordered nodes) is also written into a
    pre-allocated array.

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
    double[:, ::1] adm,
    double[:, :] dm,
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
       subtree is the node's parent. Its immediate children are the parent's parent
       and the node's sibling.

    Then it iteratively applies Eq. 2 to calculate the average distance between two
    subtrees based on the average distances from the children of one of the subtrees
    to the other.

        d(A, B) = (|B_1| * d(A, B_1) + |B_2| * d(A, B_2)) / |B|

    """
    cdef Py_ssize_t i, j, k
    cdef Py_ssize_t curr, parent, sibling, p_size, s_size
    cdef Py_ssize_t a, a_size, a_taxon
    cdef Py_ssize_t a1 = 0, a2 = 0, a1_size = 0, a2_size = 0
    cdef Py_ssize_t b, b_size, b1, b2

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
            # size of the slice is taxon count * 2 - 2. *
            k = tree[sibling, 7]
            for j in range(k - tree[sibling, 4] * 2 + 2, k + 1):
                b = postodr[j]
                b_size = tree[b, 4]

                # If both a and b are tips, take the original taxon-to-taxon distance
                # (A4.1 (a) i).
                if b_size == 1 and a_taxon != 0:
                    adm[a, b] = adm[b, a] = dm[a_taxon, tree[b, 1]]

                # If a is an internal node, and b is either (a tip or an internal node),
                # calculate the average distance based on the two child subtrees of a
                # (A4.1 (a) ii).
                elif a_taxon == 0:
                    adm[a, b] = adm[b, a] = (
                        a1_size * adm[a1, b] + a2_size * adm[a2, b]
                    ) / a_size

                # If a is a tip, and b is an internal node, calculate the average
                # distance based on the two child subtrees of b (A4.1 (a) iii).
                else:
                    b1, b2 = tree[b, 0], tree[b, 1]
                    adm[a, b] = adm[b, a] = (
                        tree[b1, 4] * adm[a, b1] + tree[b2, 4] * adm[a, b2]
                    ) / b_size

            curr = parent

    # Step 2: Calculate subtree to root (taxon 0) distances (A4.1 (b)).

    # This is done through a postorder traversal.
    for i in range(n - 1):
        a = postodr[i]
        a_size = tree[a, 4]
        if a_size == 1:
            a_taxon = tree[a, 1]
            adm[a, 0] = adm[0, a] = dm[0, a_taxon]
        else:
            a1, a2 = tree[a, 0], tree[a, 1]
            adm[a, 0] = adm[0, a] = (
                tree[a1, 4] * adm[a1, 0] + tree[a2, 4] * adm[a2, 0]
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
        k = tree[a, 7]
        for j in range(k - tree[a, 4] * 2 + 2, k):
            b = postodr[j]
            adm[a, b] = adm[b, a] = (
                s_size * adm[b, sibling] + p_size * adm[b, parent]
            ) / a_size


def _bal_avgdist_matrix(
    double[:, ::1] adm,
    double[:, :] dm,
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

    cdef Py_ssize_t n = 2 * tree[0, 4] - 1

    # Step 1: Calculate non-nested subtree to subtree distances.
    for i in range(n - 1):
        a = postodr[i]
        if tree[a, 0] == 0:
            a_taxon = tree[a, 1]
        else:
            a_taxon = 0
            a1, a2 = tree[a, 0], tree[a, 1]
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
                    adm[a, b] = adm[b, a] = dm[a_taxon, tree[b, 1]]
                elif a_taxon == 0:
                    adm[a, b] = adm[b, a] = 0.5 * (adm[a1, b] + adm[a2, b])
                else:
                    adm[a, b] = adm[b, a] = 0.5 * (
                        adm[a, tree[b, 0]] + adm[a, tree[b, 1]]
                    )
            curr = parent

    # Step 2: Calculate subtree to root distances.
    for i in range(n - 1):
        a = postodr[i]
        if tree[a, 0] == 0:
            adm[a, 0] = adm[0, a] = dm[0, tree[a, 1]]
        else:
            adm[a, 0] = adm[0, a] = 0.5 * (adm[tree[a, 0], 0] + adm[tree[a, 1], 0])

    # Step 3: Calculate nested subtree to subtree distances.
    for i in range(1, n):
        a = preodr[i]
        parent, sibling, k = tree[a, 2], tree[a, 3], tree[a, 7]
        for j in range(k - tree[a, 4] * 2 + 2, k):
            b = postodr[j]
            adm[a, b] = adm[b, a] = 0.5 * (adm[b, sibling] + adm[b, parent])


def _avgdist_taxon(
    double[:, ::1] adk,
    Py_ssize_t taxon,
    double[:, :] dm,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] postodr,
):
    """Calculate average distances between a new taxon and existing subtrees.

    This function will update adk, a float array of (n, 2) in which columns 0 and 1
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

    # Calculate the distance between taxon and the lower subtree of each node.
    for i in range(n):
        node = postodr[i]
        left, right = tree[node, 0], tree[node, 1]
        if left == 0:
            adk[node, 0] = dm[taxon, right]
        else:
            adk[node, 0] = (
                tree[left, 4] * adk[left, 0] + tree[right, 4] * adk[right, 0]
            ) / tree[node, 4]

    # Assign upper distance of root.
    adk[0, 1] = dm[taxon, 0]

    # Calculate the distance between taxon k and the upper subtree of each node.
    for i in range(1, n):
        node = preodr[i]
        parent, sibling = tree[node, 2], tree[node, 3]
        adk[node, 1] = (
            (m - tree[parent, 4]) * adk[parent, 1] + tree[sibling, 4] * adk[sibling, 0]
        ) / (m - tree[node, 4])


def _bal_avgdist_taxon(
    double[:, ::1] adk,
    Py_ssize_t taxon,
    double[:, :] dm,
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
    for i in range(n):
        node = postodr[i]
        left, right = tree[node, 0], tree[node, 1]
        if left == 0:
            adk[node, 0] = dm[taxon, right]
        else:
            adk[node, 0] = 0.5 * (adk[left, 0] + adk[right, 0])
    adk[0, 1] = dm[taxon, 0]
    for i in range(1, n):
        node = preodr[i]
        adk[node, 1] = 0.5 * (adk[tree[node, 2], 1] + adk[tree[node, 3], 0])


def _ols_lengths(
    double[::1] lens,
    double[:, ::1] adm,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
):
    r"""Calculate branch lengths of a tree based on the OLS framework.

    Using an average distance matrix between all pairs of subtrees.

    Implemented according to Eqs. 3 & 4 of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, left, right, parent, sibling
    cdef Py_ssize_t l_size, r_size, p_size, s_size
    cdef double lambda_
    
    # total numbers of taxa and nodes in the tree
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3

    for i in range(1, n):
        node = preodr[i]
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
            lambda_ = <double>(p_size * r_size + s_size * l_size) / (
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
    double[::1] lens,
    double[:, ::1] ad2,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
):
    r"""Calculate branch lengths of a tree based on an OLS framework.

    Using only average distances between pairs of distant-2 subtrees.

    This function produces the same result as `_ols_lengths`. The latter
    relies on Eq. 3 of Desper and Gascuel (2002), which involves distances between
    distant-3 subtrees. Therefore, we need to modify the equation such that it only
    takes distances between distant-2 subtrees as input.

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
    cdef Py_ssize_t i
    cdef Py_ssize_t node, left, sibling
    
    # total numbers of taxa and nodes in the tree
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3

    for i in range(1, n):
        node = preodr[i]
        left = tree[node, 0]

        # External (terminal) branch: based on Eq. 4 of the paper.
        if left == 0:
            lens[node] = 0.5 * (ad2[node, 0] + ad2[node, 1] - ad2[tree[node, 3], 0])

        # Internal branch: based on the equation discussed above.
        else:
            sibling = tree[node, 3]
            lens[node] = 0.5 * (
                ad2[left, 0]
                + ad2[tree[node, 1], 0]
                + ad2[node, 0]
                + ad2[node, 1]
                - ad2[sibling, 0]
                - ad2[left, 1]
            ) - (
                tree[sibling, 4]
                * ad2[node, 1]
                + (m - tree[tree[node, 2], 4])
                * ad2[node, 0]
            ) / (m - tree[node, 4])

    # root branch
    left = tree[0, 0]
    lens[0] = 0.5 * (ad2[left, 0] + ad2[tree[0, 1], 0] - ad2[left, 1])


def _bal_lengths(
    double[::1] lens,
    double[:, ::1] adm,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
):
    r"""Calculate branch lengths of a tree based on the balanced framework.

    Using a balanced average distance matrix between all pairs of subtrees.

    This function resembles :func:`_ols_lengths` but it uses the balanced
    framework.

    Implemented according to Eqs. 3 & 4 and the description on top of pg. 691 of Desper
    and Gascuel (2002).

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, left, right, parent, sibling

    for i in range(1, 2 * tree[0, 4] - 1):
        node = preodr[i]
        left = tree[node, 0]
        parent = tree[node, 2]
        sibling = tree[node, 3]
        if left == 0:
            lens[node] = 0.5 * (
                adm[parent, node] + adm[sibling, node] - adm[parent, sibling]
            )
        else:
            right = tree[node, 1]
            lens[node] = 0.25 * (
                adm[parent, left]
                + adm[sibling, right]
                + adm[parent, right]
                + adm[sibling, left]
            ) - 0.5 * (adm[parent, sibling] + adm[left, right])
    left, right = tree[0, 0], tree[0, 1]
    lens[0] = 0.5 * (adm[left, 0] + adm[right, 0] - adm[left, right])


def _ols_min_branch_d2(
    double[::1] lens,
    double[:, ::1] ad2,
    double[:, ::1] adk,
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
    cdef double L, numerator, lambda_0, lambda_1

    cdef Py_ssize_t min_node = 0
    cdef double min_len = 0
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3

    lens[min_node] = min_len
    for i in range(1, n):
        node = preodr[i]
        parent = tree[node, 2]
        sibling = tree[node, 3]
        size = tree[node, 4]
        p_size = m - tree[parent, 4]
        s_size = tree[sibling, 4]

        numerator = s_size + size * p_size
        lambda_0 = numerator / (s_size + size) / (p_size + 1)
        lambda_1 = numerator / (s_size + p_size) / (size + 1)

        lens[node] = L = lens[parent] + 0.5 * (
            (lambda_0 - lambda_1) * (adk[sibling, 0] + ad2[node, 0])
            + (lambda_1 - 1) * (ad2[sibling, 1] + adk[parent, 1])
            + (1 - lambda_0) * (ad2[sibling, 0] + adk[node, 0])
        )

        if L < min_len:
            min_len, min_node = L, node

    return min_node


def _bal_min_branch(
    double[::1] lens,
    double[:, ::1] adm,
    double[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
):
    """Find the branch with the minimum length change after inserting a new taxon.

    This function resembles :func:`_ols_min_branch_d2` but it 1) uses the
    balanced framework and 2) calculates based on the entire matrix. See also
    the important note of the latter.

    Implemented according to Eq. 10 of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t i
    cdef Py_ssize_t node, parent, sibling
    cdef double L

    cdef Py_ssize_t min_node = 0
    cdef double min_len = 0

    lens[min_node] = min_len
    for i in range(1, 2 * tree[0, 4] - 1):
        node = preodr[i]
        parent = tree[node, 2]
        sibling = tree[node, 3]
        lens[node] = L = lens[parent] + 0.25 * (
            adm[sibling, parent] + adk[node, 0] - adm[sibling, node] - adk[parent, 1]
        )
        if L < min_len:
            min_len, min_node = L, node
    return min_node


def _avgdist_d2_insert(
    double[:, ::1] ad2,
    Py_ssize_t target,
    double[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
):
    r"""Update average distances between distant-2 subtrees after taxon insertion.

    This function will update ad2, a float array of (n, 2) representing pairwise
    distances between all distant-2 subtrees in the tree. Here, `distant-2 subtrees`
    refer to subtrees that are two branches away from each other. Specifically, there
    are two scenarios:

    - Column 0: Distance between the lower subtree of the current node and the upper
      subtree of its parent.
    - Column 1: Distance between the lower subtree of the current node and the lower
      subtree of its sibling.

    This function assumes that the taxon will be inserted into the branch connecting
    the target node and its parent. After insertion, the taxon will become the sibling
    of the target.

               parent
                /  \
             link  sibling
             /  \
        target  taxon (k)

    This function should be executed *before* calling :func:`_insert_taxon`, which will
    mutate the tree.

    Implemented according to Eq. 8 of Desper and Gascuel (2002).

    """
    cdef Py_ssize_t i, ii
    cdef Py_ssize_t left, right, parent, sibling, size, curr, p_size, size_1

    # dimensions and positions
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    ###### Special case: insert into the root branch. ######

    if target == 0:
        # k (lower) to root (parent, upper): pre-calculated
        ad2[tip, 0] = adk[0, 1]

        # k (lower) to link (sibling, lower): equals to k to root (lower).
        ad2[link, 1] = ad2[tip, 1] = adk[0, 0]

        # Link (lower) to root (parent, upper): de novo calculation according to the
        # equation in A4.1(b). It is basically the distance between the upper and lower
        # subtrees of the root itself.
        left, right = tree[0, 0], tree[0, 1]
        ad2[link, 0] = (
            tree[left, 4] * ad2[left, 0] + tree[right, 4] * ad2[right, 0]
        ) / tree[0, 4]

        # Calculate all node (lower) to parent (upper, containing k) distances. These
        # parents include the new link.
        for node in range(1, n):
            p_size = m - tree[tree[node, 2], 4]
            ad2[node, 0] = (adk[node, 0] + ad2[node, 0] * p_size) / (p_size + 1)

        return

    ###### Regular case: insert into any other branch. ######

    parent, sibling, size = tree[target, 2], tree[target, 3], tree[target, 4]
    size_1 = size + 1

    # Temporarily copy the distances of target to link (will edit later).
    ad2[link, 0] = ad2[target, 0]
    ad2[link, 1] = ad2[target, 1]

    # Distance between k (lower) and link (parent, upper) equals to that between k and
    # the upper subtree of target.
    ad2[tip, 0] = adk[target, 1]

    # Distance between target (lower) and link (parent, upper) needs to be calculated
    # using the equation in A4.1(c). Basically, it is the distance between the lower
    # and upper subtrees of the same target.
    ad2[target, 0] = (
        tree[sibling, 4] * ad2[target, 1] + (m - tree[parent, 4]) * ad2[target, 0]
    ) / (m - size)

    # Transfer the pre-calculated distance between target (lower) and k (sibling,
    # lower).
    ad2[target, 1] = ad2[tip, 1] = adk[target, 0]

    # Within the clade below target, calculate the distance between each node (lower)
    # and its parent (upper, containing k).
    ii = tree[target, 6]
    for i in range(ii + 1, ii + size * 2 - 1):
        node = preodr[i]
        p_size = m - tree[tree[node, 2], 4]
        ad2[node, 0] = (adk[node, 0] + ad2[node, 0] * p_size) / (p_size + 1)

    # Iterate over the ancestors of target, starting from link and ending at root.
    curr = link
    while curr:
        # Calculate the distance between each pair of lower (containing k) and upper
        # ancestors.
        ad2[curr, 0] = (adk[parent, 1] + ad2[curr, 0] * size) / size_1

        # Calculate the distance between each ancestor (lower, containing k) and its
        # sibling (lower).
        ad2[curr, 1] = ad2[sibling, 1] = (
            adk[sibling, 0] + ad2[curr, 1] * size
        ) / size_1

        # Within the clade below each sibling, calculate the distance between each node
        # (lower) and its parent (upper, containing k).
        ii = tree[sibling, 6]
        for i in range(ii + 1, ii + tree[sibling, 4] * 2 - 1):
            node = preodr[i]
            p_size = m - tree[tree[node, 2], 4]
            ad2[node, 0] = (adk[node, 0] + ad2[node, 0] * p_size) / (p_size + 1)

        curr = parent
        parent, sibling, size = tree[curr, 2], tree[curr, 3], tree[curr, 4]
        size_1 = size + 1


def _bal_avgdist_insert(
    double[:, ::1] adm,
    Py_ssize_t target,
    double[:, ::1] adk,
    Py_ssize_t[:, ::1] tree,
    Py_ssize_t[::1] preodr,
    Py_ssize_t[::1] postodr,
    double[::1] powers,
    Py_ssize_t[::1] stack,
):
    r"""Update balanced average distance matrix after taxon insertion.

    This function resembles :func:`_avgdist_d2_insert` but it 1) uses the balanced
    framework and 2) updates the entire matrix.

    Two additional parameters are provided: `powers` is a pre-calculated array of
    2^(-l) powers (l is the depth difference between two nodes). `stack` is an
    integer array to store ancestral nodes of target.

    """
    cdef Py_ssize_t i, j, ii, jj, anc_i
    cdef Py_ssize_t parent, sibling, depth
    cdef Py_ssize_t curr, anc, cousin, depth_1
    cdef Py_ssize_t a, b
    cdef double power, diff

    # dimensions and positions
    cdef Py_ssize_t m = tree[0, 4] + 1
    cdef Py_ssize_t n = 2 * m - 3
    cdef Py_ssize_t link = n
    cdef Py_ssize_t tip = n + 1

    ###### Special case: insert into the root branch. ######

    if target == 0:
        # Transfer distance between k and root (upper).
        adm[0, tip] = adm[tip, 0] = adk[0, 1]

        # k to link: equals to k to root (lower).
        adm[link, tip] = adm[tip, link] = adk[0, 0]

        # Root to link: de novo calculation according to the equation in A4.1(b). It is
        # basically the distance between the upper and lower subtrees of the root.
        a1, a2 = tree[0, 0], tree[0, 1]
        adm[0, link] = adm[link, 0] = 0.5 * (adm[a1, 0] + adm[a2, 0])

        # Iterate over all nodes but the root.
        for a in range(1, n):

            # Transfer distances between the node (lower) and k.
            adm[a, tip] = adm[tip, a] = adk[a, 0]
            
            # Calculate the distance between the node (lower) and link (upper, with two
            # taxa (0 and k) added) using Eq. 8.
            adm[a, link] = adm[link, a] = 0.5 * (adk[a, 0] + adm[a, 0])

            # Calculate the distances between the node (upper, containing k) and each
            # of its descendant (lower).
            ii = tree[a, 7]
            power = powers[tree[a, 5] + 1]
            for i in range(ii - tree[a, 4] * 2 + 2, ii):
                b = postodr[i]
                adm[a, b] = adm[b, a] = adm[a, b] + power * (adk[b, 0] - adm[0, b])

        return

    ###### Regular case: insert into any other branch. ######

    parent, sibling, depth = tree[target, 2], tree[target, 3], tree[target, 5]
    depth_1 = depth + 1

    ### Step 1: Distances within the clade below target. ###

    # Locate the clade below target (including target)
    ii = tree[target, 7]
    for i in range(ii - tree[target, 4] * 2 + 2, ii + 1):
        a = postodr[i]

        # Transfer pre-calculated distance between k (lower) and any node within the clade
        # (including target, lower).
        adm[a, tip] = adm[tip, a] = adk[a, 0]

        # Distance from any descendant (lower) to link (upper) equals to that to target.
        # (The last assignment: target to link is incorrect. It will be fixed below.)
        adm[a, link] = adm[link, a] = adm[a, target]

        # Within the clade, find all ancestor (a) - descendant (b) pairs, and calculate the
        # distance between the upper subtree of a (containing k) and the lower subtree of b.
        jj = tree[a, 7]
        power = powers[tree[a, 5] - depth + 1]
        for j in range(jj - tree[a, 4] * 2 + 2, jj):
            b = postodr[j]
            adm[a, b] = adm[b, a] = adm[a, b] + power * (adk[b, 0] - adm[target, b])

    ### Step 2: Distances around the insertion point. ###

    # Distance between k (lower) and link (upper) equals to that between k and the
    # upper subtree of target.
    adm[tip, link] = adm[link, tip] = adk[target, 1]

    # Distance between target (lower) and link (upper) needs to be calculated using the
    # equation in A4.1(c). Basically, it is the distance between the lower and upper
    # subtrees of the same target.
    adm[target, link] = adm[link, target] = 0.5 * (
        adm[target, sibling] + adm[target, parent]
    )

    ### Step 3: Distances among nodes outside the clade. ###

    # Iterate over ancestors of target in ascending order.
    anc_i = 0
    curr = target
    while curr:
        stack[anc_i] = anc = tree[curr, 2]

        # Transfer the pre-calculated distance between k and the ancestor (upper).
        adm[anc, tip] = adm[tip, anc] = adk[anc, 1]

        # Calculate the distance between link (lower, containing k) and the ancestor
        # (upper).
        adm[anc, link] = adm[link, anc] = 0.5 * (adk[anc, 1] + adm[anc, target])

        # Calculate the distance between each previous ancestor (lower, containing k)
        # and the current ancestor (upper).
        diff = adk[anc, 1] - adm[target, anc]
        for i in range(anc_i):
            a = stack[i]
            power = powers[depth_1 - tree[a, 5]]
            adm[anc, a] = adm[a, anc] = adm[anc, a] + power * diff

        # Identify the sibling clade descending from the ancestor.
        cousin = tree[curr, 3]
        ii = tree[cousin, 7]
        for i in range(ii - tree[cousin, 4] * 2 + 2, ii + 1):
            a = postodr[i]

            # Transfer the pre-calculated distances between k and each descendant
            # (lower).
            adm[a, tip] = adm[tip, a] = adk[a, 0]

            # Calculate the distance between link (lower, containing k) and each
            # descendant (lower).
            adm[a, link] = adm[link, a] = 0.5 * (adk[a, 0] + adm[a, target])

            # Calculate the distance between each previous ancestor (lower, containing
            # k) and each descendant (lower).
            diff = adk[a, 0] - adm[a, target]
            for j in range(anc_i):
                b = stack[j]
                adm[a, b] = adm[b, a] = adm[a, b] + powers[depth_1 - tree[b, 5]] * diff

            # Iterate over descendants of each member of the clade, and calculate the
            # distance between the former (upper, containing k) and the latter (lower).
            jj = tree[a, 7]
            power = powers[depth + tree[a, 5] - 2 * tree[anc, 5]]
            for j in range(jj - tree[a, 4] * 2 + 2, jj):
                b = postodr[j]
                adm[a, b] = adm[b, a] = adm[a, b] + power * (
                    adk[b, 0] - adm[b, target]
                )

        curr = anc
        anc_i += 1
