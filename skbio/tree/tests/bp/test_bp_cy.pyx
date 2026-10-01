# cython: language_level=3
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy.testing as npt
import numpy as np
cimport numpy as cnp
from numpy cimport uint8_t as BOOL_t

cdef extern from "Python.h":
    cdef Py_ssize_t PY_SSIZE_T_MAX

from libc.math cimport ceil, log as ln, pow, log2

from skbio.tree import BPTree
from skbio.tree.bp._bp_cy cimport _BPKernel
from skbio.tree.bp._bp_binary_tree cimport (
    bt_node_from_left, bt_left_child, bt_right_child)

fig1_B = np.array([1, 1, 1, 0, 1, 0, 1, 1 ,0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0,
                   0, 0, 0], dtype=np.uint8)


def get_test_obj():
    return BPTree(fig1_B)._kernel


def test_rank():
    cdef _BPKernel obj = get_test_obj()
    counts_1 = fig1_B.cumsum()
    counts_0 = (1 - fig1_B).cumsum()
    for exp, t in zip((counts_1, counts_0), (1, 0)):
        for idx, e in enumerate(exp):
            npt.assert_equal(obj.rank(t, idx), e)


def test_select():
    cdef _BPKernel obj = get_test_obj()
    pos_1 = np.unique(fig1_B.cumsum(), return_index=True)[1] #- 1
    pos_0 = np.unique((1 - fig1_B).cumsum(), return_index=True)[1]

    for exp, t in zip((pos_1, pos_0), (1, 0)):
        for k in range(1, len(exp)):
            npt.assert_equal(obj.select(t, k), exp[k])


def test_rank_property():
    cdef _BPKernel obj = get_test_obj()
    for i in range(len(fig1_B)):
        npt.assert_equal(obj.rank(1, i) + obj.rank(0, i), i+1)


def test_rank_select_property():
    cdef _BPKernel obj = get_test_obj()
    pos_1 = np.unique(fig1_B.cumsum(), return_index=True)[1] #- 1
    pos_0 = np.unique((1 - fig1_B).cumsum(), return_index=True)[1]
    for t, pos in zip((0, 1), (pos_0, pos_1)):
        for k in range(len(pos)):
            # needed +t on expectation, unclear at this time why.
            npt.assert_equal(obj.rank(t, obj.select(t, k)), k + t)


def test_excess():
    cdef _BPKernel obj = get_test_obj()
    # from fig 2
    exp = [1, 2, 3, 2, 3, 2, 3, 4, 3, 2, 1, 2, 1, 2, 3, 4, 3, 4, 3, 2, 1, 0]
    for idx, e in enumerate(exp):
        npt.assert_equal(obj.excess(idx), e)


def test_depth():
    cdef _BPKernel obj = get_test_obj()
    # from fig 2
    exp = [1, 2, 3, 2, 3, 2, 3, 4, 3, 2, 1, 2, 1, 2, 3, 4, 3, 4, 3, 2, 1, 0]
    for idx, e in enumerate(exp):
        npt.assert_equal(obj.depth(idx), e)


def test_close():
    cdef _BPKernel obj = get_test_obj()
    exp = [21, 10, 3, 5, 9, 8, 12, 20, 19, 16, 18]
    for i, e in zip(np.argwhere(fig1_B == 1).squeeze(), exp):
        npt.assert_equal(obj.close(i), e)
        npt.assert_equal(obj.excess(obj.close(i)), obj.excess(i) - 1)


def test_open():
    cdef _BPKernel obj = get_test_obj()
    exp = [2, 4, 7, 6, 1, 11, 15, 17, 14, 13, 0]
    for i, e in zip(np.argwhere(fig1_B == 0).squeeze(), exp):
        npt.assert_equal(obj.open(i), e)
        npt.assert_equal(obj.excess(obj.open(i)) - 1,
                         obj.excess(i))


def test_enclose():
    cdef _BPKernel obj = get_test_obj()
    # i > 0 and i < (len(B) - 1)
    exp = [0, 1, 1, 1, 1, 1, 6, 6, 1, 0, 0, 0, 0, 13, 14, 14, 14, 14, 13, 0]
    for i, e in zip(range(1, len(fig1_B) - 1), exp):
        npt.assert_equal(obj.enclose(i), e)


def test_parent():
    cdef _BPKernel obj = get_test_obj()
    exp = [-1, 0, 1, 1, 1, 1, 1, 6, 6, 1, 0, 0, 0, 0, 13, 14, 14, 14, 14, 13,
           0, -1]
    for i, e in zip(range(len(fig1_B)), exp):
        npt.assert_equal(obj.parent(i), e)


def test_root():
    cdef _BPKernel obj = get_test_obj()
    npt.assert_equal(obj.root(), 0)


def test_is_tip():
    cdef _BPKernel obj = get_test_obj()

    exp = [0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0]
    for i, e in enumerate(exp):
        npt.assert_equal(obj.is_tip(i), e)


def test_first_child():
    cdef _BPKernel obj = get_test_obj()
    exp = [1, 2, 0, 0, 0, 0, 7, 0, 0, 7, 2, 0, 0, 14, 15, 0, 0, 0, 0, 15, 14,
           1]
    for i, e in enumerate(exp):
        npt.assert_equal(obj.first_child(i), e)


def test_last_child():
    cdef _BPKernel obj = get_test_obj()
    exp = [obj.preorder_select(7),
           obj.preorder_select(4),
           0,
           0,
           0,
           0,
           obj.preorder_select(5),
           0,
           0,
           obj.preorder_select(5),
           obj.preorder_select(4),
           0,
           0,
           obj.preorder_select(8),
           obj.preorder_select(10),
           0,
           0,
           0,
           0,
           obj.preorder_select(10),
           obj.preorder_select(8),
           obj.preorder_select(7)]
    for i, e in enumerate(exp):
        npt.assert_equal(obj.last_child(i), e)


def test_next_sibling():
    cdef _BPKernel obj = get_test_obj()
    exp = [0, 11, 4, 4, 6, 6, 0, 0, 0, 0, 11, 13, 13, 0, 0, 17, 17, 0, 0, 0, 0,
           0]
    for i, e in enumerate(exp):
        npt.assert_equal(obj.next_sibling(i), e)


def test_previous_sibling():
    cdef _BPKernel obj = get_test_obj()
    exp = [0, 0, 0, 0, 2, 2, 4, 0, 0, 4, 0, 1, 1, 11, 0, 0, 0, 15, 15, 0, 11,
           0]
    for i, e in enumerate(exp):
        npt.assert_equal(obj.previous_sibling(i), e)


def test_fwdsearch():
    cdef _BPKernel obj = get_test_obj()
    exp = {(0, 0): 10,   # close of first child
           (3, -2): 21,  # close of root
           (11, 2): 15}  # from one tip to the next

    for (i, d), e in exp.items():
        npt.assert_equal(obj.fwdsearch(i, d), e)


def test_bwdsearch():
    cdef _BPKernel obj = get_test_obj()
    exp = {(3, 0): 1,  # open of parent
           (21, 4): 17,  # nested tip
           (9, 2): 7}  # open of the node

    for (i, d), e in exp.items():
        npt.assert_equal(obj.bwdsearch(i, d), e)


def test_fwdsearch_more():
    cdef _BPKernel bp
    from skbio.tree.bp import parse_newick
    bp = parse_newick('((a,b,(c)),d,((e,f)));')._kernel

    # simulating close so only testing open parentheses. A "close" on a closed
    # parenthesis does not make sense, so the result is not useful.
    # In practice, an "close" method should ensure it is operating on a closed
    # parenthesis.
    # [(open_idx, close_idx), ...]
    exp = [(1, 10), (0, 21), (2, 3), (4, 5), (6, 9), (7, 8), (11, 12),
           (13, 20), (14, 19), (15, 16), (17, 18)]

    for open_, exp_close in exp:
        obs_close = bp.fwdsearch(open_, -1)
        assert obs_close == exp_close

    # slightly modified version of fig2 with an extra child forcing a test
    # of the direct sibling check with negative partial excess

    # this translates into:
    # 012345678901234567890123
    # ((()()(()))()((()()())))
    bp = parse_newick('((a,b,(c)),d,((e,f,g)));')._kernel

    # [(open_idx, close_idx), ...]
    exp = [(0, 23), (1, 10), (2, 3), (4, 5), (6, 9), (7, 8), (11, 12),
           (13, 22), (14, 21), (15, 16), (17, 18), (19, 20)]

    for open_, exp_close in exp:
        obs_close = bp.fwdsearch(open_, -1)
        assert obs_close == exp_close


def test_bwdsearch_more():
    cdef _BPKernel bp
    from skbio.tree.bp import parse_newick
    bp = parse_newick('((a,b,(c)),d,((e,f)));')._kernel

    # simulating open so only testing closed parentheses.
    # [(close_idx, open_idx), ...]
    exp = [(21, 0), (8, 7), (9, 6), (10, 1), (3, 2), (5, 4), (12, 11),
           (16, 15), (20, 13), (19, 14), (18, 17)]

    for close_, exp_open in exp:
        obs_open = bp.bwdsearch(close_, 0) + 1
        assert obs_open == exp_open

    # slightly modified version of fig2 with an extra child forcing a test
    # of the direct sibling check with negative partial excess

    # this translates into:
    # 012345678901234567890123
    # ((()()(()))()((()()())))
    bp = parse_newick('((a,b,(c)),d,((e,f,g)));')._kernel

    # [(close_idx, open_idx), ...]
    exp = [(23, 0), (10, 1), (3, 2), (5, 4), (9, 6), (8, 7), (12, 11),
           (22, 13), (21, 14), (16, 15), (18, 17), (20, 19)]

    for close_, exp_open in exp:
        obs_open = bp.bwdsearch(close_, 0) + 1
        assert obs_open == exp_open


def test_scan_block_forward():
    cdef _BPKernel bp
    from skbio.tree.bp import parse_newick
    bp = parse_newick('((a,b,(c)),d,((e,f)));')._kernel

    # [(open, close), ...]
    b = 4
    d = -1
    exp_b_4 = [(0, ((0, -1), (1, -1), (2, 3), (3, -1))),
               (1, ((4, 5), (5, -1), (6, -1), (7, -1))),
                   # 8 and 9 are nonsensical from finding a "close" perspective
               (2, ((8, 9), (9, 10), (10, -1), (11, -1))),
               (3, ((12, -1), (13, -1), (14, -1), (15, -1))),
                   # 16 and 18 are nonsensical from a "close" perspective
               (4, ((16, 19), (17, 18), (18, 19), (19, -1))),
                   # 20 is nonsensical from finding a "close" perspective
               (5, ((20, 21), (21, -1)))]

    for k, exp_results in exp_b_4:
        for idx, exp_result in exp_results:
            obs_result = bp.scan_block_forward(idx, k, b, bp.excess(idx) + d)
            assert obs_result == exp_result

    b = 8
    exp_b_8 = [(0, ((0, -1), (1, -1), (2, 3), (3, -1),
                    (4, 5), (5, -1), (6, -1), (7, -1))),
               (1, ((8, 9), (9, 10), (10, -1), (11, 12),
                    (12, -1), (13, -1), (14, -1), (15, -1))),
               (2, ((16, 19), (17, 18), (18, 19), (19, 20),
                    (20, 21), (21, -1)))]

    for k, exp_results in exp_b_8:
        for idx, exp_result in exp_results:
            obs_result = bp.scan_block_forward(idx, k, b, bp.excess(idx) + d)
            assert obs_result == exp_result


def test_scan_block_backward():
    cdef _BPKernel bp
    from skbio.tree.bp import parse_newick
    bp = parse_newick('((a,b,(c)),d,((e,f)));')._kernel

    # adding +1 to simluate "open" so calls on open parentheses are weird
    # [(open, close), ...]
    b = 4
    d = 0
    exp_b_4 = [(0, ((0, 0), (1, 0), (2, 0), (3, 2))),
               (1, ((4, 0), (5, 4), (6, 5), (7, 0))),
               (2, ((8, 0), (9, 0), (10, 0), (11, 10))),
               (3, ((12, 0), (13, 12), (14, 0), (15, 0))),
               (4, ((16, 0), (17, 16), (18, 17), (19, 0))),
               (5, ((20, 0), (21, 0)))]

    for k, exp_results in exp_b_4:
        for idx, exp_result in exp_results:
            obs_result = bp.scan_block_backward(idx, k, b, bp.excess(idx) + d)
            obs_result += 1  # simulating open
            assert obs_result == exp_result

    b = 8
    exp_b_8 = [(0, ((0, 0), (1, 0), (2, 0), (3, 2),
                    (4, 3), (5, 4), (6, 5), (7, 0))),
               (1, ((8, 0), (9, 0), (10, 0), (11, 10),
                    (12, 11), (13, 12), (14, 9), (15, 8))),
               (2, ((16, 0), (17, 16), (18, 17), (19, 0),
                    (20, 0), (21, 0)))]

    for k, exp_results in exp_b_8:
        for idx, exp_result in exp_results:
            obs_result = bp.scan_block_backward(idx, k, b, bp.excess(idx) + d)
            obs_result += 1  # simulating open
            assert obs_result == exp_result


def test_rmm():
    from skbio.tree.bp import parse_newick
    # test tree is ((a,b,(c)),d,((e,f)));
    # this is from fig 2 of Cordova and Navarro:
    # http://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf
    tree = parse_newick('((a,b,(c)),d,((e,f)));')
    exp = np.array([[0, 1, 0, 1, 1, 0, 0, 1, 2, 1, 1, 2, 0],   # m
                    [4, 4, 4, 4, 4, 4, 0, 3, 4, 3, 4, 4, 1]],  # M
                   dtype=np.intp)
    npt.assert_equal(tree._m, exp[0])
    npt.assert_equal(tree._M, exp[1])

    # and the scan-based reference construction agrees
    ref = reference_index(tree.data)
    npt.assert_equal(ref['m'], exp[0])
    npt.assert_equal(ref['M'], exp[1])


def reference_index(cnp.ndarray[BOOL_t, ndim=1] B):
    """Scan-based construction of the BP navigation index.

    A direct port of the original compiled construction (the rmM tree of
    Navarro and Sadakane, http://www.dcc.uchile.cl/~gnavarro/ps/talg12.pdf, and
    the excess and select indexes), kept as the oracle for the vectorized
    ``skbio.tree.bp._bp._build_index``.
    """
    cdef:
        Py_ssize_t B_size = B.shape[0]
        Py_ssize_t b, n_tip, height, n_internal, n_total
        Py_ssize_t i, j, k, lvl, pos, node, lchild, rchild, offset
        Py_ssize_t lower_limit, upper_limit
        Py_ssize_t min_, max_, excess = 0, r = 0, rank
        Py_ssize_t[:, ::1] mM
        Py_ssize_t[::1] rr, e_index

    b = <Py_ssize_t>ceil(ln(<double> B_size) * ln(ln(<double> B_size)))
    if b < 1:
        b = 1
    n_tip = <Py_ssize_t>ceil(B_size / <double> b)
    height = <Py_ssize_t>ceil(log2(n_tip))
    n_internal = <Py_ssize_t>(pow(2, height)) - 1
    n_total = n_tip + n_internal

    mM = np.zeros((n_total, 2), dtype=np.intp)
    rr = np.zeros(n_total, dtype=np.intp)

    i = 0
    while i < B_size:
        offset = i // b
        lower_limit = i
        upper_limit = min(i + b, B_size)
        min_ = PY_SSIZE_T_MAX
        max_ = 0

        rr[offset + n_internal] = r
        for j in range(lower_limit, upper_limit):
            excess += -1 + (2 * B[j])
            r += B[j]
            if excess < min_:
                min_ = excess
            if excess > max_:
                max_ = excess

        mM[offset + n_internal, 0] = min_
        mM[offset + n_internal, 1] = max_
        i += b

    for lvl in range(height - 1, -1, -1):
        for pos in range(<Py_ssize_t>pow(2, lvl)):
            node = bt_node_from_left(pos, lvl)
            lchild = bt_left_child(node)
            rchild = bt_right_child(node)

            if lchild >= n_total:
                continue
            elif rchild >= n_total:
                mM[node, 0] = mM[lchild, 0]
                mM[node, 1] = mM[lchild, 1]
            else:
                mM[node, 0] = min(mM[lchild, 0], mM[rchild, 0])
                mM[node, 1] = max(mM[lchild, 1], mM[rchild, 1])

            rr[node] = rr[lchild]

    # excess via rank, as the original _excess: 2 * rank(1, i) - i - 1, where
    # rank(1, i) is the block's starting rank plus a scan within the block
    e_index = np.empty(B_size, dtype=np.intp)
    for i in range(B_size):
        k = i // b
        rank = rr[bt_node_from_left(k, height)]
        for j in range(k * b, min((k + 1) * b, B_size, i + 1)):
            rank += B[j]
        e_index[i] = 2 * rank - i - 1

    step = B.astype(bool)
    step[0] = True
    k_index_1 = np.flatnonzero(step).astype(np.intp)
    step = (B == 0)
    step[0] = True
    k_index_0 = np.flatnonzero(step).astype(np.intp)

    return {'e_index': np.asarray(e_index), 'k_index_0': k_index_0,
            'k_index_1': k_index_1, 'm': np.asarray(mM[:, 0]),
            'M': np.asarray(mM[:, 1]), 'r': np.asarray(rr), 'b': b,
            'height': height}


def kernel_index_op(tree, str op, Py_ssize_t a, Py_ssize_t b=0):
    """Call an index operation of the Cython engine that Python cannot reach.

    The oracle of the Numba engine's parity tests (``test_bp_numba``) for the
    ``cdef`` methods of ``_BPKernel``: ``rank(t, i)``, ``select(t, k)``,
    ``excess(i)``, ``fwdsearch(i, d)``, ``bwdsearch(i, d)``, ``open(i)`` and
    ``enclose(i)``.
    """
    cdef _BPKernel k = tree._kernel
    if op == "rank":
        return k.rank(a, b)
    elif op == "select":
        return k.select(a, b)
    elif op == "excess":
        return k.excess(a)
    elif op == "fwdsearch":
        return k.fwdsearch(a, b)
    elif op == "bwdsearch":
        return k.bwdsearch(a, b)
    elif op == "open":
        return k.open(a)
    elif op == "enclose":
        return k.enclose(a)
    raise ValueError(op)
