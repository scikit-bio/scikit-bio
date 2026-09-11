# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import lru_cache
from itertools import product
from io import StringIO
import unittest

import numpy as np
import numpy.testing as npt

from skbio import DNA, RNA, Protein, Sequence, TreeNode, SubstitutionMatrix
from skbio.alignment import AlignPath, multi_align, pair_align, align_score
from skbio.alignment._utils import encode_sequences
from skbio.alignment._pair import _encode_path
from skbio.alignment._multi import (
    _align_profiles,
    _ProfileWorkspace,
    _merge_profiles,
    _multi_distances,
    _fd_dist,
    _align_pair,
    _align_pair_roll,
    _guide_merges,
)


@lru_cache(None)
def _paths(p, q):
    """Enumerate complete paths without a DP recurrence or production helpers."""
    if p == q == 0:
        return ("",)
    result = []
    for move, di, dj in [("D", 1, 1), ("X", 1, 0), ("Y", 0, 1)]:
        if p >= di and q >= dj:
            result.extend(move + x for x in _paths(p - di, q - dj))
    return tuple(result)


def _score_path(path, scores, gap, free):
    """Score columns and maximal move runs independently of affine DP states."""
    o, e = (0, gap) if np.isscalar(gap) else gap
    p, q = scores.shape
    i = j = t = 0
    score = 0.0
    while t < len(path):
        move = path[t]
        if move == "D":
            score += float(scores[i, j])
            i += 1
            j += 1
            t += 1
            continue
        stop = t + 1
        while stop < len(path) and path[stop] == move:
            stop += 1
        length = stop - t
        terminal = j in (0, q) if move == "X" else i in (0, p)
        if not (free and terminal):
            score -= o + e * length
        if move == "X":
            i += length
        else:
            j += length
        t = stop
    assert (i, j) == (p, q)
    return score


def _moves(indices):
    return "".join("Y" if a < 0 else "X" if b < 0 else "D" for a, b in indices.T)


def _profile(rows, order, dtype=np.float64):
    bits = np.array([[x == "-" for x in row] for row in rows])
    counts = np.array(
        [[sum(row[i] == x for row in rows) for x in "AC"] for i in range(len(rows[0]))],
        dtype=dtype,
    )
    return counts, bits, order


def _tree(n):
    """Build an unbalanced guide with default row IDs."""
    tree = TreeNode(name="0")
    for i in range(1, n):
        tree = TreeNode(children=[tree, TreeNode(name=str(i))])
    return tree


class ProfileKernelTests(unittest.TestCase):
    def test_exhaustive_paths(self):
        # Includes boundary insertions, direction switches, zero costs, long runs,
        # and both fused numeric types. All scores are exactly representable.
        rng = np.random.default_rng(27)
        for dtype, p, q, gap, free, method in product(
            [np.float32, np.float64],
            range(1, 4),
            range(1, 4),
            [0, 2, (0, 2), (5, 0), (5, 2), (0.5, 0.25)],
            [False, True],
            ["full", "rolling"],
        ):
            scores = rng.choice([-20.0, -1.0, 0.0, 0.25, 2.0], (p, q)).astype(dtype)
            o, e = (0, gap) if np.isscalar(gap) else gap
            indices, score = _align_profiles(
                scores, dtype(o), dtype(e), free, method=method
            )
            expected = max(_score_path(x, scores, gap, free) for x in _paths(p, q))
            with self.subTest(dtype=dtype, p=p, q=q, gap=gap, free=free):
                self.assertEqual(score, expected)
                self.assertEqual(
                    _score_path(_moves(indices), scores, gap, free), expected
                )
                npt.assert_array_equal(indices[0, indices[0] >= 0], np.arange(p))
                npt.assert_array_equal(indices[1, indices[1] >= 0], np.arange(q))
                self.assertFalse((indices == -1).all(axis=0).any())

    def test_backend_agreement_and_workspace(self):
        rng = np.random.default_rng(32)
        for dtype, gap, free, atol in product(
            [np.float32, np.float64],
            [(0, 0.2), (1.1, 0.3), (3, 0)],
            [False, True],
            [0, 1e-5, 0.1],
        ):
            full, rolling = _ProfileWorkspace(), _ProfileWorkspace()
            for m, n in [(1, 5), (7, 2), (2, 9), (3, 3), (1, 1)]:
                scores = rng.normal(size=(m, n)).astype(dtype)
                costs = tuple(map(dtype, gap))
                a = _align_profiles(scores, *costs, free, full, "full", dtype(atol))
                b = _align_profiles(
                    scores, *costs, free, rolling, "rolling", dtype(atol)
                )
                self.assertEqual(a[1], b[1])
                npt.assert_array_equal(a[0], b[0])
                self.assertFalse(np.shares_memory(a[0], full.buffers["path"]))
            before = dict(rolling.buffers)
            _align_profiles(
                np.ones((1, 1), dtype=dtype),
                *costs,
                free,
                rolling,
                "rolling",
                dtype(atol),
            )
            for name, buffer in before.items():
                self.assertIs(buffer, rolling.buffers[name])
        ws = _ProfileWorkspace()
        a = ws.get("test", (3, 4), np.float32)
        b = ws.get("test", (2, 5), np.float32)
        self.assertTrue(np.shares_memory(a, b))
        self.assertTrue(b.flags.c_contiguous)
        self.assertEqual(ws.get("test", (1, 1), np.float64).dtype, np.float64)

    def test_linear_affine_equivalence(self):
        scores = np.array([[2.0, -3.0], [-1.0, 2.0], [0.0, 0.0]])
        for free in [False, True]:
            a = _align_profiles(scores, 0.0, 2.0, free)
            optimum = max(_score_path(x, scores, (0, 2), free) for x in _paths(3, 2))
            self.assertEqual(a[1], optimum)

    def test_direction_switch_and_ties(self):
        scores = np.array([[-20.0]])
        indices, score = _align_profiles(scores, 5.0, 2.0, False)
        self.assertEqual(score, -14.0)
        self.assertEqual(_moves(indices), "YX")  # terminal X wins
        for o in [0.0, 5.0]:
            indices, score = _align_profiles(scores, o, 2.0, True)
            self.assertEqual(score, 0.0)
            self.assertEqual(_moves(indices), "YX")
        for o in [0.0, 2.0]:
            indices, _ = _align_profiles(np.zeros((3, 3)), o, 0.0, True)
            self.assertEqual(_moves(indices), "YYYXXX")

    def test_small_difference_is_not_a_tie(self):
        for dtype in [np.float32, np.float64]:
            scores = np.array([[np.finfo(dtype).eps]], dtype=dtype)
            for o in [0.0, 1.0]:
                indices, score = _align_profiles(scores, dtype(o), dtype(0), True)
                self.assertEqual(_moves(indices), "D")
                self.assertEqual(score, scores[0, 0])

    def test_hand_scoring(self):
        scores = np.full((3, 2), 2.0)
        self.assertEqual(_score_path("XDD", scores, (5, 2), False), -3)
        self.assertEqual(_score_path("XDD", scores, (5, 2), True), 4)
        self.assertEqual(_score_path("DXD", scores, (5, 2), True), -3)
        for gap, cost in [(2, 6), ((5, 2), 11), ((5, 0), 5)]:
            self.assertEqual(
                _score_path("XXXDD", np.full((5, 2), 2.0), gap, False), 4 - cost
            )

    def test_profile_average_and_row_counts(self):
        matrix = np.array([[2.0, -1.0], [-1.0, 2.0]])
        a = _profile(["A", "A", "-", "-"], [0, 1, 2, 3])
        b = _profile(["A", "C"], [4, 5])
        # The isolated columns can contain gap-only rows here: real child profiles
        # have other columns with residues. Test the local averaging independently.
        c, bits, order = _merge_profiles(a, b, matrix, 10.0, 1.0, False)
        npt.assert_array_equal(c, [[3, 1]])
        npt.assert_array_equal(bits[:, 0], [0, 0, 1, 1, 0, 0])
        self.assertEqual(order, list(range(6)))
        score = (
            sum(
                0 if x == "-" or y == "-" else matrix["AC".index(x), "AC".index(y)]
                for x in "AA--"
                for y in "AC"
            )
            / 8
        )
        self.assertEqual(score, 0.25)
        a = _profile(["A"], [0])
        b = _profile(["C", "C", "C"], [1, 2, 3])
        c, _, _ = _merge_profiles(a, b, matrix, 10.0, 1.0, False)
        npt.assert_array_equal(c / 4, [[0.25, 0.75]])

    def test_merge_against_row_pair_oracle(self):
        matrix = np.array([[2.0, -1.0], [-1.0, 2.0]])
        for rows_a, rows_b, gap, free in product(
            [["A-C", "ACC"], ["AC-", "-CA", "A--"]],
            [["CA", "-A"], ["AAC", "C--"]],
            [(0, 2), (5, 2), (5, 0)],
            [False, True],
        ):
            r, s = len(rows_a), len(rows_b)
            p, q = len(rows_a[0]), len(rows_b[0])
            scores = np.array(
                [
                    [
                        sum(
                            0
                            if a[i] == "-" or b[j] == "-"
                            else matrix["AC".index(a[i]), "AC".index(b[j])]
                            for a in rows_a
                            for b in rows_b
                        )
                        / (r * s)
                        for j in range(q)
                    ]
                    for i in range(p)
                ]
            )
            candidates = []
            for path in _paths(p, q):
                i = j = 0
                merged = [""] * (r + s)
                for move in path:
                    for k, row in enumerate(rows_a):
                        merged[k] += "-" if move == "Y" else row[i]
                    for k, row in enumerate(rows_b, r):
                        merged[k] += "-" if move == "X" else row[j]
                    i += move != "Y"
                    j += move != "X"
                candidates.append((_score_path(path, scores, gap, free), merged))
            best = max(x[0] for x in candidates)
            a = _profile(rows_a, list(range(r)))
            b = _profile(rows_b, list(range(r, r + s)))
            counts, bits, _ = _merge_profiles(a, b, matrix, *gap, free)
            ungapped = [x.replace("-", "") for x in rows_a + rows_b]
            merged = AlignPath.from_bits(bits).to_aligned(ungapped)
            self.assertTrue(
                any(
                    abs(score - best) < 1e-12 and rows == merged
                    for score, rows in candidates
                )
            )
            npt.assert_array_equal(counts, _profile(merged, [])[0])
            for old, mask in [(rows_a, bits[:r]), (rows_b, bits[r:])]:
                npt.assert_array_equal(mask[:, ~mask.all(axis=0)], _profile(old, [])[1])

    def test_old_gaps_do_not_make_internal_insertions_free(self):
        # Profile A's second row ends at column 1, but an insertion into A after
        # column 1 is internal to the *profile*, even with free ends.
        a = _profile(["AC", "A-"], [0, 1])
        b = _profile(["AAC"], [2])
        matrix = np.array([[2.0, -20.0], [-20.0, 2.0]])
        scores = (a[0] / 2) @ matrix @ b[0].T
        self.assertEqual(_score_path("DYD", scores, (5, 2), True), -4)
        indices, score = _align_profiles(scores, 5.0, 2.0, True)
        self.assertGreater(score, -4)
        self.assertNotEqual(_moves(indices), "DYD")


class MultiAlignTests(unittest.TestCase):
    def test_automatic(self):
        seqs = ["ACGT", "AGT", "ACGT"]
        path = multi_align(seqs, free_ends=False)
        self.assertIs(type(path), AlignPath)
        self.assertEqual(path.to_aligned(seqs), ["ACGT", "A-GT", "ACGT"])
        self.assertEqual(align_score((path, seqs), free_ends=False), 6)
        for free, gap in product([False, True], [2, (5, 2), (5, 0)]):
            path = multi_align(seqs, gap_cost=gap, free_ends=free)
            self.assertEqual([s.replace("-", "") for s in path.to_aligned(seqs)], seqs)

    def test_types(self):
        seqs = ["ACGT", "AGT", "ACGT"]
        for cls in [str, Sequence, DNA, RNA, Protein]:
            data = [cls(x.replace("T", "U") if cls is RNA else x) for x in seqs]
            path = multi_align(data, free_ends=False)
            npt.assert_array_equal(
                path.to_bits(), [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
            )
        for data in [
            [[1, 2, 3], [1, 3], [1, 2, 3]],
            [["cat", "dog", "bird"], ["cat", "bird"], ["cat", "dog", "bird"]],
            ["αβγ", "αγ", "αβγ"],
        ]:
            path = multi_align(iter(data), guide_tree=_tree(3), free_ends=False)
            npt.assert_array_equal(path.to_bits()[1], [0, 1, 0])
        matrix = SubstitutionMatrix(
            "ACGT", [[2, -1, -1, -1], [-1, 2, -1, -1], [-1, -1, 2, -1], [-1, -1, -1, 2]]
        )
        path = multi_align(seqs, sub_score=matrix, free_ends=False)
        self.assertEqual(path.to_aligned(seqs)[1], "A-GT")
        path = multi_align(
            [Protein("MKT"), Protein("MT"), Protein("MKT")],
            sub_score="BLOSUM62",
            guide_tree=_tree(3),
            free_ends=False,
        )
        self.assertEqual(path.shape[0], 3)

    def test_tree_ids_order_and_no_mutation(self):
        seqs = ["ACGT", "AGT", "ACGT", "ACT"]
        ids = ["a", "b", "c", "d"]
        tree = TreeNode.read(["((d,b),(c,a));"])
        before = str(tree)
        path = multi_align(seqs, ids=iter(ids), guide_tree=tree, free_ends=False)
        self.assertEqual([x.replace("-", "") for x in path.to_aligned(seqs)], seqs)
        self.assertEqual(str(tree), before)
        self.assertEqual(ids, ["a", "b", "c", "d"])
        npt.assert_array_equal(path.starts, [0, 0, 0, 0])

    def test_metadata_ids(self):
        from skbio.io import read

        seqs = list(
            read(
                StringIO(">a\nACGT\n>b\nAGT\n>c\nACGT\n"),
                format="fasta",
                constructor=DNA,
            )
        )
        tree = TreeNode.read(["((c,a),b);"])
        metadata = [dict(x.metadata) for x in seqs]
        path = multi_align(seqs, guide_tree=tree, free_ends=False)
        self.assertEqual(path.to_aligned(seqs), ["ACGT", "A-GT", "ACGT"])
        self.assertEqual([x.metadata for x in seqs], metadata)
        multi_align(seqs, free_ends=False)
        # No metadata IDs: use positional IDs, including other metadata fields.
        plain = [DNA("AC", metadata={"description": "test"}), DNA("AC"), DNA("AC")]
        multi_align(plain, guide_tree=_tree(3))
        for values in [["a", "a", "c"], ["a", None, "c"], ["a", 1, "c"]]:
            data = [DNA("AC", metadata={"id": x}) for x in values]
            with self.assertRaisesRegex(ValueError, "unique strings"):
                multi_align(data)
            multi_align(data, ids=["0", "1", "2"], guide_tree=_tree(3))
        partial = [DNA("AC", metadata={"id": "a"}), DNA("AC"), DNA("AC")]
        with self.assertRaisesRegex(ValueError, "every sequence or none"):
            multi_align(partial)
        multi_align(partial, ids=["0", "1", "2"], guide_tree=_tree(3))
        # Presence means key membership, not truthiness: an empty string is an ID.
        multi_align([DNA("AC", metadata={"id": ""}), DNA("AC", metadata={"id": "b"})])

    def test_two_sequences(self):
        for seqs, gap, free in product(
            [["AAAA", "CCCC"], ["ACGT", "AGT"]], [2, (5, 2)], [False, True]
        ):
            expected = pair_align(*seqs, gap_cost=gap, free_ends=free).score
            path = multi_align(seqs, gap_cost=gap, free_ends=free)
            self.assertIs(type(path), AlignPath)
            self.assertEqual(
                align_score((path, seqs), gap_cost=gap, free_ends=free), expected
            )
        self.assertEqual(multi_align(["A", "C"], guide_tree=_tree(2)).shape[0], 2)

    def test_duplicates_and_packing(self):
        for n in [3, 7, 8, 9, 16, 17]:
            seqs = ["AAA"] * n
            with self.assertRaisesRegex(ValueError, "normalization is not positive"):
                multi_align(seqs)
            path = multi_align(seqs, guide_tree=_tree(n))
            self.assertEqual(path.to_aligned(seqs), seqs)
            seqs = ["ACGT" if i % 2 else "AGT" for i in range(n)]
            path = multi_align(seqs, guide_tree=_tree(n), free_ends=False)
            self.assertEqual([x.replace("-", "") for x in path.to_aligned(seqs)], seqs)
            self.assertFalse(path.to_bits().all(axis=0).any())

    def test_unrelated_sequences_with_tree(self):
        seqs = ["AAAA", "CCCC", "GGGG"]
        with self.assertRaisesRegex(ValueError, "supply a guide tree"):
            multi_align(seqs, free_ends=False)
        path = multi_align(seqs, guide_tree=_tree(3))
        self.assertEqual([s.replace("-", "") for s in path.to_aligned(seqs)], seqs)

    def test_backends(self):
        seqs = ["ACGT", "AGT", "ACGT", "ACT"]
        for gap, free, atol in product([2, (5, 2)], [True, False], [0, 1e-5]):
            a = multi_align(
                seqs, gap_cost=gap, free_ends=free, atol=atol, method="full"
            )
            b = multi_align(
                seqs, gap_cost=gap, free_ends=free, atol=atol, method="rolling"
            )
            self.assertEqual(a.to_aligned(seqs), b.to_aligned(seqs))
        with self.assertRaisesRegex(ValueError, "method"):
            multi_align(seqs, method="invalid")
        for atol in [-1, np.nan, np.inf, (0, 1)]:
            with self.assertRaisesRegex(ValueError, "atol"):
                multi_align(seqs, atol=atol)
        with np.errstate(over="ignore"):
            with self.assertRaisesRegex(ValueError, "scoring dtype"):
                multi_align(seqs, atol=1e100)

    def test_invalid_inputs(self):
        for seqs in [[], ["A"]]:
            with self.assertRaisesRegex(ValueError, "At least two"):
                multi_align(seqs)
        for seqs in [["", "A"], ["A", ""], ["", "", ""]]:
            with self.assertRaisesRegex(ValueError, "length of zero"):
                multi_align(seqs)
        for value in ["A-", "A.", Sequence("A-"), DNA("A."), RNA("A-")]:
            with self.assertRaisesRegex(ValueError, "ungapped"):
                multi_align([value, value])
        with self.assertRaisesRegex(TypeError, "decoded"):
            multi_align([b"AC", b"AC"])
        for free in [0, "yes", (True, False), None]:
            with self.assertRaisesRegex(TypeError, "Boolean"):
                multi_align(["A", "A"], free_ends=free)
        multi_align(["A", "A"], free_ends=np.bool_(True))
        for gap in [-1, np.inf, np.nan, (1, -1), (-1, 1), (np.inf, 1)]:
            with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
                multi_align(["A", "A"], gap_cost=gap)
        for score in [
            (np.inf, -1),
            (1, np.nan),
            SubstitutionMatrix("AC", [[1, 0], [-1, 1]]),
        ]:
            with self.assertRaisesRegex(ValueError, "finite and symmetric"):
                multi_align(["AC", "AC"], sub_score=score)
        with self.assertRaises(TypeError):
            multi_align(["AC", DNA("AC")])
        with self.assertRaises(ValueError):
            multi_align(
                ["AC", "AT"], sub_score=SubstitutionMatrix("AC", [[1, -1], [-1, 1]])
            )
        for ids in [[], ["x"], ["x", "x"], [0, 1], ["x", []]]:
            with self.assertRaisesRegex(ValueError, "unique strings"):
                multi_align(["A", "A"], ids=ids)
        with self.assertRaisesRegex(TypeError, "TreeNode"):
            multi_align(["A", "A"], guide_tree="(0,1);")
        for newick in ["(0,0);", "(0,2);", "0;", "((0,1),2);"]:
            with self.assertRaisesRegex(ValueError, "tips must match"):
                multi_align(["A", "A"], guide_tree=TreeNode.read([newick]))
        for newick in ["(0,1,2);", "((0), (1,2));", "(((0,1),2));"]:
            with self.assertRaisesRegex(ValueError, "must be binary"):
                multi_align(["A"] * 3, guide_tree=TreeNode.read([newick]))


class PairWorkspaceTests(unittest.TestCase):
    def test_pairwise_agreement(self):
        # Changing both dimensions catches stale boundaries and row-stride reuse.
        rng = np.random.default_rng(81)
        for dtype, gap, free, atol in product(
            [np.float32, np.float64],
            [(0, 0.2), (1.1, 0.3), (3, 0)],
            [False, True],
            [0, 1e-5, 0.1],
        ):
            matrix = SubstitutionMatrix(
                "ACGT",
                np.array(
                    [[2.3 if i == j else -1.1 for j in range(4)] for i in range(4)],
                    dtype=dtype,
                ),
            )
            workspace = _ProfileWorkspace()
            costs = tuple(map(dtype, gap))
            for m, n in [(1, 5), (7, 2), (2, 9), (3, 3), (1, 1)]:
                seqs = ["".join(rng.choice(list("ACGT"), length)) for length in (m, n)]
                encoded, submat, _ = encode_sequences(seqs, matrix)
                moves, score = _align_pair(
                    submat[encoded[0]], encoded[1], *costs, free, workspace, dtype(atol)
                )
                path = _encode_path(moves, 0, m, 0, n)
                expected = pair_align(
                    *seqs, sub_score=matrix, gap_cost=gap, free_ends=free, atol=atol
                )
                self.assertEqual(score, expected.score)
                npt.assert_array_equal(path.to_bits(), expected.paths[0].to_bits())
                _, score_only = _align_pair(
                    submat[encoded[0]],
                    encoded[1],
                    *costs,
                    free,
                    workspace,
                    dtype(atol),
                    traceback=False,
                )
                self.assertEqual(score_only, score)
                self.assertEqual(workspace.buffers["dp0"].dtype, dtype)

    def test_rolling_helper(self):
        for dtype, gap, free in product(
            [np.float32, np.float64], [(0, 0.3), (1.1, 0.3)], [False, True]
        ):
            scores = np.array([[2, -1, 2], [-1, 2, -1]], dtype=dtype)
            workspace = _ProfileWorkspace()
            args = (*map(dtype, gap), free, workspace, dtype(1e-5))
            moves, score = _align_pair_roll(scores, *args)
            self.assertTrue(np.shares_memory(moves, workspace.buffers["path"]))
            saved = moves.copy()
            full, expected = _align_pair(scores, np.arange(3), *args)
            npt.assert_array_equal(saved, full)
            self.assertEqual(score, expected)

    def test_self_score_without_traceback(self):
        # A shifted self-alignment beats the ungapped diagonal for this matrix.
        matrix = np.array([[-1, 2], [2, -1]], dtype=np.float32)
        seq = np.array([0, 1], dtype=np.intp)
        workspace = _ProfileWorkspace()
        moves, score = _align_pair(
            matrix[seq],
            seq,
            np.float32(0),
            np.float32(0),
            False,
            workspace,
            np.float32(0),
            traceback=False,
        )
        self.assertIsNone(moves)
        self.assertEqual(score, 2)
        self.assertNotIn("path", workspace.buffers)

    def test_guide_indices(self):
        tree = TreeNode.read(["((d,b),(a,c));"])
        npt.assert_array_equal(
            _guide_merges(tree, list("abcd")), [[3, 1], [0, 2], [4, 5]]
        )


class MultiDistanceTests(unittest.TestCase):
    def test_hand_calculation(self):
        # AC vs AG: S=0, S_max=2, S_rand=(1-1-1-1)/2=-1, d=ln(3).
        matrix = np.full((3, 3), -1.0, dtype=np.float64)
        np.fill_diagonal(matrix, 1.0)
        encoded = [np.array([0, 1]), np.array([0, 2]), np.array([0, 1])]
        dm = _multi_distances(
            encoded, matrix, 0.0, 2.0, False, list("abc"), _ProfileWorkspace(), 0.0
        )
        npt.assert_allclose(dm, [np.log(3), 0, np.log(3)])
        self.assertEqual(dm.dtype, np.float64)

    def test_gap_statistics(self):
        # AC vs A: a single terminal gap, L=2, expected substitution sum / L=0.
        # Penalized: S=-1, S_rand=-2, S_max=1.5 -> d=-ln(1/3.5).
        # Free: S=1, S_rand=0, S_max=1.5 -> d=-ln(1/1.5).
        matrix = np.array([[1.0, -1.0], [-1.0, 1.0]])
        encoded = [np.array([0, 1]), np.array([0])]
        for free, expected in [(False, np.log(3.5)), (True, np.log(1.5))]:
            dm = _multi_distances(
                encoded, matrix, 1.0, 1.0, free, ["0", "1"], _ProfileWorkspace(), 0.0
            )
            self.assertAlmostEqual(dm[0], expected)

    def test_condensed_distances_and_guide(self):
        from scipy.spatial.distance import squareform
        from skbio import DistanceMatrix
        from skbio.tree import upgma

        seqs = ["ACGTACGT", "ACGTCGT", "ACGTACCT", "ACGTACG", "ACGTACGT"]
        ids = list("abcde")
        for dtype, gap, free in product(
            [np.float32, np.float64], [(0, 0.3), (1.1, 0.3)], [False, True]
        ):
            submat = SubstitutionMatrix(
                "ACGT",
                np.array(
                    [[2.3 if i == j else -1.1 for j in range(4)] for i in range(4)],
                    dtype=dtype,
                ),
            )
            encoded, matrix, _ = encode_sequences(seqs, submat)
            costs = tuple(map(dtype, gap))
            observed = _multi_distances(
                encoded, matrix, *costs, free, ids, _ProfileWorkspace(), dtype(1e-5)
            )
            kwargs = dict(sub_score=submat, gap_cost=gap, free_ends=free)
            selfs = [pair_align(s, s, max_paths=0, **kwargs).score for s in seqs]
            expected = np.zeros((len(seqs), len(seqs)))
            for i in range(len(seqs)):
                for j in range(i + 1, len(seqs)):
                    pair = pair_align(seqs[i], seqs[j], **kwargs)
                    # Explicit residue pairs avoid reusing the count-matrix formula.
                    composition = sum(
                        float(matrix[a, b]) for a in encoded[i] for b in encoded[j]
                    )
                    expected[i, j] = expected[j, i] = _fd_dist(
                        pair.score,
                        pair.paths[0],
                        composition,
                        (selfs[i] + selfs[j]) / 2,
                        *map(float, costs),
                        free,
                        8 * np.finfo(dtype).eps,
                    )
            npt.assert_allclose(observed, squareform(expected), rtol=1e-14, atol=1e-14)
            old_tree = upgma(DistanceMatrix(expected, ids))
            a = multi_align(seqs, ids=ids, **kwargs)
            b = multi_align(seqs, ids=ids, guide_tree=old_tree, **kwargs)
            npt.assert_array_equal(a.to_bits(), b.to_bits())

    def test_distance_domains(self):
        cases = [
            (["AAAA", "CCCC", "AAAA"], (1, -1), 2, "random baseline"),
            (["AC", "CA", "AC"], (0, 0), 2, "not positive"),
            (
                ["AC", "CA", "AC"],
                SubstitutionMatrix("AC", [[-10, 2], [2, -10]]),
                0,
                "exceeds",
            ),
        ]
        for seqs, matrix, gap, message in cases:
            with self.assertRaisesRegex(ValueError, message):
                multi_align(seqs, sub_score=matrix, gap_cost=gap, free_ends=False)
            # The same data remain alignable when a guide is supplied.
            path = multi_align(
                seqs,
                sub_score=matrix,
                gap_cost=gap,
                free_ends=False,
                guide_tree=_tree(3),
            )
            self.assertEqual([x.replace("-", "") for x in path.to_aligned(seqs)], seqs)
        path = multi_align(["AAAA", "CCCC", "GGGG"], free_ends=True)
        self.assertEqual(path.shape[0], 3)

    def test_roundoff(self):
        path = pair_align("A", "A", free_ends=False).paths[0]
        for dtype in [np.float32, np.float64]:
            eps = np.finfo(dtype).eps
            # Direct numerical inputs isolate the documented boundary at one.
            for ratio in [1, 1 + 2 * eps, 1 + 8 * eps]:
                self.assertEqual(_fd_dist(ratio, path, 0, 1, 0, 0, False, 8 * eps), 0)
            self.assertAlmostEqual(
                _fd_dist(0.5, path, 0, 1, 0, 0, False, 8 * eps), np.log(2)
            )
            with self.assertRaisesRegex(ValueError, "exceeds"):
                _fd_dist(1 + 16 * eps, path, 0, 1, 0, 0, False, 8 * eps)
        for score in [np.inf, np.nan]:
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                _fd_dist(score, path, 0, 1, 0, 0, False, 0)

    def test_custom_self_alignment(self):
        # Non-diagonal self-alignments score 2, rather than the diagonal's -20.
        # An incorrect diagonal shortcut would make normalization nonpositive.
        sm = SubstitutionMatrix("AC", [[-10.0, 2.0], [2.0, -10.0]])
        self.assertEqual(
            pair_align("AC", "AC", sub_score=sm, gap_cost=0, free_ends=False).score, 2
        )
        with self.assertRaisesRegex(ValueError, "exceeds the self-score"):
            multi_align(["AC", "CA", "AC"], sub_score=sm, gap_cost=0, free_ends=False)


if __name__ == "__main__":
    unittest.main()
