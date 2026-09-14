# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import lru_cache
from itertools import product
import unittest

import numpy as np
import numpy.testing as npt

from skbio import DNA, RNA, Protein, Sequence, TabularMSA, TreeNode, SubstitutionMatrix
from skbio.io import read as sk_read
from skbio.util._array import ArrayWorkspace
from skbio.alignment import AlignPath, multi_align, pair_align, align_score
from skbio.alignment._utils import encode_sequences
from skbio.alignment._pair import _encode_path
from skbio.alignment._multi import (
    _align_profiles,
    _merge_profiles,
    _multi_distances,
    _fd_dist,
    _align_pair,
    _align_pair_roll,
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


class MultiAlignTests(unittest.TestCase):
    def test_multi_align_nucl(self):
        """A simple case of three nucleotide sequences."""
        seqs = list(map(DNA, [
            "CAGCTATATATCGCTACG",
            "CTGCTTATATCCCTAGG",
            "AAGCTATACATCCAACATG",
        ]))

        # Default parameters (match = 1, mismatch = -1, gap = -2, free ends)
        path = multi_align(seqs)
        self.assertIs(type(path), AlignPath)
        obs = path.to_aligned(seqs)
        exp = [
            "CAGCTATATATCGCTACG--",
            "CTGCT-TATATCCCTAGG--",
            "AAGCTATACATC-CAACATG",
        ]
        self.assertListEqual(obs, exp)

        # Construct tabular MSA
        msa = TabularMSA.from_path_seqs(path, seqs)
        obs = [str(seq) for seq in msa]
        self.assertListEqual(obs, exp)

        # Restore original sequences by removing gaps
        obs = [str(seq.degap()) for seq in msa]
        exp = [str(seq) for seq in seqs]
        self.assertListEqual(obs, exp)

        # Calculate alignment score (sum-of-pairs; SP)
        # NOTE: `align_score` has the same default parameters.
        obs = align_score((path, seqs))
        exp = 16.0
        self.assertEqual(obs, exp)

    def test_multi_align_prot(self):
        """Align protein sequences."""
        seqs = [Protein("MKT"), Protein("MT"), Protein("MKT")]
        # params = dict(sub_score="BLOSUM62", free_ends=False)
        path = multi_align(
            seqs,
            sub_score="BLOSUM62",
            guide_tree=_tree(3),
            free_ends=False,
        )
        self.assertEqual(path.shape[0], 3)

    def test_multi_align_p53(self):
        """Align P53 transactivation motif sequences (protein)."""
        # - human:  NP_000537.3:6-30    [Homo sapiens]
        # - monkey: NP_001040616.1:6-30 [Macaca mulatta]
        # - mouse:  NP_035770.2:9-31    [Mus musculus]
        # - rat:    NP_112251.2:6-30    [Rattus norvegicus]
        # - dog:    NP_001376147.1:6-30 [Canis lupus familiaris]
        # - pig:    NP_998989.3:6-30    [Sus scrofa]
        # - cat:    NP_001009294.1:6-30 [Felis catus]
        # - horse:  NP_001189334.1:7-30 [Equus caballus]
        # - frog:   NP_001001903.1:5-28 [Xenopus tropicalis]
        # - trout:  NP_001118164.1:4-23 [Oncorhynchus mykiss]
        p53_faa = "\n".join((
            ">human",
            "SDPSVEPPLSQETFSDLWKLLPENN",
            ">monkey",
            "SDPSIEPPLSQETFSDLWKLLPENN",
            ">mouse",
            "SDISLELPLSQETFSGLWKLLPP",
            ">rat",
            "SDMSIELPLSQETFSCLWKLLPPDD",
            ">dog",
            "SELNIDPPLSQETFSELWNLLPENN",
            ">pig",
            "SELGVEPPLSQETFSDLWKLLPENN",
            ">cat",
            "LELTIEPPLSQETFSELWNLLPENN",
            ">horse",
            "ELGIEPPLSQETFSDLWKLLPENN",
            ">frog",
            "SETGMEPPLSQETFEDLWSLLPDP",
            ">trout",
            "LAENVSLPLSQESFEDLWKM",
        ))
        seqs = list(sk_read([p53_faa], format="fasta", constructor=Protein))
        params = dict(sub_score="BLOSUM62", gap_cost=(11, 1))  # BLASTP
        path = multi_align(seqs, **params)
        obs = path.to_aligned(seqs)
        exp = [
            "SDPSVEPPLSQETFSDLWKLLPENN",
            "SDPSIEPPLSQETFSDLWKLLPENN",
            "SDISLELPLSQETFSGLWKLLPP--",
            "SDMSIELPLSQETFSCLWKLLPPDD",
            "SELNIDPPLSQETFSELWNLLPENN",
            "SELGVEPPLSQETFSDLWKLLPENN",
            "LELTIEPPLSQETFSELWNLLPENN",
            "-ELGIEPPLSQETFSDLWKLLPENN",
            "SETGMEPPLSQETFEDLWSLLPDP-",
            "LAENVSLPLSQESFEDLWKM-----",
        ]
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score((path, seqs), **params), 3971.0)

        msa = TabularMSA.from_path_seqs(path, seqs)
        obs = [str(seq.degap()) for seq in msa]
        exp = [str(seq) for seq in seqs]
        self.assertListEqual(obs, exp)

        # TODO: Let `from_path_seqs` incorporate sequence IDs.
        # obs = [seq.metadata["id"] for seq in msa]
        # exp = [seq.metadata["id"] for seq in seqs]
        # self.assertListEqual(obs, exp)

    def test_multi_align_trna(self):
        """Align human mitochondrial tRNA sequences (nucleotide)."""
        # Source: NC_012920.1
        # - tRNA-Trp : 5512-5579
        # - tRNA-Ala : c5655-5587
        # - tRNA-Asn : c5729-5657
        # - tRNA-Cys : c5826-5761
        # - tRNA-Tyr : c5891-5826
        trna_frn = "\n".join((
            ">Trp",
            "AGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTGCAATACTTAATTTCTG",
            ">Ala",
            "AAGGGCTTAGCTTAATTAAAGTGGCTGATTTGCGTTCAGTTGATGCAGAGTGGGGTTTTGCAGTCCTTA",
            ">Asn",
            "TAGATTGAAGCCAGTTGATTAGGGTGCTTAGCTGTTAACTAAGTGTTTGTGGGTTTAAGTCCCATTGGTCTAG",
            ">Cys",
            "AGCTCCGAGGTGATTTTCATATTGAATTGCAAATTCGAAGAAGCAGCTTCAAACCTGCCGGGGCTT",
            ">Tyr",
            "GGTAAAATGGCTGAGTGAAGCATTGGACTGTAAATCTAAAGACAGGGGTTAGGCCTCTTTTTACCA",
        ))
        seqs = list(sk_read([trna_frn], format='fasta', constructor=DNA))

        # NOTE: tRNA sequences are too diverse and `free_ends=True` (default mode) will
        # produce a poorly overlapped alignment.
        params = dict(sub_score=(2, -3), gap_cost=(5, 2), free_ends=False)  # BLASTN
        path = multi_align(seqs, **params)
        obs = path.to_aligned(seqs)
        exp = [
            "AGAAATTTAGGTTAAATACAGACCAAGA-----GCCTTCA------AAGCCC---TCAGT"
            "AAGTT--GCAATACTTAATTTCT-G",
            "AAGGGCTTAGCTTAATTAAAGTGGCTGATTT--GCGTT---CAGTTGATGCAGAGT--G-"
            "GGGTTTTGCAGTCC--TT-----A-",
            "TAGATTGAAGCC-------AGT---TGATTAGGGTGCTTAGCTGTTAACTAAGTGTTTGT"
            "GGGTTTA--AGTCCCATTGGTCTAG",
            "--------AGCTCCGA-------GGTGATTTTCATATTGAATTGCAAATTCGAAGA-AGC"
            "AGCTT---CAAACCTGCCGGGGCTT",
            "--GG------------TAAAATGGCTGAGTGAAGCATTGGACTGTAAATCTAAAGACAG-"
            "GGGTTAGGCC-TCT--TTTTACCA-",
        ]
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score((path, seqs), **params), -1183.0)

        obs = [seq.replace("-", "") for seq in obs]
        exp = [str(seq) for seq in seqs]
        self.assertListEqual(obs, exp)

    def test_multi_align_params(self):
        """Test that common parameter settings work (don't break)."""
        seqs = ["ACGT", "AGT", "ACGT"]
        path = multi_align(seqs, free_ends=False)
        self.assertIs(type(path), AlignPath)
        self.assertEqual(path.to_aligned(seqs), ["ACGT", "A-GT", "ACGT"])
        self.assertEqual(align_score((path, seqs), free_ends=False), 6)
        for free, gap in product([False, True], [2, (5, 2), (5, 0)]):
            path = multi_align(seqs, gap_cost=gap, free_ends=free)
            self.assertEqual([s.replace("-", "") for s in path.to_aligned(seqs)], seqs)

    def test_multi_align_input(self):
        """Test various types of input sequences and parameters."""
        seqs = ["ACGT", "AGT", "ACGT"]

        # input sequence class
        for cls in [str, Sequence, DNA, RNA, Protein]:
            data = [cls(x.replace("T", "U") if cls is RNA else x) for x in seqs]
            path = multi_align(data, free_ends=False)
            npt.assert_array_equal(
                path.to_bits(), [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
            )

        # input sequence type
        for data in [
            [[1, 2, 3], [1, 3], [1, 2, 3]],  # integers
            [["cat", "dog", "bird"], ["cat", "bird"], ["cat", "dog", "bird"]],  # words
            ["αβγ", "αγ", "αβγ"],  # unicode
        ]:
            path = multi_align(iter(data), guide_tree=_tree(3), free_ends=False)
            npt.assert_array_equal(path.to_bits()[1], [0, 1, 0])

        # custom substitution matrix
        submat = SubstitutionMatrix("ACGT", [
            [2, -1, -1, -1],
            [-1, 2, -1, -1],
            [-1, -1, 2, -1],
            [-1, -1, -1, 2],
        ])
        path = multi_align(seqs, sub_score=submat, free_ends=False)
        self.assertEqual(path.to_aligned(seqs)[1], "A-GT")

    def test_multi_align_two(self):
        """Test only two sequences (will fallback to pairwise alignment)."""
        seq1 = DNA("GAATTC")
        seq2 = DNA("AGATCT")
        obs = multi_align((seq1, seq2))

        # Output type is AlignPath nor PairAlignPath
        self.assertIs(type(obs), AlignPath)

        # Path matches `pair_align`
        res = pair_align(seq1, seq2)
        exp = res.paths[0]
        npt.assert_array_equal(obs.to_bits(), exp.to_bits())

        # Score matches `pair_align`
        obs = align_score((obs, (seq1, seq2)))
        exp = res.score
        self.assertEqual(obs, exp)

        # Custom parameters
        params = dict(sub_score=(5, -4), gap_cost=(5, 2))
        obs = multi_align((seq1, seq2), **params)
        exp = pair_align(seq1, seq2, **params).paths[0]
        npt.assert_array_equal(obs.to_bits(), exp.to_bits())

    def test_multi_align_tree(self):
        """Test guide tree's impact on alignment."""
        seqs = [DNA("CATTAACGT"),
                DNA("CGTTACGGT"),
                DNA("AGTTAACGG")]

        # Auto-determine merge order
        obs = multi_align(seqs).to_aligned(seqs)
        exp = ["CA-TTAACGT-",
               "-CGTTA-CGGT",
               "-AGTTAACGG-"]
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score(obs), 7.0)

        # Merging seqs 0 and 2 first produces the same output.
        tree = TreeNode.read(["((0,2),1);"])
        obs = multi_align(seqs, guide_tree=tree).to_aligned(seqs)
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score(obs), 7.0)

        # Merge seqs 0 and 1 first produces a less-aligned alignment, likely because
        # they are less similar, although the alignment score is the same.
        tree = TreeNode.read(["((0,1),2);"])
        obs = multi_align(seqs, guide_tree=tree).to_aligned(seqs)
        exp = ["CATTAACGT------",
               "------CGTTACGGT",
               "AGTTAACGG------"]
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score(obs), 7.0)

        # Merge seqs 1 and 2 first produces an alternative and even better alignment.
        tree = TreeNode.read(["((1,2),0);"])
        obs = multi_align(seqs, guide_tree=tree).to_aligned(seqs)
        exp = ["CATTAACGT-",
               "CGTTA-CGGT",
               "AGTTAACGG-"]
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score(obs), 9.0)

        msg = "Tree must be strictly bifurcating."
        tree = TreeNode.read(["(0,1,2);"])
        with self.assertRaises(ValueError) as cm:
            multi_align(seqs, guide_tree=tree)
        self.assertEqual(str(cm.exception), msg)

        # Supply sequence IDs
        for id_, seq in zip("abc", seqs):
            seq.metadata["id"] = id_
        tree = TreeNode.read(["((b,c),a);"])
        obs = multi_align(seqs, guide_tree=tree).to_aligned(seqs)
        self.assertListEqual(obs, exp)

    def test_multi_align_ids(self):
        """Test of custom sequence IDs."""
        # FASTA parsing automatically assign sequence IDs.
        fasta = ">a\nACGT\n>b\nAGT\n>c\nACGT\n"
        seqs = list(sk_read([fasta], format="fasta", constructor=DNA))

        # Confirm sequence IDs are assigned.
        obs = [seq.metadata["id"] for seq in seqs]
        exp = list("abc")
        self.assertListEqual(obs, exp)

        # Sequence IDs are not used.
        obs = multi_align(seqs, free_ends=False).to_aligned(seqs)
        exp = ["ACGT", "A-GT", "ACGT"]
        self.assertListEqual(obs, exp)

        # Sequence IDs match guide tree.
        tree = TreeNode.read(["((c,a),b);"])
        obs = multi_align(seqs, guide_tree=tree, free_ends=False).to_aligned(seqs)
        self.assertListEqual(obs, exp)

        # Sequence IDs don't match guide tree
        nwks = ("(a,b);",          # fewer tips
                "((a,b),(b,c));",  # more tips (with duplicate)
                "((a,b),x);",      # one tip not matching any sequence ID
                "((d,e),f);",      # all tips not match sequence IDs
                "((0,1),2);")      # numbers (invalid when sequences have IDs)
        for nwk in nwks:
            tree = TreeNode.read([nwk])
            self.assertRaises(ValueError, multi_align, seqs, guide_tree=tree)

        # TODO: If tree is not supplied then this is fine.
        msg = "Sequence IDs must be unique."
        tree = TreeNode.read(["((a,b),c);"])
        seqs[1].metadata["id"] = "a"
        with self.assertRaises(ValueError) as cm:
            multi_align(seqs, guide_tree=tree)
        self.assertEqual(str(cm.exception), msg)

        msg = "Metadata 'id' must be present in every sequence or none."
        del seqs[1].metadata["id"]
        with self.assertRaises(ValueError) as cm:
            multi_align(seqs, guide_tree=tree)
        self.assertEqual(str(cm.exception), msg)

        msg = "Tip name 'a' is absent from taxa."
        del seqs[0].metadata["id"]
        del seqs[2].metadata["id"]
        with self.assertRaises(ValueError) as cm:
            multi_align(seqs, guide_tree=tree)
        self.assertEqual(str(cm.exception), msg)

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

    # def test_invalid_inputs(self):
    #     for seqs in [[], ["A"]]:
    #         with self.assertRaisesRegex(ValueError, "At least two"):
    #             multi_align(seqs)
    #     for seqs in [["", "A"], ["A", ""], ["", "", ""]]:
    #         with self.assertRaisesRegex(ValueError, "length of zero"):
    #             multi_align(seqs)
    #     for value in ["A-", "A.", Sequence("A-"), DNA("A."), RNA("A-")]:
    #         with self.assertRaisesRegex(ValueError, "ungapped"):
    #             multi_align([value, value])
    #     with self.assertRaisesRegex(TypeError, "decoded"):
    #         multi_align([b"AC", b"AC"])
    #     for free in [0, "yes", (True, False), None]:
    #         with self.assertRaisesRegex(TypeError, "Boolean"):
    #             multi_align(["A", "A"], free_ends=free)
    #     multi_align(["A", "A"], free_ends=np.bool_(True))
    #     for gap in [-1, np.inf, np.nan, (1, -1), (-1, 1), (np.inf, 1)]:
    #         with self.assertRaisesRegex(ValueError, "finite and nonnegative"):
    #             multi_align(["A", "A"], gap_cost=gap)
    #     for score in [
    #         (np.inf, -1),
    #         (1, np.nan),
    #         SubstitutionMatrix("AC", [[1, 0], [-1, 1]]),
    #     ]:
    #         with self.assertRaisesRegex(ValueError, "finite and symmetric"):
    #             multi_align(["AC", "AC"], sub_score=score)
    #     with self.assertRaises(TypeError):
    #         multi_align(["AC", DNA("AC")])
    #     with self.assertRaises(ValueError):
    #         multi_align(
    #             ["AC", "AT"], sub_score=SubstitutionMatrix("AC", [[1, -1], [-1, 1]])
    #         )
    #     for ids in [[], ["x"], ["x", "x"], [0, 1], ["x", []]]:
    #         with self.assertRaisesRegex(ValueError, "unique strings"):
    #             multi_align(["A", "A"], guide_tree=_tree(2), ids=ids)
    #     with self.assertRaisesRegex(TypeError, "TreeNode"):
    #         multi_align(["A", "A"], guide_tree="(0,1);")
    #     for newick in ["(0,0);", "(0,2);", "0;", "((0,1),2);"]:
    #         with self.assertRaisesRegex(ValueError, "tips must match"):
    #             multi_align(["A", "A"], guide_tree=TreeNode.read([newick]))
    #     for newick in ["(0,1,2);", "((0), (1,2));", "(((0,1),2));"]:
    #         with self.assertRaisesRegex(ValueError, "must be binary"):
    #             multi_align(["A"] * 3, guide_tree=TreeNode.read([newick]))


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
            full, rolling = ArrayWorkspace(), ArrayWorkspace()
            for m, n in [(1, 5), (7, 2), (2, 9), (3, 3), (1, 1)]:
                scores = rng.normal(size=(m, n)).astype(dtype)
                costs = tuple(map(dtype, gap))
                a = _align_profiles(scores, *costs, free, full, "full", dtype(atol))
                b = _align_profiles(
                    scores, *costs, free, rolling, "rolling", dtype(atol)
                )
                self.assertEqual(a[1], b[1])
                npt.assert_array_equal(a[0], b[0])
                self.assertFalse(np.shares_memory(a[0], full.arrays["path"]))
            before = dict(rolling.arrays)
            _align_profiles(
                np.ones((1, 1), dtype=dtype),
                *costs,
                free,
                rolling,
                "rolling",
                dtype(atol),
            )
            for name, buffer in before.items():
                self.assertIs(buffer, rolling.arrays[name])
        ws = ArrayWorkspace()
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
            workspace = ArrayWorkspace()
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
                self.assertEqual(workspace.arrays["dp0"].dtype, dtype)

    def test_rolling_helper(self):
        for dtype, gap, free in product(
            [np.float32, np.float64], [(0, 0.3), (1.1, 0.3)], [False, True]
        ):
            scores = np.array([[2, -1, 2], [-1, 2, -1]], dtype=dtype)
            workspace = ArrayWorkspace()
            args = (*map(dtype, gap), free, workspace, dtype(1e-5))
            moves, score = _align_pair_roll(scores, *args)
            self.assertTrue(np.shares_memory(moves, workspace.arrays["path"]))
            saved = moves.copy()
            full, expected = _align_pair(scores, np.arange(3), *args)
            npt.assert_array_equal(saved, full)
            self.assertEqual(score, expected)

    def test_self_score_without_traceback(self):
        # A shifted self-alignment beats the ungapped diagonal for this matrix.
        matrix = np.array([[-1, 2], [2, -1]], dtype=np.float32)
        seq = np.array([0, 1], dtype=np.intp)
        workspace = ArrayWorkspace()
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
        self.assertNotIn("path", workspace.arrays)


class MultiDistanceTests(unittest.TestCase):
    def test_hand_calculation(self):
        # AC vs AG: S=0, S_max=2, S_rand=(1-1-1-1)/2=-1, d=ln(3).
        matrix = np.full((3, 3), -1.0, dtype=np.float64)
        np.fill_diagonal(matrix, 1.0)
        encoded = [np.array([0, 1]), np.array([0, 2]), np.array([0, 1])]
        dm = _multi_distances(
            encoded, matrix, 0.0, 2.0, False, list("abc"), ArrayWorkspace(), 0.0
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
                encoded, matrix, 1.0, 1.0, free, ["0", "1"], ArrayWorkspace(), 0.0
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
                encoded, matrix, *costs, free, ids, ArrayWorkspace(), dtype(1e-5)
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
