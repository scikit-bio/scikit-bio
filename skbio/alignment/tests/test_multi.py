# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import product
import unittest

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

from skbio import DNA, RNA, Protein, Sequence, TabularMSA, TreeNode, SubstitutionMatrix
from skbio.io import read as sk_read
from skbio.tree import TreeNode
from skbio.stats.distance import DistanceMatrix
from skbio.util._array import ArrayWorkspace
from skbio.alignment import AlignPath, multi_align, pair_align, align_score
from skbio.alignment._multi import MultiAlignResult, _merge_align, _score_dists


class MultiAlignTests(unittest.TestCase):
    def test_multi_align_nucl(self):
        """A simple case of three nucleotide sequences."""
        seqs = list(map(DNA, [
            "CAGCTATATATCGCTACG",
            "CTGCTTATATCCCTAGG",
            "AAGCTATACATCCAACATG",
        ]))

        # Default parameters (match = 1, mismatch = -1, gap = -2, free ends)
        res = multi_align(seqs)
        self.assertIs(type(res), MultiAlignResult)
        path = res[0]
        self.assertIs(type(path), AlignPath)
        self.assertIsNone(res[1])
        self.assertIsNone(res[2])
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

        # Return tuple with guide tree and distance matrix
        path, tree, dm = multi_align(seqs, keep_tree=True, keep_distmat=True)
        self.assertIs(type(path), AlignPath)
        self.assertIs(type(tree), TreeNode)
        self.assertIs(type(dm), DistanceMatrix)

        # Examine bits, guide tree and distance matrix
        exp = np.array(
            [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
             [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]],
            dtype=np.uint8)
        npt.assert_array_equal(path.to_bits(), exp)
        exp = TreeNode.read(["((0,1),2);"])
        self.assertEqual(tree.compare_rfd(exp), 0.0)
        self.assertTupleEqual(dm.ids, tuple("012"))
        exp = np.array([[0.,      0.36617, 0.51669],
                        [0.36617, 0.,      1.13943],
                        [0.51669, 1.13943, 0.     ]])
        npt.assert_array_equal(dm.data.round(5), exp)

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
        path = multi_align(seqs, **params).path
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
        path = multi_align(seqs, **params).path
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
        path = multi_align(seqs, free_ends=False).path
        self.assertIs(type(path), AlignPath)
        self.assertEqual(path.to_aligned(seqs), ["ACGT", "A-GT", "ACGT"])
        self.assertEqual(align_score((path, seqs), free_ends=False), 6)
        for free, gap in product([False, True], [2, (5, 2), (5, 0)]):
            path = multi_align(seqs, gap_cost=gap, free_ends=free).path
            self.assertEqual([s.replace("-", "") for s in path.to_aligned(seqs)], seqs)

    def test_multi_align_input(self):
        """Test various types of input sequences and parameters."""
        seqs = ["ACGT", "AGT", "ACGT"]

        # Alternative input sequence class; input as list; output as object
        for cls in [str, Sequence, DNA, RNA, Protein]:
            data = [cls(x.replace("T", "U") if cls is RNA else x) for x in seqs]
            path = multi_align(data, free_ends=False).path
            npt.assert_array_equal(
                path.to_bits(), [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]]
            )

        # Alternative input sequence type; input as iterator; output as tuple
        for data in [
            [[1, 2, 3], [1, 3], [1, 2, 3]],  # integers
            [["cat", "dog", "bird"], ["cat", "bird"], ["cat", "dog", "bird"]],  # words
            ["äëï", "äï", "äëï"],  # unicode
        ]:
            path, *_ = multi_align(iter(data), free_ends=False)
            npt.assert_array_equal(path.to_bits()[1], [0, 1, 0])

        # Custom substitution matrix
        submat = SubstitutionMatrix("ACGT", [
            [2, -1, -1, -1],
            [-1, 2, -1, -1],
            [-1, -1, 2, -1],
            [-1, -1, -1, 2],
        ])
        path, *_ = multi_align(seqs, sub_score=submat, free_ends=False)
        self.assertEqual(path.to_aligned(seqs)[1], "A-GT")

    def test_multi_align_error(self):
        """Test invalid input sequences and parameters."""

        for seqs in (["ACGT"], []):
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
            with self.assertRaisesRegex(ValueError, "finite and non-negative"):
                multi_align(["A", "A"], gap_cost=gap)

        for score in [
            (np.inf, -1), (1, np.nan), SubstitutionMatrix("AC", [[1, 0], [-1, 1]]),
        ]:
            with self.assertRaisesRegex(ValueError, "finite and symmetric"):
                multi_align(["AC", "AC"], sub_score=score)

    def test_multi_align_two(self):
        """Test only two sequences (will fallback to pairwise alignment)."""
        seq1 = DNA("GAATTC")
        seq2 = DNA("AGATCT")
        obs, *_ = multi_align((seq1, seq2))

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
        obs, *_ = multi_align((seq1, seq2), **params)
        exp = pair_align(seq1, seq2, **params).paths[0]
        npt.assert_array_equal(obs.to_bits(), exp.to_bits())

    def test_multi_align_tree(self):
        """Test guide tree's impact on alignment."""
        seqs = [DNA("CATTAACGT"),
                DNA("CGTTACGGT"),
                DNA("AGTTAACGG")]

        # Auto-determine merge order
        obs = multi_align(seqs).path.to_aligned(seqs)
        exp = ["CA-TTAACGT-",
               "-CGTTA-CGGT",
               "-AGTTAACGG-"]
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score(obs), 7.0)

        # Merging seqs 0 and 2 first produces the same output.
        tree = TreeNode.read(["((0,2),1);"])
        obs = multi_align(seqs, guide_tree=tree).path.to_aligned(seqs)
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score(obs), 7.0)

        # Merge seqs 0 and 1 first produces a less-aligned alignment, likely because
        # they are less similar, although the alignment score is the same.
        tree = TreeNode.read(["((0,1),2);"])
        obs = multi_align(seqs, guide_tree=tree).path.to_aligned(seqs)
        exp = ["CATTAACGT------",
               "------CGTTACGGT",
               "AGTTAACGG------"]
        self.assertListEqual(obs, exp)
        self.assertEqual(align_score(obs), 7.0)

        # Merge seqs 1 and 2 first produces an alternative and even better alignment.
        tree = TreeNode.read(["((1,2),0);"])
        obs = multi_align(seqs, guide_tree=tree).path.to_aligned(seqs)
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
        obs = multi_align(seqs, guide_tree=tree).path.to_aligned(seqs)
        self.assertListEqual(obs, exp)

        # Retain guide tree in output
        res = multi_align(seqs, guide_tree=tree, keep_tree=True)
        self.assertIs(res.tree, tree)

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
        obs = multi_align(seqs, free_ends=False).path.to_aligned(seqs)
        exp = ["ACGT", "A-GT", "ACGT"]
        self.assertListEqual(obs, exp)

        # Sequence IDs match guide tree.
        tree = TreeNode.read(["((c,a),b);"])
        obs = multi_align(seqs, guide_tree=tree, free_ends=False).path.to_aligned(seqs)
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
        path = multi_align(seqs, ids=iter(ids), guide_tree=tree, free_ends=False).path
        self.assertEqual([x.replace("-", "") for x in path.to_aligned(seqs)], seqs)
        self.assertEqual(str(tree), before)
        self.assertEqual(ids, ["a", "b", "c", "d"])
        npt.assert_array_equal(path.starts, [0, 0, 0, 0])


class MergeAlignTests(unittest.TestCase):
    """Tests of alignment merging operations."""

    def test_merge_align_walk(self):
        # This test uses a 5-sequence case to demonstrate the progressive alignment
        # procedures step-by-step.
        seqs = ["AGGCATG",
                "GCCATG",
                "AGCCAGGC",
                "TGACAATC",
                "TGAAATGC"]

        # Substitution matrix (A, C, G, T)
        submat = np.array([[ 1, -1, -1, -1],
                           [-1,  1, -1, -1],
                           [-1, -1,  1, -1],
                           [-1, -1, -1,  1]], dtype=np.float32)
        
        # Sequences encoded as indices in the substitution matrix
        encoded = [np.array([0, 2, 2, 1, 0, 3, 2]),
                   np.array([2, 1, 1, 0, 3, 2]),
                   np.array([0, 2, 1, 1, 0, 2, 2, 1]),
                   np.array([3, 2, 0, 1, 0, 0, 3, 1]),
                   np.array([3, 2, 0, 0, 0, 3, 2, 1])]

        # Linkage matrix defining merging order
        lnkmat = np.array([[3, 4],   # merge sequences 3 and 4
                           [0, 1],   # merge sequences 0 and 1
                           [2, 6],   # merge sequence 2 and alignment (0, 1)
                           [5, 7]])  # merge alignments (3, 4) and (0, 1, 2)

        # Re-usable memory space
        works = ArrayWorkspace()

        # Alignment parameters (the default)
        params = dict(
            submat=submat,
            gap_o=np.float32(0.0),
            gap_e=np.float32(2.0),
            free=False,
            works=works,
            atol=np.float32(0.0),
        )

        # For one-hot encoding
        eye = np.eye(len(submat))

        # A helper to convert a sequence into an alignment profile
        def _seq_profile(i):
            return (eye[enc := encoded[i]], np.zeros((1, len(enc)), dtype=bool), [i])

        # Test on one sequence
        aln3 = _seq_profile(3)
        exp_cnts = np.array([[0., 0., 0., 1.],
                             [0., 0., 1., 0.],
                             [1., 0., 0., 0.],
                             [0., 1., 0., 0.],
                             [1., 0., 0., 0.],
                             [1., 0., 0., 0.],
                             [0., 0., 0., 1.],
                             [0., 1., 0., 0.]], dtype=np.float32)
        npt.assert_array_equal(aln3[0], exp_cnts)
        exp_bits = np.array([[0, 0, 0, 0, 0, 0, 0, 0]], dtype=bool)
        npt.assert_array_equal(aln3[1], exp_bits)
        exp_order = [3]
        self.assertListEqual(aln3[2], exp_order)

        # A helper to construct aligned sequences
        def _get_aligned(bits, order):
            return AlignPath.from_bits(bits).to_aligned([seqs[i] for i in order])

        # Step 1: Merge sequences 3 and 4
        aln4 = _seq_profile(4)
        obs_cnts, obs_bits, obs_odr = aln34 = _merge_align(aln3, aln4, **params)

        exp_cnts = np.array([[0., 0., 0., 2.],
                             [0., 0., 2., 0.],
                             [2., 0., 0., 0.],
                             [0., 1., 0., 0.],
                             [2., 0., 0., 0.],
                             [2., 0., 0., 0.],
                             [0., 0., 0., 2.],
                             [0., 0., 1., 0.],
                             [0., 2., 0., 0.]], dtype=np.float32)
        self.assertEqual(obs_cnts.dtype, np.float32)
        npt.assert_array_equal(obs_cnts, exp_cnts)
        exp_bits = np.array([[0, 0, 0, 0, 0, 0, 0, 1, 0],
                             [0, 0, 0, 1, 0, 0, 0, 0, 0]], dtype=bool)
        self.assertEqual(obs_bits.dtype, np.bool_)
        npt.assert_array_equal(obs_bits, exp_bits)
        exp_odr = [3, 4]
        self.assertListEqual(obs_odr, exp_odr)
        obs_aln = _get_aligned(obs_bits, obs_odr)
        exp_aln = ["TGACAAT-C",
                   "TGA-AATGC"]
        self.assertListEqual(obs_aln, exp_aln)

        # Result should match pairwise alignment (including tie breaking)
        res = pair_align(seqs[3], seqs[4], free_ends=False, atol=0)
        exp_aln = res.paths[0].to_aligned((seqs[3], seqs[4]))
        self.assertListEqual(obs_aln, exp_aln)

        # Step 2: Merge sequences 0 and 1
        aln0 = _seq_profile(0)
        aln1 = _seq_profile(1)
        obs_cnts, obs_bits, obs_odr = aln01 = _merge_align(aln0, aln1, **params)
        obs_aln = _get_aligned(obs_bits, obs_odr)
        exp_cnts = np.array([[1., 0., 0., 0.],
                             [0., 0., 2., 0.],
                             [0., 1., 1., 0.],
                             [0., 2., 0., 0.],
                             [2., 0., 0., 0.],
                             [0., 0., 0., 2.],
                             [0., 0., 2., 0.]], dtype=np.float32)
        exp_bits = np.array([[0, 0, 0, 0, 0, 0, 0],
                             [1, 0, 0, 0, 0, 0, 0]], dtype=bool)
        exp_odr = [0, 1]
        exp_aln = ["AGGCATG",
                   "-GCCATG"]
        npt.assert_array_equal(obs_cnts, exp_cnts)
        npt.assert_array_equal(obs_bits, exp_bits)
        self.assertListEqual(obs_odr, exp_odr)
        self.assertListEqual(obs_aln, exp_aln)

        res = pair_align(seqs[0], seqs[1], free_ends=False, atol=0)
        exp_aln = res.paths[0].to_aligned((seqs[0], seqs[1]))
        self.assertListEqual(obs_aln, exp_aln)

        # Step 3: Merge sequence 2 and alignment (0, 1). This adds a gap to the end of
        # the alignment.
        aln2 = _seq_profile(2)
        obs_cnts, obs_bits, obs_odr = aln201 = _merge_align(aln2, aln01, **params)
        obs_aln = _get_aligned(obs_bits, obs_odr)
        exp_cnts = np.array([[2., 0., 0., 0.],
                             [0., 0., 3., 0.],
                             [0., 2., 1., 0.],
                             [0., 3., 0., 0.],
                             [3., 0., 0., 0.],
                             [0., 0., 1., 2.],
                             [0., 0., 3., 0.],
                             [0., 1., 0., 0.]], dtype=np.float32)
        exp_bits = np.array([[0, 0, 0, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 0, 0, 1],
                             [1, 0, 0, 0, 0, 0, 0, 1]], dtype=bool)
        exp_odr = [2, 0, 1]
        exp_aln = ["AGCCAGGC",
                   "AGGCATG-",
                   "-GCCATG-"]
        npt.assert_array_equal(obs_cnts, exp_cnts)
        npt.assert_array_equal(obs_bits, exp_bits)
        self.assertListEqual(obs_odr, exp_odr)
        self.assertListEqual(obs_aln, exp_aln)

        # Step 4: Merge alignments (3, 4) and (0, 1, 2). This introduces a gap in the
        # latter alignment.
        obs_cnts, obs_bits, obs_odr = _merge_align(aln34, aln201, **params)
        obs_aln = _get_aligned(obs_bits, obs_odr)
        exp_cnts = np.array([[2., 0., 0., 2.],
                             [0., 0., 5., 0.],
                             [2., 2., 1., 0.],
                             [0., 4., 0., 0.],
                             [5., 0., 0., 0.],
                             [2., 0., 0., 0.],
                             [0., 0., 1., 4.],
                             [0., 0., 4., 0.],
                             [0., 3., 0., 0.]], dtype=np.float32)
        exp_bits = np.array([[0, 0, 0, 0, 0, 0, 0, 1, 0],
                             [0, 0, 0, 1, 0, 0, 0, 0, 0],
                             [0, 0, 0, 0, 0, 1, 0, 0, 0],
                             [0, 0, 0, 0, 0, 1, 0, 0, 1],
                             [1, 0, 0, 0, 0, 1, 0, 0, 1]], dtype=bool)
        exp_odr = [3, 4, 2, 0, 1]
        exp_aln = ["TGACAAT-C",
                   "TGA-AATGC",
                   "AGCCA-GGC",
                   "AGGCA-TG-",
                   "-GCCA-TG-"]
        npt.assert_array_equal(obs_cnts, exp_cnts)
        npt.assert_array_equal(obs_bits, exp_bits)
        self.assertListEqual(obs_odr, exp_odr)
        self.assertListEqual(obs_aln, exp_aln)


class ScoreDistsTests(unittest.TestCase):

    def test_score_dists(self):
        seqs = ["ACGTACGT",
                "ACGTCGT",
                "ACGTACCT",
                "ACGTACG",
                "ACGTACGT"]
        encoded = [np.array(["ACGT".index(x) for x in seq], dtype=np.intp)
                   for seq in seqs]
        works = ArrayWorkspace()

        # normal case
        submat = np.full((4, 4), -1, dtype=np.float32)
        np.fill_diagonal(submat, 1)
        gap_o, gap_e = np.float32(0), np.float32(2)
        free = False
        atol = np.float32(0)
        
        obs = _score_dists(encoded, submat, gap_o, gap_e, free, works, atol)
        self.assertEqual(obs.dtype, np.float64)
        exp = np.array([[ 0.     ,  0.21357,  0.18232,  0.21357, -0.     ],
                        [ 0.21357,  0.     ,  0.42488,  0.43693,  0.21357],
                        [ 0.18232,  0.42488,  0.     ,  0.42488,  0.18232],
                        [ 0.21357,  0.43693,  0.42488,  0.     ,  0.21357],
                        [-0.     ,  0.21357,  0.18232,  0.21357,  0.     ]])
        npt.assert_array_equal(squareform(obs).round(5), exp)

        # Output is valid input for SciPy's linkage
        obs = linkage(obs, method="average")
        exp = np.array([[0, 4],
                        [2, 5],
                        [1, 6],
                        [3, 7]], dtype=np.intp)
        npt.assert_array_equal(obs[:, :2].astype(np.intp), exp)

        # Different parameters
        np.fill_diagonal(submat, 3)
        gap_o = np.float32(5)
        free = True
        obs = _score_dists(encoded, submat, gap_o, gap_e, free, works, atol)
        exp = np.array([[ 0.     ,  0.33987,  0.18232,  0.06899, -0.     ],
                        [ 0.33987,  0.     ,  0.55118,  0.43937,  0.33987],
                        [ 0.18232,  0.55118,  0.     ,  0.2803 ,  0.18232],
                        [ 0.06899,  0.43937,  0.2803 ,  0.     ,  0.06899],
                        [-0.     ,  0.33987,  0.18232,  0.06899,  0.     ]])
        npt.assert_array_equal(squareform(obs).round(5), exp)

        # Edge case: aligning AC, AG, AC
        # AC vs AC: d=0 (identical sequences)
        # AC vs AG: S=0, S_max=2, S_rand=(1-1-1-1)/2=-1, d=ln(3).
        submat = np.full((3, 3), -1.0, dtype=np.float64)
        np.fill_diagonal(submat, 1.0)
        encoded = [np.array([0, 1]), np.array([0, 2]), np.array([0, 1])]
        obs = _score_dists(encoded, submat, 0, 2, False, works, 0)
        npt.assert_allclose(obs, [np.log(3), 0, np.log(3)])

        # Edge case: AC vs A
        # Single terminal gap, L=2, expected substitution sum / L=0.
        submat = np.array([[1.0, -1.0], [-1.0, 1.0]])
        encoded = [np.array([0, 1]), np.array([0])]

        # Penalized gap: S=-1, S_rand=-2, S_max=1.5 -> d=-ln(1/3.5).
        obs = _score_dists(encoded, submat, 1, 1, False, works, 0)
        npt.assert_allclose(obs, [np.log(3.5)])

        # Free gap: S=1, S_rand=0, S_max=1.5 -> d=-ln(1/1.5).
        obs = _score_dists(encoded, submat, 1, 1, True, works, 0)
        npt.assert_allclose(obs, [np.log(1.5)])

        # Edge case: AAA vs AAA
        # S = S_max = S_rand -> d is undefined. The current implementation returns 0.
        encoded = [np.array([0, 0, 0]), np.array([0, 0, 0])]
        obs = _score_dists(encoded, submat, 1, 1, True, works, 0)
        npt.assert_allclose(obs, [0])

        # Edge case: match negative, mismatch positive
        # S=2, S_max=-3, S_rand=-4 -> S_eff=6. The current implementation clips at 1.
        submat = np.array([[-5.0, 1.0], [1.0, -5.0]])
        encoded = [np.array([0, 1]), np.array([1, 0])]
        obs = _score_dists(encoded, submat, 0, 2, False, works, 0)
        npt.assert_allclose(obs, [0])

        # Edge case: mismatch neutral
        # S=0, S_max=2, S_rand=1 -> S_eff=-1. The current implementation clips at 1e-6.
        submat = np.array([[1.0, 0.0], [0.0, 1.0]])
        obs = _score_dists(encoded, submat, 0, 2, False, works, 0)
        npt.assert_allclose(obs, [-np.log(1e-6)])


if __name__ == "__main__":
    unittest.main()
