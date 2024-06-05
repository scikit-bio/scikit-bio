# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd
from scipy.spatial.distance import euclidean

from skbio.sequence import DNA, Protein
from skbio import DistanceMatrix, OrdinationResults
from skbio.embedding._protein import ProteinVector

from skbio.embedding._embedding import (
    Embedding,
    SequenceEmbedding,
    SequenceVector,
    embed_vec_to_numpy,
    embed_vec_to_dataframe,
    embed_vec_to_distances,
    embed_vec_to_ordination
)


class EmbeddingTests(TestCase):
    def setUp(self):
        self.emb = np.random.randn(62, 10)
        self.seq = "IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQQFVANVEEEEAWINEKMTLVASED"

    def test_id(self):
        emb, s = self.emb, self.seq
        p_emb = Embedding(emb, list(s))
        npt.assert_array_equal(p_emb.ids, np.array(list(s)))

    def test_embedding(self):
        emb, s = self.emb, self.seq
        p_emb = Embedding(emb, s)
        self.assertTupleEqual(p_emb.embedding.shape, (62, 10))

    def test_str(self):
        with self.assertRaises(NotImplementedError):
            Embedding(self.emb, self.seq).__str__()

    def test_assert_length(self):
        msg = "The embedding (62) must have the same length as the ids (63)."
        with self.assertRaises(ValueError) as cm:
            Embedding(self.emb, self.seq + "A")
        self.assertEqual(str(cm.exception), msg)


class SequenceEmbeddingTests(TestCase):
    def setUp(self):
        self.emb = np.random.randn(62, 10)
        self.seq = "IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQQFVANVEEEEAWINEKMTLVASED"

    def test_repr(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)
        rstr = repr(p_emb)
        self.assertIn("SequenceEmbedding", rstr)
        self.assertIn("62", rstr)
        self.assertIn("10", rstr)
        self.assertIn("IGKEEIQQRL", rstr)

    def test_str(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)

        self.assertEqual(p_emb.__str__(), s)
        self.assertEqual(p_emb.sequence, s)
        self.assertEqual(str(p_emb.ids.tobytes().decode("ascii")), s)

    def test_bytes(self):
        emb, s = self.emb, self.seq
        p_emb = SequenceEmbedding(emb, s)
        res = p_emb.bytes()
        res_str = str(res.tobytes().decode("ascii"))
        self.assertEqual(res_str, s)

    def test_init(self):
        emb, s = self.emb, self.seq

        # sequence as string
        p_emb = SequenceEmbedding(emb, s)
        self.assertTupleEqual(p_emb.embedding.shape, (62, 10))

        # sequence as bytes
        p_emb = SequenceEmbedding(emb, s.encode("ascii"))
        self.assertTupleEqual(p_emb.embedding.shape, (62, 10))

        # sequence as skbio.Sequence
        p_emb = SequenceEmbedding(emb, Protein(s))
        self.assertTupleEqual(p_emb.embedding.shape, (62, 10))

    def test_assert_length(self):
        msg = "The embedding (62) must have the same length as the ids (63)."
        with self.assertRaises(ValueError) as cm:
            SequenceEmbedding(self.emb, self.seq + "A")
        self.assertEqual(str(cm.exception), msg)


class SequenceVectorTests(TestCase):

    def setUp(self):
        # Create some sample SequenceVector objects for testing
        self.vector1 = np.array([1, 2, 3])
        self.vector2 = np.array([4, 5, 6])
        self.vector3 = np.array([7, 8, 9])
        self.bad_vector = np.array([7, 8])
        self.seq_vectors = [
            SequenceVector(self.vector1, "ACGT"),
            SequenceVector(self.vector2, "GCTA"),
            SequenceVector(self.vector3, "TTAG")
        ]

    def test_init(self):
        vec = np.array([1, 2, 3])
        seq = "ACGT"

        # sequence as string
        obs = SequenceVector(vec, seq)
        npt.assert_array_equal(obs.vector, vec)
        npt.assert_array_equal(obs.embedding, vec.reshape(1, -1))
        npt.assert_array_equal(obs.ids, np.array([b"ACGT"]))

        # sequence as bytes
        obs = SequenceVector(vec, seq.encode("ascii"))
        npt.assert_array_equal(obs.vector, vec)

        # sequence as skbio.Sequence
        obs = SequenceVector(vec, DNA(seq))
        npt.assert_array_equal(obs.vector, vec)

        # input is a matrix, not a vector
        vec2d = np.vstack([vec, vec])
        msg = "Only one vector per sequence is allowed."
        with self.assertRaisesRegex(ValueError, msg):
            SequenceVector(vec2d, seq)

    def test_vector(self):
        # Test if the vector attribute is set correctly
        for i, vector in enumerate([self.vector1, self.vector2, self.vector3]):
            npt.assert_array_equal(self.seq_vectors[i].vector, vector)

    def test_sequence(self):
        # Test if the sequence attribute is set correctly
        for i, sequence in enumerate(["ACGT", "GCTA", "TTAG"]):
            self.assertEqual(self.seq_vectors[i].sequence, sequence)

    def test_repr(self):
        # Test if the __repr__ method returns the correct string
        for seq_vector in self.seq_vectors:
            self.assertTrue(seq_vector.__repr__().startswith("SequenceVector"))
            self.assertIn("vector", seq_vector.__repr__())

            # check latent dimension
            self.assertIn("4", seq_vector.__repr__())

    def test_str(self):
        # Test if the __str__ method returns the correct string
        for seq_vector in self.seq_vectors:
            self.assertEqual(str(seq_vector), seq_vector.sequence)


class EmbedVecUtilityTests(TestCase):

    def setUp(self):
        self.vector1 = np.array([1, 2, 3])
        self.vector2 = np.array([4, 5, 6])
        self.vector3 = np.array([7, 8, 9])
        self.bad_vector = np.array([7, 8])
        self.seq_vectors = [
            SequenceVector(self.vector1, "ACGT"),
            SequenceVector(self.vector2, "GCTA"),
            SequenceVector(self.vector3, "TTAG")
        ]

    def test_embed_vec_to_numpy(self):
        # Test if to_numpy returns the correct numpy array
        exp = np.array([self.vector1, self.vector2, self.vector3])
        obs = embed_vec_to_numpy(self.seq_vectors)
        npt.assert_array_equal(obs, exp)

        # skip validation
        obs = embed_vec_to_numpy(self.seq_vectors, validate=False)
        npt.assert_array_equal(obs, exp)

    def test_embed_vec_to_numpy_raises(self):
        # input contains non-vector
        lst = [SequenceVector(self.vector1, "ACGT"),
               SequenceEmbedding(np.vstack([self.vector2, self.vector3]), "AT")]
        msg = "Input iterable contains objects that do not subclass EmbeddingVector."
        with self.assertRaisesRegex(ValueError, msg):
            embed_vec_to_numpy(lst)

        # mixed sequence types
        lst = [SequenceVector(self.vector1, "ACGT"),
               ProteinVector(self.vector2, "MKRPL")]
        msg = "All objects must be of the same type."
        with self.assertRaisesRegex(ValueError, msg):
            embed_vec_to_numpy(lst)

        # lengths are not equal
        lst = [SequenceVector(self.vector1, "ACGT"),
               SequenceVector(self.vector2, "GCTA"),
               SequenceVector(self.bad_vector, "TTAG")]
        msg = "All vectors must have the same length."
        with self.assertRaisesRegex(ValueError, msg):
            embed_vec_to_numpy(lst)

    def test_embed_vec_to_distances(self):
        # Test if to_distances returns a DistanceMatrix object
        obs = embed_vec_to_distances(self.seq_vectors)
        self.assertIsInstance(obs, DistanceMatrix)
        self.assertTupleEqual(obs.shape, (3, 3))
        self.assertTrue(all(isinstance(d, float) for d in obs.condensed_form()))

        d12 = euclidean(self.vector1, self.vector2)
        d13 = euclidean(self.vector1, self.vector3)
        d23 = euclidean(self.vector2, self.vector3)
        exp = DistanceMatrix([[0, d12, d13],
                              [d12, 0, d23],
                              [d13, d23, 0]],
                             ids=["ACGT", "GCTA", "TTAG"])
        npt.assert_allclose(obs.data, exp.data)
        self.assertEqual(obs.ids, exp.ids)

        obs = embed_vec_to_distances(self.seq_vectors, validate=False)
        self.assertIsInstance(obs, DistanceMatrix)

    def test_embed_vec_to_ordination(self):
        # Test if to_ordination returns an OrdinationResults object
        obs = embed_vec_to_ordination(self.seq_vectors)
        self.assertIsInstance(obs, OrdinationResults)
        self.assertEqual(obs.samples.shape, (3, 3))
        self.assertEqual(obs.features.shape, (3, 3))
        reconstructed = (obs.samples.values @ obs.features.values.T)
        npt.assert_allclose(
            reconstructed, embed_vec_to_numpy(self.seq_vectors)
        )

        obs = embed_vec_to_ordination(self.seq_vectors, validate=False)
        self.assertIsInstance(obs, OrdinationResults)

    def test_embed_vec_to_dataframe(self):
        # Test if to_dataframe returns a pandas DataFrame object
        obs = embed_vec_to_dataframe(self.seq_vectors)
        self.assertIsInstance(obs, pd.DataFrame)
        self.assertTupleEqual(obs.shape, (3, 3))

        exp = pd.DataFrame([self.vector1, self.vector2, self.vector3],
                            index=["ACGT", "GCTA", "TTAG"])
        pd.testing.assert_frame_equal(obs, exp)

        obs = embed_vec_to_dataframe(self.seq_vectors, validate=False)
        self.assertIsInstance(obs, pd.DataFrame)


if __name__ == "__main__":
    main()
