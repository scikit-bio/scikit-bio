# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from skbio.sequence import Sequence
from skbio._base import SkbioObject
from skbio.stats.ordination import OrdinationResults
from skbio.metadata._mixin import MetadataMixin
from skbio.diversity import beta_diversity
from scipy.linalg import svd
from typing import List


def _repr_helper(rstr, org_name, new_name, dim_name,
                 regex_match, shape):
    rstr = rstr.replace(org_name, new_name)
    n_indent = 4  # see Sequence.__repr__
    indent = " " * n_indent
    rstr = rstr.replace(
        regex_match,
        f"{dim_name} dimension: {shape}\n{indent}has gaps",
    )
    return rstr


class Embedding(SkbioObject):
    r"""Store embeddings for a biological object."""

    @property
    def embedding(self):
        r""" The embedding tensor. """
        return self._embedding

    @property
    def ids(self):
        r""" IDs corresponding to each row of the embedding. """
        # each embedding row corresponds to an id
        return self._ids

    def __init__(self, embedding, ids, **kwargs):
        r"""
        Parameters
        ----------
        embedding : array_like
           Embedding matrix where the first axis is indexed by `ids`
        ids : array_like
           List of ids
        """

        # make sure that the embedding has the same length as the sequence
        ids_len = len(ids)
        if embedding.shape[0] != ids_len:
            raise ValueError(
                f"The embedding ({embedding.shape[0]}) must have the "
                f"same length as the ids ({ids_len})."
            )

        self._embedding = np.array(embedding)
        self._ids = np.asarray(ids)

    def __str__(self):
        raise NotImplementedError("This method should be implemented by subclasses.")

    def bytes(self):
        r""" Bytes representation of string encoding"""
        seq = np.frombuffer(str(self).encode("ascii"),
                            dtype=np.uint8)
        return seq


class SequenceEmbedding(Embedding):
    r"""Store embeddings for a biological sequence."""

    def __init__(self, embedding, sequence, **kwargs):

        if isinstance(sequence, Sequence):
            sequence = str(sequence)
        if isinstance(sequence, str):
            sequence = sequence.encode("ascii")
        seq = np.frombuffer(sequence, dtype=np.uint8)
        super(SequenceEmbedding, self).__init__(embedding, seq, **kwargs)

    def __str__(self):
        r""" String representation of the underlying sequence """
        return str(self._ids.tobytes().decode('ascii'))

    @property
    def sequence(self):
        r""" String representation of the underlying sequence """
        return str(self)

    def __repr__(self):
        r"""
        Return a string representation of the ProteinEmbedding object.

        Returns
        -------
        str
            A string representation of the ProteinEmbedding object.

        See Also
        --------
        Protein
        """
        seq = Sequence(self.sequence)
        rstr = _repr_helper(
            repr(seq), "Sequence", "SequenceEmbedding", "embedding",
            regex_match="length", shape=self.embedding.shape[1]
        )
        return rstr


class EmbeddingVector(Embedding):
    r"""Store a vector representation for a biological entity."""
    def __init__(self, vector, obj, **kwargs):
        super(EmbeddingVector, self).__init__(vector, obj, **kwargs)

    def __str__(self):
        return self._ids[0].decode('ascii')

    @property
    def vector(self):
        return self._embedding.squeeze()

    @property
    def embedding(self):
        r""" The embedding tensor. """
        return self._embedding.reshape(1, -1)


class SequenceVector(EmbeddingVector):
    r"""Store a vector representation for a biological sequence."""
    def __init__(self, vector, sequence, **kwargs):

        if isinstance(sequence, Sequence):
            sequence = str(sequence)
        if isinstance(sequence, str):
            sequence = sequence.encode("ascii")

        vector = np.atleast_2d(vector)
        if vector.shape[0] != 1:
            raise ValueError("Only 1 vector per sequence is allowed.")

        seq = np.array([sequence], dtype='O')
        super(SequenceVector, self).__init__(vector, seq, **kwargs)

    @property
    def sequence(self):
        r""" String representation of the underlying sequence """
        return str(self)

    def __repr__(self):
        r"""
        Return a string representation of the SequenceEmbedding object.

        Returns
        -------
        str
            A string representation of the SequenceEmbedding object.

        See Also
        --------
        skbio.sequence.Protein
        """
        seq = Sequence(str(self))
        rstr = _repr_helper(
            repr(seq), "Sequence", "SequenceVector", "vector",
            regex_match="length", shape=self.embedding.shape[1]
        )
        return rstr


def embed_vec_to_numpy(embedding_vectors, validate=True):
    r""" Convert an iterable of EmbeddingVector objects to a numpy array.

    Parameters
    ----------
    embedding_vectors : iterable of EmbeddingVector objects
        An iterable of EmbeddingVector objects, or objects that
        subclass EmbeddingVector.
    validate : bool, optional
        If ``True``, validate that all vectors have the same length
        and are valid types

    Returns
    -------
    np.ndarray
        A numpy array of shape (n_sequences, n_features) where n_features
        corresponds to the dimensionality of the latent space

    Raises
    ------
    ValueError
        If the vectors do not have the same length.

    Notes
    -----

    """
    if validate:
        subcls = [issubclass(type(ev), EmbeddingVector) for ev in embedding_vectors]
        if not all(subcls):
            raise ValueError("Input iterable contains objects that "
                             "do not subclass EmbeddingVector.")

        types = [type(ev) for ev in embedding_vectors]
        if not all(t == types[0] for t in types):
            raise ValueError("All objects must be of the same type.")

        lens = [len(ev.vector) for ev in embedding_vectors]
        if not all(ln == lens[0] for ln in lens):
            raise ValueError("All vectors must have the same length.")
    data = np.vstack([ev.vector for ev in embedding_vectors])
    return data


def embed_vec_to_distances(embedding_vectors, metric='euclidean',
                           validate=True):
    r""" Convert EmbeddingVector objects to a DistanceMatrix object.

    Parameters
    ----------
    embedding_vectors : iterable of EmbeddingVector objects
        An iterable of EmbeddingVector objects, or objects that
        subclass EmbeddingVector.
    metric : str or Callable, optional
        The distance metric to use. Must be a valid metric for
        `scipy.spatial.distance.pdist`.
    validate : bool, optional
        If ``True``, validate that all vectors have the same length
        and are valid types

    Returns
    -------
    DistanceMatrix
        A DistanceMatrix object.

    See Also
    --------
    DistanceMatrix
    """
    data = embed_vec_to_numpy(embedding_vectors)
    ids = [str(sv) for sv in embedding_vectors]
    return beta_diversity(metric, data, ids)


def embed_vec_to_ordination(embedding_vectors, validate=True):
    r""" Convert iterable of EmbeddingVector objects to an Ordination object.

    A singular value decomposition (SVD) is performed on the data.

    Parameters
    ----------
    embedding_vectors : iterable of EmbeddingVector objects
        An iterable of EmbeddingVector objects, or objects that
        subclass EmbeddingVector.
    validate : bool, optional
        If ``True``, validate that all vectors have the same length
        and are valid types

    Returns
    -------
    OrdinationResults
        An Ordination object.

    See Also
    --------
    OrdinationResults
    """
    data = embed_vec_to_numpy(embedding_vectors)
    u, s, vh = svd(data, full_matrices=False)
    eigvals = s ** 2
    short_name = "SVD"
    long_name = "Singular Value Decomposition"
    # note that we are moving half of the singular values
    # in the eigvals to the samples and the other half to the features
    # this is to help with the interpretation of the ordination
    # if visualizing with biplots
    ordr = OrdinationResults(
        short_method_name = short_name,
        long_method_name = long_name,
        eigvals = eigvals,
        proportion_explained = eigvals / eigvals.sum(),
        samples=pd.DataFrame(
            u @ np.diag(s), index=[str(sv) for sv in embedding_vectors]),
        features=pd.DataFrame(vh.T, index=range(data.shape[1])),
    )
    return ordr

def embed_vec_to_dataframe(embedding_vectors, validate=True):
    r""" Convert a list of SequenceVector objects to a pandas DataFrame.

    Parameters
    ----------
    embedding_vectors : iterable of EmbeddingVector objects
        An iterable of EmbeddingVector objects, or objects that
        subclass EmbeddingVector.
    validate : bool, optional
        If ``True``, validate that all vectors have the same length
        and are valid types

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the embedding vector as rows.

    See Also
    --------
    pd.DataFrame
    """
    data = embed_vec_to_numpy(embedding_vectors)
    df = pd.DataFrame(data, index=[str(sv) for sv in embedding_vectors])
    return df
