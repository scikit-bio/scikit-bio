# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
from scipy.linalg import svd

from skbio.sequence import Sequence
from skbio._base import SkbioObject
from skbio.stats.ordination import OrdinationResults
from skbio.diversity import beta_diversity


def _repr_helper(rstr, org_name, new_name, dim_name, regex_match, shape):
    rstr = rstr.replace(org_name, new_name)
    n_indent = 4  # see Sequence.__repr__
    indent = " " * n_indent
    rstr = rstr.replace(
        regex_match,
        f"{dim_name} dimension: {shape}\n{indent}has gaps",
    )
    return rstr


class Embedding(SkbioObject):
    r"""Embedding for a biological object.

    Parameters
    ----------
    embedding : array_like
        Embedding matrix where the first axis is indexed by `ids`.
    ids : array_like
        List of biological object IDs.
    """

    @property
    def embedding(self):
        r"""The embedding tensor."""
        return self._embedding

    @property
    def ids(self):
        r"""IDs corresponding to each row of the embedding."""
        # each embedding row corresponds to an id
        return self._ids

    def __init__(self, embedding, ids, **kwargs):
        # make sure that the embedding has the same length as the sequence
        ids_len = len(ids)
        if embedding.shape[0] != ids_len:
            raise ValueError(
                f"The embedding ({embedding.shape[0]}) must have the "
                f"same length as the ids ({ids_len})."
            )

        self._embedding = np.asarray(embedding)
        self._ids = np.asarray(ids)

    def __str__(self):
        raise NotImplementedError("This method should be implemented by subclasses.")

    def bytes(self):
        r"""Bytes representation of string encoding."""
        seq = np.frombuffer(str(self).encode("ascii"), dtype=np.uint8)
        return seq


class SequenceEmbedding(Embedding):
    r"""Embedding for a biological sequence."""

    def __init__(self, embedding, sequence, **kwargs):
        if isinstance(sequence, Sequence):
            sequence = str(sequence)
        if isinstance(sequence, str):
            sequence = sequence.encode("ascii")
        seq = np.frombuffer(sequence, dtype=np.uint8)
        super(SequenceEmbedding, self).__init__(embedding, seq, **kwargs)

    def __str__(self):
        r"""String representation of the underlying sequence."""
        return str(self._ids.tobytes().decode("ascii"))

    @property
    def sequence(self):
        r"""String representation of the underlying sequence."""
        return str(self)

    def __repr__(self):
        r"""Return a string representation of the SequenceEmbedding object.

        Returns
        -------
        str
            String representation of the SequenceEmbedding object.

        See Also
        --------
        skbio.sequence.Protein

        """
        seq = Sequence(self.sequence)
        rstr = _repr_helper(
            repr(seq),
            "Sequence",
            "SequenceEmbedding",
            "embedding",
            regex_match="length",
            shape=self.embedding.shape[1],
        )
        return rstr


class EmbeddingVector(Embedding):
    r"""Vector representation for a biological entity."""

    def __init__(self, vector, obj, **kwargs):
        super(EmbeddingVector, self).__init__(vector, obj, **kwargs)

    def __str__(self):
        return self._ids[0].decode("ascii")

    @property
    def vector(self):
        r"""Vector representation for the biological entity."""
        return self._embedding.squeeze()

    @property
    def embedding(self):
        r"""The embedding tensor."""
        return self._embedding.reshape(1, -1)


class SequenceVector(EmbeddingVector):
    r"""Vector representation for a biological sequence."""

    def __init__(self, vector, sequence, **kwargs):
        if isinstance(sequence, Sequence):
            sequence = str(sequence)
        if isinstance(sequence, str):
            sequence = sequence.encode("ascii")

        vector = np.atleast_2d(vector)
        if vector.shape[0] != 1:
            raise ValueError("Only one vector per sequence is allowed.")

        seq = np.array([sequence], dtype="O")
        super(SequenceVector, self).__init__(vector, seq, **kwargs)

    @property
    def sequence(self):
        r"""String representation of the underlying sequence."""
        return str(self)

    def __repr__(self):
        r"""Return a string representation of the SequenceVector object.

        Returns
        -------
        str
            A string representation of the SequenceVector object.

        See Also
        --------
        skbio.sequence.Protein

        """
        seq = Sequence(str(self))
        rstr = _repr_helper(
            repr(seq),
            "Sequence",
            "SequenceVector",
            "vector",
            regex_match="length",
            shape=self.embedding.shape[1],
        )
        return rstr


def embed_vec_to_numpy(vectors, validate=True):
    r"""Convert an iterable of EmbeddingVector objects to a NumPy array.

    Parameters
    ----------
    vectors : iterable of EmbeddingVector objects
        An iterable of EmbeddingVector objects, or objects that
        subclass EmbeddingVector.
    validate : bool, optional
        If ``True``, validate that all vectors have the same length
        and are valid types.

    Returns
    -------
    ndarray of shape (n_objects, n_features)
        A NumPy array of shape (n_sequences, n_features) where n_features
        corresponds to the dimensionality of the latent space.

    Raises
    ------
    ValueError
        If the vectors do not have the same length.

    """
    if validate:
        subcls = [issubclass(type(ev), EmbeddingVector) for ev in vectors]
        if not all(subcls):
            raise ValueError(
                "Input iterable contains objects that "
                "do not subclass EmbeddingVector."
            )

        types = [type(ev) for ev in vectors]
        if not all(t == types[0] for t in types):
            raise ValueError("All objects must be of the same type.")

        lens = [len(ev.vector) for ev in vectors]
        if not all(ln == lens[0] for ln in lens):
            raise ValueError("All vectors must have the same length.")
    data = np.vstack([ev.vector for ev in vectors])
    return data


def embed_vec_to_distances(vectors, metric="euclidean", validate=True):
    r"""Convert EmbeddingVector objects to a DistanceMatrix object.

    Parameters
    ----------
    vectors : iterable of EmbeddingVector objects
        An iterable of EmbeddingVector objects, or objects that
        subclass EmbeddingVector.
    metric : str or callable, optional
        The distance metric to use. Must be a valid metric for
        ``scipy.spatial.distance.pdist``.
    validate : bool, optional
        If ``True``, validate that all vectors have the same length
        and are valid types.

    Returns
    -------
    DistanceMatrix
        A DistanceMatrix object.

    See Also
    --------
    skbio.stats.distance.DistanceMatrix
    """
    data = embed_vec_to_numpy(vectors)
    ids = [str(ev) for ev in vectors]
    return beta_diversity(metric, data, ids)


def embed_vec_to_ordination(vectors, validate=True):
    r"""Convert EmbeddingVector objects to an Ordination object.

    A singular value decomposition (SVD) is performed on the data.

    Parameters
    ----------
    vectors : iterable of EmbeddingVector objects
        An iterable of EmbeddingVector objects, or objects that
        subclass EmbeddingVector.
    validate : bool, optional
        If ``True``, validate that all vectors have the same length
        and are valid types.

    Returns
    -------
    OrdinationResults
        An Ordination object.

    See Also
    --------
    skbio.stats.ordination.OrdinationResults
    """
    data = embed_vec_to_numpy(vectors)
    u, s, vh = svd(data, full_matrices=False)
    eigvals = s**2
    short_name = "SVD"
    long_name = "Singular Value Decomposition"
    # note that we are moving half of the singular values
    # in the eigvals to the samples and the other half to the features
    # this is to help with the interpretation of the ordination
    # if visualizing with biplots
    ordr = OrdinationResults(
        short_method_name=short_name,
        long_method_name=long_name,
        eigvals=eigvals,
        proportion_explained=eigvals / eigvals.sum(),
        samples=pd.DataFrame(u @ np.diag(s), index=[str(ev) for ev in vectors]),
        features=pd.DataFrame(vh.T, index=range(data.shape[1])),
    )
    return ordr


def embed_vec_to_dataframe(vectors, validate=True):
    r"""Convert a list of SequenceVector objects to a pandas DataFrame.

    Parameters
    ----------
    vectors : iterable of EmbeddingVector objects
        An iterable of EmbeddingVector objects, or objects that
        subclass EmbeddingVector.
    validate : bool, optional
        If ``True``, validate that all vectors have the same length
        and are valid types.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the embedding vectors as rows.

    See Also
    --------
    pd.DataFrame

    """
    data = embed_vec_to_numpy(vectors)
    return pd.DataFrame(data, index=[str(ev) for ev in vectors])
