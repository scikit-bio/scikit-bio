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
from skbio.stats.distance import DistanceMatrix
from skbio.stats.ordination import OrdinationResults
from skbio.metadata._mixin import MetadataMixin
from scipy.spatial.distance import pdist, squareform
from typing import List


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
        self._ids = np.array(ids)

    def __str__(self):
        raise NotImplemented

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
        seq = Sequence(str(self._ids))

        rstr = repr(seq)
        rstr = rstr.replace("Sequence", "SequenceEmbedding")
        n_indent = 4  # see Sequence.__repr__
        indent = " " * n_indent
        dim = self.embedding.shape[1]
        rstr = rstr.replace(
            "has gaps",
            f"embedding dimension: {dim}\n{indent}has gaps",
        )
        return rstr

    
class SequenceVector(Embedding):
    r"""Store a vector representation for a biological sequence."""
    def __init__(self, vector, sequence, **kwargs):

        if isinstance(sequence, Sequence):
            sequence = str(sequence)
        if isinstance(sequence, str):
            sequence = sequence.encode("ascii")

        if len(vector.shape) == 1:
            vector = vector.reshape(1, -1)
        if len(vector.shape) == 2:
            assert vector.shape[0] == 1, (
                "Only 1 vector per sequence is allowed."
            )
            
        seq = np.array([sequence], dtype='O')
        super(SequenceVector, self).__init__(vector, seq, **kwargs)    

    def __str__(self):
        return str(self._ids[0].decode('ascii'))
            
    @property
    def sequence(self):
        r""" String representation of the underlying sequence """
        return str(self)

    @property
    def vector(self):
        return self.embedding
    
    def __repr__(self):
        r"""
        Return a string representation of the SequenceEmbedding object.

        Returns
        -------
        str
            A string representation of the SequenceEmbedding object.

        See Also
        --------
        Protein
        """
        seq = Sequence(str(self))

        rstr = repr(seq)
        rstr = rstr.replace("Sequence", "SequenceVector")
        n_indent = 4  # see Sequence.__repr__
        indent = " " * n_indent
        dim = self.embedding.shape[1]
        rstr = rstr.replace(
            "has gaps",
            f"embedding dimension: {dim}\n{indent}has gaps",
        )
        return rstr    

    @property
    def embedding(self):
        r""" The embedding tensor. """
        return self._embedding.reshape(1, -1)
    
    @staticmethod
    def to_numpy(sequence_vectors : List["SequenceVector"]):
        lens = [len(pv.vector) for pv in sequence_vectors]
        if not all(l == lens[0] for l in lens):
            raise ValueError("All vectors must have the same length.")
        data = np.vstack([pv.vector for pv in sequence_vectors])
        return data

    @staticmethod
    def to_distance_matrix(sequence_vectors : List["SequenceVector"],
                           metric='euclidean'):
        """
        Convert a SequenceVector object to a DistanceMatrix object.

        Parameters
        ----------
        sequence_vectors : iterable of SequenceVector objects
            An iterable of SequenceVector objects.
        metric : str, optional
            The distance metric to use. Must be a valid metric for
            `scipy.spatial.distance.pdist`.

        Returns
        -------
        DistanceMatrix
            A DistanceMatrix object.

        See Also
        --------
        DistanceMatrix
        """
        data = SequenceVector.to_numpy(sequence_vectors)
        ids = [str(pv) for pv in sequence_vectors]
        dm = squareform(pdist(data, metric))
        return DistanceMatrix(dm, ids=ids)

    @staticmethod
    def to_ordination(sequence_vectors : List["SequenceVector"]):
        """
        Convert a list of SequenceVector objects to an Ordination object.

        Parameters
        ----------
        sequence_vectors : iterable of SequenceVector objects
            An iterable of SequenceVector objects.

        Returns
        -------
        OrdinationResults
            An Ordination object.

        See Also
        --------
        OrdinationResults
        """
        data = SequenceVector.to_numpy(sequence_vectors)
        u, s, v = np.linalg.svd(data)
        eigvals = s ** 2
        ordr = OrdinationResults(
            short_method_name = 'SequenceVectors',
            long_method_name = 'SequenceVectors',
            eigvals = eigvals,
            proportion_explained = eigvals / eigvals.sum(),
            samples=pd.DataFrame(
                u * s, index=[str(pv) for pv in sequence_vectors]),
            features=pd.DataFrame(v.T * s, index=range(data.shape[1])),
        )
        return ordr

    @staticmethod
    def to_dataframe(sequence_vectors : List["SequenceVector"]):
        """
        Convert a list of SequenceVector objects to a pandas DataFrame.

        Parameters
        ----------
        sequence_vectors : iterable of SequenceVector objects
            An iterable of SequenceVector objects.

        Returns
        -------
        pd.DataFrame
            A pandas DataFrame.

        See Also
        --------
        pd.DataFrame
        """
        data = SequenceVector.to_numpy(sequence_vectors)
        df = pd.DataFrame(data, index=[str(pv) for pv in sequence_vectors])
        return df
