# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.sequence import Sequence
from skbio._base import SkbioObject
from skbio.metadata._mixin import MetadataMixin


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
