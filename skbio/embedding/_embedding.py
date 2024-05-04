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
    """Store embeddings for a biological object."""

    @property
    def embedding(self):
        return self._embedding

    @property
    def ids(self):
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

        self._embedding = np.array(embedding)
        self._ids = ids

    def __str__(self):
        return str(self._ids)


class SequenceEmbedding(Embedding):
    """Store embeddings for a biological sequence."""

    def __init__(self, embedding, sequence, **kwargs):
        super(SequenceEmbedding, self).__init__(embedding, sequence, **kwargs)

    @property
    def sequence(self):
        return str(self._ids)

    def __repr__(self):
        """
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
        rstr = rstr.replace(
            "has gaps",
            f"embedding dimension: {self.embedding.shape[1]}\n{indent}has gaps",
        )
        return rstr
