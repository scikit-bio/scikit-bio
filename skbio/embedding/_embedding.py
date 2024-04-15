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
        self._embedding = embedding
        self._ids = ids

        if not isinstance(embedding, np.ndarray):
            raise ValueError("Input `embedding` must be a numpy array.")


class SequenceEmbedding(Embedding):
    """Store embeddings for a biological sequence."""

    def __init__(self, embedding, sequence, **kwargs):
        super(SequenceEmbedding, self).__init__(embedding, sequence, **kwargs)
