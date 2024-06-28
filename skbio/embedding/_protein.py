# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd

from skbio.sequence import Protein
from skbio.stats.ordination import OrdinationResults
from skbio.embedding._embedding import SequenceEmbedding, SequenceVector, _repr_helper


def _validate_protein(sequence):
    if isinstance(sequence, Protein):
        sequence = str(sequence)
    elif isinstance(sequence, str):
        if " " in sequence:
            sequence = sequence.replace(" ", "")

        # perform a check to make sure the sequence is a valid protein sequence
        _ = Protein(sequence)
    return sequence


class ProteinEmbedding(SequenceEmbedding):
    r"""Embedding of a protein sequence.

    Parameters
    ----------
    embedding : array_like
        The embedding of the protein sequence. Row vectors correspond to
        the latent residues coordinates.
    sequence : str, Protein, or 1D ndarray
        Characters representing the protein sequence itself.
    clip_head : bool, optional
        If ``True``, then the first row of the embedding will be removed.
        Some language models specify start tokens, and this parameter can
        be used to account for this.
    clip_tail : bool, optional
        If ``True``, then the last row of the embedding will be removed.
        Some language models specify end tokens, and this parameter can
        be used to account for this.

    See Also
    --------
    skbio.sequence.Protein

    Examples
    --------
    >>> from skbio.embedding import ProteinEmbedding
    >>> import numpy as np
    >>> embedding = np.random.rand(10, 3)
    >>> sequence = "ACDEFGHIKL"
    >>> ProteinEmbedding(embedding, sequence)
    ProteinEmbedding
    --------------------------
    Stats:
        length: 10
        embedding dimension: 3
        has gaps: False
        has degenerates: False
        has definites: True
        has stops: False
    --------------------------
    0 ACDEFGHIKL

    """

    default_write_format = "embed"

    def __init__(self, embedding, sequence, clip_head=False, clip_tail=False, **kwargs):
        embedding = np.asarray(embedding)
        if clip_head:
            embedding = embedding[1:]
        if clip_tail:
            embedding = embedding[:-1]

        sequence = _validate_protein(sequence)
        super(ProteinEmbedding, self).__init__(
            embedding=embedding, sequence=sequence, **kwargs
        )

    @property
    def residues(self):
        r"""Array containing underlying residue characters."""
        return self._ids.view("|S1")

    def __repr__(self):
        r"""Return a string representation of the ProteinEmbedding object.

        Returns
        -------
        str
            String representation of the ProteinEmbedding object.

        See Also
        --------
        skbio.sequence.Protein

        """
        seq = Protein(str(self))
        rstr = _repr_helper(
            repr(seq),
            "Protein",
            "ProteinEmbedding",
            "embedding",
            regex_match="has gaps",
            shape=self.embedding.shape[1],
        )
        return rstr


example_protein_embedding = ProteinEmbedding(
    np.random.default_rng(0).normal(size=(62, 1024)),
    "IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQQFVANVEEEEAWINEKMTLVASED",
)


class ProteinVector(SequenceVector):
    r"""Vector representation for a protein sequence.

    Parameters
    ----------
    vector : 1D or 2D array_like
        The vector representation of the protein sequence. Typically a 1D array. Can
        also be a 2D array with only one row.
    sequence : str, Sequence, or 1D ndarray
        Characters representing the protein sequence itself.

    See Also
    --------
    SequenceVector
    skbio.sequence.Protein

    Examples
    --------
    >>> from skbio.embedding import ProteinVector
    >>> import numpy as np
    >>> vector = np.random.rand(10)
    >>> sequence = "ACDEFGHIKL"
    >>> ProteinVector(vector, sequence)
    ProteinVector
    --------------------------
    Stats:
        length: 10
        vector dimension: 10
        has gaps: False
        has degenerates: False
        has definites: True
        has stops: False
    --------------------------
    0 ACDEFGHIKL

    """

    default_write_format = "embed"

    def __init__(self, vector, sequence: str, **kwargs):
        sequence = _validate_protein(sequence)
        super(ProteinVector, self).__init__(vector, sequence=sequence, **kwargs)

    def __repr__(self):
        r"""Return a string representation of the ProteinVector object.

        Returns
        -------
        str
            String representation of the ProteinEmbedding object.

        See Also
        --------
        skbio.sequence.Protein

        """
        seq = Protein(str(self))
        rstr = _repr_helper(
            repr(seq),
            "Protein",
            "ProteinVector",
            "vector",
            regex_match="has gaps",
            shape=self.embedding.shape[1],
        )
        return rstr
