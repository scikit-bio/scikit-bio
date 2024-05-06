# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------
from skbio.sequence import Protein
from skbio.embedding._embedding import SequenceEmbedding
from skbio.embedding._embedding import SequenceVector
from skbio.stats.ordination import OrdinationResults
from scipy.spatial.distance import pdist, squareform
from skbio.util import get_data_path
from pathlib import Path
import pandas as pd
import numpy as np
from typing import List


def _validate_protein(sequence):
    if isinstance(sequence, Protein):
        sequence = str(sequence)
    elif isinstance(sequence, str):
        if " " in sequence:
            sequence = sequence.replace(" ", "")

        # perform a check to make sure the sequence is a valid
        # protein sequence
        Protein(sequence)
    return sequence
    

class ProteinEmbedding(SequenceEmbedding):
    r"""Stores the embeddings of the protein sequence.

    Parameters
    ----------
    sequence : str, Protein, or 1D np.ndarray
        Characters representing the protein sequence itself.
    embedding : np.ndarray
        The embedding of the protein sequence. Row vectors correspond to
        the latent residues coordinates.
    clip_head : bool, optional
        If ``True``, then the first row of the embedding will be removed.
        Some language models specify start tokens, and this parameter can
        be used to account for this.
    clip_end : bool, optional
        If ``True``, then the last row of the embedding will be removed.
        Some language models specify end tokens, and this parameter can
        be used to account for this.

    See Also
    --------
    Protein

    """

    default_write_format = "embed"

    def __init__(
        self, embedding, sequence, clip_head=False, clip_tail=False, **kwargs
    ):
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
        r""" Array containing underlying residue characters """
        return self._ids.view("|S1")

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
        seq = Protein(str(self))

        rstr = repr(seq)
        rstr = rstr.replace("Protein", "ProteinEmbedding")
        n_indent = 4  # see Sequence.__repr__
        indent = " " * n_indent
        rstr = rstr.replace(
            "has gaps",
            f"embedding dimension: {self.embedding.shape[1]}\n{indent}has gaps",
        )
        return rstr


example_protein_embedding = ProteinEmbedding(
    np.random.randn(62, 1024),
    'IGKEEIQQRLAQFVDHWKELKQLAAARGQRLEESLEYQQFVANVEEEEAWINEKMTLVASED')


class ProteinVector(SequenceVector):
    """ A vector representation of the protein sequence.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray
        Characters representing the protein sequence itself.
    vector : np.ndarray
        The vector representation of the protein sequence.

    See Also
    --------
    Protein

    """
    default_write_format = "embed"

    def __init__(
        self, vector, sequence: str, **kwargs
    ):

        sequence = _validate_protein(sequence)
        if len(vector.shape) == 1:
            vector = vector.reshape(1, -1)
        if len(vector.shape) == 2:
            assert vector.shape[0] == 1, (
                "Only 1 vector per sequence is allowed."
            )
        super(ProteinVector, self).__init__(
            vector, sequence=sequence,  **kwargs
        )

    def __repr__(self):
        """
        Return a string representation of the ProteinVector object.

        Returns
        -------
        str
            A string representation of the ProteinEmbedding object.

        See Also
        --------
        Protein
        """
        seq = Protein(str(self))

        rstr = repr(seq)
        rstr = rstr.replace("Protein", "ProteinVector")
        n_indent = 4  # see Sequence.__repr__
        indent = " " * n_indent
        rstr = rstr.replace(
            "has gaps",
            f"vector dimension: {self.embedding.shape[1]}\n{indent}has gaps",
        )
        return rstr
