# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------
from skbio.sequence import Protein
from skbio.embedding._embedding import SequenceEmbedding


class ProteinEmbedding(SequenceEmbedding):
    """Stores the embeddings of the protein sequence.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray
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

    def __init__(
        self, embedding, sequence: Protein, clip_head=False, clip_tail=False, **kwargs
    ):
        if clip_head:
            embedding = embedding[1:]
        if clip_tail:
            embedding = embedding[:-1]

        assert isinstance(sequence, Protein) or isinstance(sequence, str)
        if isinstance(sequence, str):
            sequence = Protein(sequence)

        # make sure that the embedding has the same length as the sequence
        sequence_len = len(sequence)
        if embedding.shape[0] != sequence_len:
            raise ValueError("The embedding must have the same length as the sequence.")

        super(ProteinEmbedding, self).__init__(
            embedding=embedding, sequence=sequence, **kwargs
        )

    def __str__(self):
        return str(self._ids)
