r"""Biological Embeddings (:mod:`skbio.embedding`)
=================================================

.. currentmodule:: skbio.embedding

This module provides support for storing embeddings for biological objects, such
as protein embeddings outputted from protein language models (pLMs).


Embedding types
---------------

.. autosummary::
    :toctree: generated/

    Embedding
    SequenceEmbedding
    ProteinEmbedding


Embedding vector types
----------------------

.. autosummary::
    :toctree: generated/

    SequenceVector
    ProteinVector


Embedding vector utilities
--------------------------

.. autosummary::
    :toctree: generated/

    embed_vec_to_numpy
    embed_vec_to_distances
    embed_vec_to_ordination
    embed_vec_to_dataframe

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._embedding import (
    Embedding,
    SequenceEmbedding,
    EmbeddingVector,
    SequenceVector,
    embed_vec_to_numpy,
    embed_vec_to_distances,
    embed_vec_to_ordination,
    embed_vec_to_dataframe,
)

from ._protein import ProteinEmbedding, ProteinVector, example_protein_embedding


__all__ = [
    "Embedding",
    "SequenceEmbedding",
    "EmbeddingVector",
    "SequenceVector",
    "embed_vec_to_numpy",
    "embed_vec_to_distances",
    "embed_vec_to_ordination",
    "embed_vec_to_dataframe",
    "ProteinEmbedding",
    "ProteinVector",
    "example_protein_embedding",
]
