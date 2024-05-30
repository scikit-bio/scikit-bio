r""" Biological Embeddings (:mod:`skbio.embedding`)
=================================================

.. currentmodule:: skbio.embedding

This module provides support for storing embeddings for biological objects, such
as protein embeddings outputted from protein language models (pLMs).

Embedding types
---------------

.. autosummary::
    :toctree: generated/

    ProteinEmbedding
    ProteinVector

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.embedding._protein import (
    ProteinEmbedding, ProteinVector, example_protein_embedding
)


__all__ = ["ProteinEmbedding", "ProteinVector", "example_protein_embedding"]
