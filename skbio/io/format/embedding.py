r"""Embedding format (:mod:`skbio.io.format.embed`).
====================================================

.. currentmodule:: skbio.io.format.embed

This module provides support for reading and writing embedding files that
are outputted by sequential language models such as protein language models
(pLMs).

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |generator of :mod:`skbio.sequence.ProteinEmbedding` objects    |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.ProteinEmbedding` objects                 |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
The format is a HDF5 file with the following structure:

  - `embeddings` (dataset)
  - `id` (dataset)
  - `idptr` (dataset)

The `idptr` dataset contains the cumulative sum of the sequence lengths
in the hdf5. This is used to index both the sequences and the embeddings
in the hdf5, which can be useful for iterating through the embeddings and
avoiding the need to load all of the embedding into memory.  For protein
embeddings the `id` is the original sequence used to generate the embeddings.
The `embeddings` dataset contains the embeddings for each sequence, where the
first dimension is the sequence length and the second dimension is the
embedding dimension. The row vectors in the `embeddings` correspond to the
residues of the sequence in the `id` dataset.

Examples
--------
Here we will read in an example protein embedding file and write it back out.
Note that the embedding from implicitly gets the ``.write`` method from
the IO registry. This ``ByteIO`` object can be a file path in a regular
use case.

>>> import io, skbio
>>> f = io.BytesIO()
>>> skbio.embedding.example_protein_embedding.write(f)  # doctest: +ELLIPSIS
<_io.BytesIO object at ...>
>>> roundtrip = skbio.read(f, skbio.ProteinEmbedding)
>>> roundtrip
ProteinEmbedding
--------------------------------------------------------------------
Stats:
    length: 62
    embedding dimension: 1024
    has gaps: False
    has degenerates: False
    has definites: True
    has stops: False
--------------------------------------------------------------------
0  IGKEEIQQRL AQFVDHWKEL KQLAAARGQR LEESLEYQQF VANVEEEEAW INEKMTLVAS
60 ED
"""

import numpy as np
import h5py
from skbio.io import create_format
from skbio.embedding._protein import ProteinEmbedding
from skbio.sequence import Protein

embed = create_format("embed", encoding="binary")


@embed.sniffer()
def _embed_sniffer(fh):
    # this can be buffered, in which case .peek will return the buffer
    # so slice just in case
    magic = fh.peek(8)[:8]

    # From https://en.wikipedia.org/wiki/Hierarchical_Data_Format
    # Note that Wikipedia specifies: "\211HDF\r\n\032\n" which is an ordinal form:
    # >>> ord('\211')
    # 137
    # >>> ord('\x89')
    # 137
    # >>> ord('\032')
    # 26
    # >>> ord('\x1a')
    # 26
    if magic == b"\x89HDF\r\n\x1a\n":
        with h5py.File(fh, "r") as h5file:
            if "embedding" in h5file and "id" in h5file and "idptr" in h5file:
                return True, {}

    return False, {}


@embed.reader(None)
def _embed_to_generator(
    fh,
    embed_constructor=ProteinEmbedding,
    obj_constructor=Protein,
    kwargs: dict = {},
):
    h5grp = h5py.File(fh, "r")
    embed_fh = h5grp["embedding"]
    id_fh = h5grp["id"]
    idptr_fh = h5grp["idptr"]
    j = 0
    n = idptr_fh.shape[0]
    for i in range(n):
        idptr = idptr_fh[i]
        id_ = id_fh[j:idptr]
        emb = embed_fh[j:idptr]
        string = str(id_.tobytes().decode("ascii")).replace("\x00", "")
        j = idptr
        yield embed_constructor(emb, string, **kwargs)


def _embed_to_object(
    fh,
    embed_constructor=ProteinEmbedding,
    obj_constructor=Protein,
    kwargs: dict = {},
):
    h5grp = h5py.File(fh, "r")
    embed_fh = h5grp["embedding"]
    id_fh = h5grp["id"]

    # assumes that there is only a single object in the file
    emb = embed_fh[()].squeeze()
    id_ = id_fh[()]
    string = str(id_.tobytes().decode("ascii")).replace("\x00", "")
    return embed_constructor(emb, string, **kwargs)


@embed.reader(ProteinEmbedding)
def _embed_to_protein(
    fh,
    kwargs: dict = {},
):
    return _embed_to_object(
        fh,
        embed_constructor=ProteinEmbedding,
        obj_constructor=Protein,
        kwargs=kwargs,
    )


def _objects_to_embed(objs, fh):
    with h5py.File(fh, "w") as h5grp:
        for i, obj in enumerate(objs):
            # store string representation of the object
            # that will serve as an identifier for the entire object.
            # for sequences, this could be the sequence itself
            # for molecules, this could be the SMILES string.
            # The entries in this string representation can be used
            # to index the row vectors in the embedding.
            # For sequences, this is the positional index of the sequence.
            # For molecules, this is the position index of atoms in the SMILES string.
            arr = np.frombuffer(str(obj).encode("ascii"), dtype=np.uint8)
            # Store the embedding itself. We are assuming that the
            # embbedding is a 2D numpy array
            emb = obj.embedding
            # store the pointers that keep track of the start and
            # end of the embedding for each object, as well as well as
            # the corresponding string representation
            if "idptr" in h5grp:
                idptr_fh = h5grp["idptr"]
                idptr_fh.resize((i + 1,))
                idptr_fh[i] = len(arr) + idptr_fh[i - 1]
            else:
                idptr_fh = h5grp.create_dataset(
                    "idptr", data=[len(arr)], maxshape=(None,), dtype=np.int32
                )
            if "id" in h5grp:
                id_fh = h5grp["id"]
                id_fh.resize((idptr_fh[i],))
                id_fh[idptr_fh[i - 1] : idptr_fh[i]] = arr
            else:
                id_fh = h5grp.create_dataset(
                    "id", data=arr, maxshape=(None,), dtype=np.int32
                )

            if "embedding" in h5grp:
                embed_fh = h5grp["embedding"]
                assert embed_fh.shape[1] == emb.shape[1], (
                    "Embedding dimension mismatch between objects. "
                    f"({embed_fh.shape}) and ({emb.shape})"
                )

                embed_fh.resize(idptr_fh[i], axis=0)
                embed_fh[idptr_fh[i - 1] : idptr_fh[i]] = emb
            else:
                embed_fh = h5grp.create_dataset(
                    "embedding",
                    data=emb,
                    maxshape=(None, emb.shape[1]),
                    dtype=obj.embedding.dtype,
                )


@embed.writer(None)
def _generator_to_embed(objs, fh):
    return _objects_to_embed(objs, fh)


@embed.writer(ProteinEmbedding)
def _protein_to_embed(obj, fh):
    _objects_to_embed([obj], fh)
