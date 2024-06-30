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
|Yes   |Yes   |generator of :mod:`skbio.embedding.ProteinEmbedding` objects   |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.embedding.ProteinEmbedding` objects                |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |generator of :mod:`skbio.embedding.ProteinVector` objects      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.embedding.ProteinVector` objects                   |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
The format is a HDF5 file with the following structure:

  - ``embeddings`` (dataset)
  - ``embedding_ptr`` (dataset)
  - ``id`` (dataset)
  - ``idptr`` (dataset)
  - ``format`` (attribute)
  - ``format-version`` (attribute)
  - ``dtype`` (attribute)
  - ``dim`` (attribute)


The `idptr` dataset contains the cumulative sum of the sequence lengths
in the hdf5. This is used to index both the sequences and the embeddings
in the hdf5, which can be useful for iterating through the embeddings and
avoiding the need to load all of the embedding into memory.  For protein
embeddings the `id` is the original sequence used to generate the embeddings.
The `embeddings` dataset contains the embeddings for each sequence, where the
first dimension is the sequence length and the second dimension is the
embedding dimension. The row vectors in the `embeddings` correspond to the
residues of the sequence in the `id` dataset.  The `embptr` is an optional
dataset that is used in case the `id` length is different from the `embedding`
length. This could be do to string formatting, for instance for dealing with
protein vectors, the id is a full length sequence, not a single residue.
The `emdptr` is used to separately keep track of the individual embeddings
in these scenarios.

The format attribute is a string that specifies the format of the embedding.
If the ``format`` attribute is present and has the value of `embed`, then
the file is a valid embedding file. The `format-version` attribute is a string
that specifies the version of the format. The `dtype` attribute is a string
that specifies the data type of the embeddings. Currently supported dtypes
include `float32` or `float64`.  The `dim` attribute is an integer that
specifies the dimensionality of the embeddings. The `embed` format currently
does not support storing embeddings with different dimensionality in the
same file.

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
>>> roundtrip = skbio.read(f, into=skbio.ProteinEmbedding)
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
from math import ceil
import h5py
from skbio.io import create_format
from skbio.embedding._protein import ProteinEmbedding
from skbio.embedding._protein import ProteinVector
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
    constructor=ProteinEmbedding,
    obj_constructor=Protein,
    kwargs: dict = {},
):
    h5grp = h5py.File(fh, "r")
    embed_fh = h5grp["embedding"]
    id_fh = h5grp["id"]
    idptr_fh = h5grp["idptr"]
    has_embedding_ptr = "embedding_ptr" in h5grp
    if has_embedding_ptr:
        embptr_fh = h5grp["embedding_ptr"]
    j, k = 0, 0
    n = idptr_fh.shape[0]
    for i in range(n):
        idptr = idptr_fh[i]
        id_ = id_fh[j:idptr]
        if has_embedding_ptr:
            embptr = embptr_fh[i]
            emb = embed_fh[k:embptr]
        else:
            emb = embed_fh[j:idptr]
        string = str(id_.tobytes().decode("ascii"))
        j = idptr
        k = embptr if has_embedding_ptr else None
        yield constructor(emb, string, **kwargs)


def _embed_to_object(
    fh,
    constructor=ProteinEmbedding,
    obj_constructor=Protein,
    kwargs: dict = {},
):
    h5grp = h5py.File(fh, "r")
    embed_fh = h5grp["embedding"]
    id_fh = h5grp["id"]

    # assumes that there is only a single object in the file
    emb = np.array(embed_fh[:])
    id_ = id_fh[()]
    string = str(id_.tobytes().decode("ascii"))
    return constructor(emb, string, **kwargs)


@embed.reader(ProteinEmbedding)
def _embed_to_protein(
    fh,
    kwargs: dict = {},
):
    return _embed_to_object(
        fh,
        constructor=ProteinEmbedding,
        obj_constructor=Protein,
        kwargs=kwargs,
    )


@embed.reader(ProteinVector)
def _vector_to_protein(
    fh,
    kwargs: dict = {},
):
    return _embed_to_object(
        fh,
        constructor=ProteinVector,
        obj_constructor=Protein,
        kwargs=kwargs,
    )


def _objects_to_embed(objs, fh, include_embedding_pointer=True):
    with h5py.File(fh, "w") as h5grp:
        h5grp.attrs["format"] = "embedding"
        h5grp.attrs["format-version"] = "1.0"
        max_idsize = 1
        max_embsize = 1
        resize_by = 1.38
        resize = False
        for i, obj in enumerate(objs):
            # store string representation of the object
            # that will serve as an identifier for the entire object.
            # for sequences, this could be the sequence itself
            # for molecules, this could be the SMILES string.
            # The entries in this string representation can be used
            # to index the row vectors in the embedding.
            # For sequences, this is the positional index of the sequence.
            # For molecules, this is the position index of atoms in the SMILES string.
            arr = obj.bytes()
            # Store the embedding itself. We are assuming that the
            # embbedding is a 2D numpy array
            emb = obj.embedding
            dtype = emb.dtype
            if "dtype" not in h5grp.attrs:
                h5grp.attrs["dtype"] = dtype.name
            if "dim" not in h5grp.attrs:
                h5grp.attrs["dim"] = emb.shape[1]

            # resize if necessary
            if i > 0:
                if include_embedding_pointer:
                    if (len(arr) + idptr_fh[i - 1]) > max_idsize or (
                        emb.shape[0] + embptr_fh[i - 1]
                    ) > max_embsize:
                        max_idsize = ceil(len(arr) + idptr_fh[i - 1]) * resize_by
                        max_embsize = ceil(emb.shape[0] + embptr_fh[i - 1]) * resize_by
                        resize = True
                else:
                    if len(arr) + idptr_fh[i - 1] > max_idsize:
                        max_idsize = ceil(len(arr) + idptr_fh[i - 1] * resize_by)
                        resize = True

            # store the pointers that keep track of the start and
            # end of the embedding for each object, as well as well as
            # the corresponding string representation
            if "idptr" in h5grp:
                idptr_fh = h5grp["idptr"]
                if resize:
                    idptr_fh.resize((ceil(i * resize_by),))
                idptr_fh[i] = len(arr) + idptr_fh[i - 1]
            else:
                idptr_fh = h5grp.create_dataset(
                    "idptr",
                    data=[len(arr)],
                    maxshape=(None,),
                    dtype=np.int32,
                    compression="gzip",
                )
            if "id" in h5grp:
                id_fh = h5grp["id"]
                if resize:
                    id_fh.resize((max_idsize,))
                id_fh[idptr_fh[i - 1] : idptr_fh[i]] = arr
            else:
                id_fh = h5grp.create_dataset(
                    "id", data=arr, maxshape=(None,), dtype=np.uint8, compression="gzip"
                )

            if include_embedding_pointer:
                if "embedding_ptr" in h5grp:
                    embptr_fh = h5grp["embedding_ptr"]
                    if resize:
                        embptr_fh.resize((ceil(i * resize_by),))
                    embptr_fh[i] = emb.shape[0] + embptr_fh[i - 1]

                else:
                    embptr_fh = h5grp.create_dataset(
                        "embedding_ptr",
                        data=[emb.shape[0]],
                        maxshape=(None,),
                        dtype=np.int32,
                    )
            if "embedding" in h5grp:
                embed_fh = h5grp["embedding"]
                assert embed_fh.shape[1] == emb.shape[1], (
                    "Embedding dimension mismatch between objects. "
                    f"({embed_fh.shape}) and ({emb.shape})"
                )
                if resize:
                    embed_fh.resize(max_embsize, axis=0)

                if include_embedding_pointer:
                    embed_fh[embptr_fh[i - 1] : embptr_fh[i]] = emb
                else:
                    embed_fh[idptr_fh[i - 1] : idptr_fh[i]] = emb

            else:
                embed_fh = h5grp.create_dataset(
                    "embedding",
                    data=emb,
                    maxshape=(None, emb.shape[1]),
                    dtype=obj.embedding.dtype,
                    compression="gzip",
                )

            resize = False

        # resize the datasets to the actual number of objects
        max_idsize = idptr_fh[i]
        max_embsize = embptr_fh[i] if include_embedding_pointer else max_idsize

        id_fh.resize((max_idsize,))
        idptr_fh.resize((i,))
        embed_fh.resize(max_embsize, axis=0)
        if include_embedding_pointer:
            embptr_fh.resize((i,))


@embed.writer(None)
def _generator_to_embed(objs, fh):
    return _objects_to_embed(objs, fh)


@embed.writer(ProteinEmbedding)
def _protein_to_embed(obj, fh):
    _objects_to_embed([obj], fh, include_embedding_pointer=False)


@embed.writer(ProteinVector)
def _protein_to_vector(obj, fh):
    _objects_to_embed([obj], fh, include_embedding_pointer=True)
