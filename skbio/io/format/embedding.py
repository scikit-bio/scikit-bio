"""Embedding format (:mod:`skbio.io.format.embed`)."""

import numpy as np
import h5py
from skbio.io import create_format, FASTAFormatError
from skbio.io.format.fasta import _sequence_to_fasta
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
    obj_kwargs: dict = {},
    embed_kwargs: dict = {},
):
    h5grp = h5py.File(fh, "r")
    embed_fh = h5grp["embedding"]
    id_fh = h5grp["id"]
    idptr_fh = h5grp["idptr"]

    n = embed_fh.shape[0]
    for i in range(n):
        emb = embed_fh[i]
        idptr = idptr_fh[i]
        id_ = id_fh[idptr : idptr_fh[i + 1]]
        string = str(id_.tobytes().decode("ascii"))
        yield embed_constructor(emb, string, **embed_kwargs)


def _embed_to_object(
    fh,
    embed_constructor=ProteinEmbedding,
    obj_constructor=Protein,
    obj_kwargs: dict = {},
    embed_kwargs: dict = {},
):
    h5grp = h5py.File(fh, "r")
    embed_fh = fh["embedding"]
    id_fh = fh["id"]

    # assumes that there is only a single object in the file
    emb = embed_fh[()]
    id_ = id_fh[()]

    string = str(id_.tobytes().decode("ascii"))
    return embed_constructor(emb, string, **embed_kwargs)


@embed.reader(ProteinEmbedding)
def _embed_to_protein(
    fh,
    obj_kwargs: dict = {},
    embed_kwargs: dict = {},
):
    return _embed_to_object(
        fh,
        embed_constructor=ProteinEmbedding,
        obj_constructor=Protein,
        obj_kwargs=obj_kwargs,
        embed_kwargs=embed_kwargs,
    )


def _object_to_embed(obj, fh, mode="w"):
    h5grp = h5py.File(fh, mode)

    # Store the embedding itself. We are assuming that the
    # embbedding is a 2D numpy array
    emb = obj.embedding
    emb = emb.reshape(1, emb.shape[0], emb.shape[1])

    if "embedding" in h5grp:
        embed_fh = h5grp["embedding"]
        n = emb.shape[0]
        embed_fh.resize(n, axis=0)
        embed_fh[-1] = emb

    else:
        embed_fh = h5grp.create_dataset(
            "embedding",
            data=emb,
            maxshape=(None, obj.embedding.shape[0], obj.embedding.shape[1]),
            dtype=obj.embedding.dtype,
        )

    # store string representation of the object
    # that will serve as an identifier for the entire object.
    # for sequences, this could be the sequence itself
    # for molecules, this could be the SMILES string.
    # The entries in this string representation can be used
    # to index the row vectors in the embedding.
    # For sequences, this is the positional index of the sequence.
    # For molecules, this is the position index of atoms in the SMILES string.
    arr = np.frombuffer(str(obj).encode("ascii"), dtype=np.uint8)
    if "id" in h5grp:
        id_fh = h5grp["id"]
        m = len(id)
        id_fh.resize(m + len(arr))
        id_fh[m : m + len(arr)] = arr
    else:
        id_fh = h5grp.create_dataset("id", data=arr, maxshape=(None,), dtype=np.int32)
    if "idptr" in h5grp:
        idptr = h5grp["idptr"]
        n = embed_eh.shape[0]
        idptr_fh = h5grp["ids"]
        idptr.resize(n)
        idptr_fh[n] = len(arr)
    else:
        idptr_fh = h5grp.create_dataset(
            "idptr", data=[len(arr)], maxshape=(None,), dtype=np.int32
        )


@embed.writer(None)
def _generator_to_embed(obj, fh):
    for it in obj:
        _object_to_embed(it, fh, mode="a")


@embed.writer(ProteinEmbedding)
def _protein_to_embed(obj, fh, mode="w"):
    _object_to_embed(it, fh, mode=mode)
