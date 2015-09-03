# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import collections

import numpy as np

from skbio._base import SkbioObject
from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.util import find_duplicates, OperationError, UniqueError
from skbio.util._decorator import experimental
from skbio.util._misc import resolve_key


_Shape = collections.namedtuple('Shape', ['sequence', 'position'])


class TabularMSA(SkbioObject):
    @property
    @experimental(as_of='0.4.0-dev')
    def dtype(self):
        return self._dtype

    @property
    @experimental(as_of='0.4.0-dev')
    def shape(self):
        return self._shape

    @property
    @experimental(as_of='0.4.0-dev')
    def keys(self):
        if self._keys is None:
            raise OperationError(
                "Keys do not exist. Use `reindex` to set them.")
        return self._keys

    @experimental(as_of='0.4.0-dev')
    def __init__(self, iterable, key=None):
        iterable = iter(iterable)
        seq = next(iterable, None)

        dtype = None
        length = 0
        seqs = []
        if seq is not None:
            seqs.append(seq)
            dtype = type(seq)
            if not issubclass(dtype, IUPACSequence):
                raise TypeError(
                    "`iterable` must contain scikit-bio sequence objects that "
                    "have an alphabet, not type %r" % dtype.__name__)
            length = len(seq)

            for seq in iterable:
                if type(seq) is not dtype:
                    raise TypeError(
                        "`iterable` cannot contain mixed types. Type %r does "
                        "not match type %r" %
                        (type(seq).__name__, dtype.__name__))
                if len(seq) != length:
                    raise ValueError(
                        "`iterable` must contain sequences of the same "
                        "length: %r != %r" % (len(seq), length))
                seqs.append(seq)

        self._seqs = seqs
        self._dtype = dtype
        self._shape = _Shape(sequence=len(seqs), position=length)
        self.reindex(key=key)

    @experimental(as_of='0.4.0-dev')
    def __bool__(self):
        # It is impossible to have 0 sequences and >0 positions.
        return self.shape.position > 0

    # Python 2 compatibility.
    __nonzero__ = __bool__

    @experimental(as_of='0.4.0-dev')
    def __len__(self):
        return self.shape.sequence

    @experimental(as_of='0.4.0-dev')
    def __iter__(self):
        return iter(self._seqs)

    @experimental(as_of='0.4.0-dev')
    def __reversed__(self):
        return reversed(self._seqs)

    @experimental(as_of='0.4.0-dev')
    def __str__(self):
        # TODO implement me!
        return super(TabularMSA, self).__str__()

    @experimental(as_of='0.4.0-dev')
    def __eq__(self, other):
        if not isinstance(other, TabularMSA):
            return False

        # Use np.array_equal instead of (a == b).all():
        #   http://stackoverflow.com/a/10580782/3776794
        return ((self._seqs == other._seqs) and
                np.array_equal(self._keys, other._keys))

    @experimental(as_of='0.4.0-dev')
    def __ne__(self, other):
        return not (self == other)

    @experimental(as_of='0.4.0-dev')
    def has_keys(self):
        return self._keys is not None

    @experimental(as_of='0.4.0-dev')
    def reindex(self, key=None):
        keys = None
        if key is not None:
            keys = [resolve_key(seq, key) for seq in self._seqs]
            duplicates = find_duplicates(keys)
            if duplicates:
                raise UniqueError(
                    "Keys must be unique. Duplicate keys: %r" % duplicates)
            keys = np.asarray(keys)
        self._keys = keys
