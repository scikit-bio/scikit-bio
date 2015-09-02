# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio._base import SkbioObject
from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.util import find_duplicates, OperationError, UniqueError
from skbio.util._misc import resolve_key


class TabularMSA(SkbioObject):
    @property
    def keys(self):
        if self._keys is None:
            raise OperationError("")
        return self._keys

    def __init__(self, iterable, key=None):
        iterable = iter(iterable)
        seq = next(iterable, None)

        self._seqs = []
        if seq is not None:
            self._seqs.append(seq)
            dtype = type(seq)
            if not issubclass(dtype, IUPACSequence):
                raise TypeError("")
            length = len(seq)

            for seq in itertable:
                if type(seq) is not dtype:
                    raise TypeError("")
                if len(seq) != length:
                    raise ValueError("")
                self._seqs.append(seq)

        self._keys = None
        if key is not None:
            keys = [resolve_key(seq, key) for seq in self._seqs]
            repeats = find_duplicates(keys)
            if repeats:
                raise UniqueError(repeats)
            self._keys = np.asarray(keys)

    def __str__(self):
        return super(TabularMSA, self).__str__()
