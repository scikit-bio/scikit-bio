# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re

from collections import defaultdict
import itertools

import numpy as np

from skbio.util._decorator import classproperty
from skbio._base import SkbioObject
from skbio.sequence import Protein, InvalidCodonError, GeneticCodeInitError

class GeneticCode(object):

    @classmethod
    def from_id(cls, id):
        return StandardGeneticCode

    @property
    def values(self):
        return self._values.copy()

    @classproperty
    def reading_frames(cls):
        return [1, 2, 3, -1, -2, -3]

    def __init__(self, amino_acids, starts):
        amino_acids = Protein(amino_acids)
        starts = Protein(starts)
        self._starts = set(str(amino_acids[starts.values == b'M']))

        self._values = np.zeros(64, dtype=np.uint8)
        codons = (''.join(x) for x in itertools.product('UCAG', repeat=3))
        for key, value in zip(codons, amino_acids):
            data = np.fromstring(key, dtype=np.uint8)
            data[data == ord('U')] = 0
            data[data == ord('C')] = 1
            data[data == ord('A')] = 2
            data[data == ord('G')] = 3
            index = (data * np.asarray([16, 4, 1], dtype=np.uint8)).sum()
            self._values[index] = ord(value.values[0])


StandardGeneticCode = GeneticCode(
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "---M---------------M---------------M----------------------------"
)
