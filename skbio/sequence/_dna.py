# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range
from future.utils import viewitems, with_metaclass
from future.standard_library import hooks
from six import string_types

import re
import collections
import numbers
from abc import ABCMeta, abstractmethod
from itertools import product
from contextlib import contextmanager

import numpy as np
from scipy.spatial.distance import hamming

from skbio._base import SkbioObject
from skbio.sequence import SequenceError
from skbio.util import classproperty, overrides
from skbio.util._misc import reprnator
from ._nucleotide_sequence import NucleotideSequence
from ._iupac_sequence import IUPACSequence

with hooks():
    from itertools import zip_longest


class DNA(NucleotideSequence):
    """Base class for DNA sequences.

    A `DNA` is a `NucelotideSequence` that is restricted to only
    containing characters used in IUPAC DNA lexicon.

    See Also
    --------
    NucleotideSequence
    Sequence

    Notes
    -----
    All uppercase and lowercase IUPAC DNA characters are supported.

    """

    @classproperty
    @overrides(NucleotideSequence)
    def complement_map(cls):
        """Return the mapping of characters to their complements.

        The complement of a gap character is itself.

        Returns
        -------
        dict
            Mapping of characters to their complements.

        """
        comp_map = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y',
            'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }

        comp_map.update({c: c for c in cls.gap_chars})
        return comp_map

    @classproperty
    @overrides(IUPACSequence)
    def nondegenerate_chars(cls):
        """Return the non-degenerate IUPAC DNA characters.

        Returns
        -------
        set
            Non-degenerate IUPAC DNA characters.

        """
        return set("ACGT")

    @classproperty
    @overrides(IUPACSequence)
    def degenerate_map(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate DNA character to the set of
            non-degenerate IUPAC DNA characters it represents.

        """
        return {
            "R": set("AG"), "Y": set("CT"), "M": set("AC"), "K": set("TG"),
            "W": set("AT"), "S": set("GC"), "B": set("CGT"), "D": set("AGT"),
            "H": set("ACT"), "V": set("ACG"), "N": set("ACGT")
        }
