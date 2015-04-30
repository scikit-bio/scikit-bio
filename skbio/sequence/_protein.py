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
from ._iupac_sequence import IUPACSequence

with hooks():
    from itertools import zip_longest


class Protein(IUPACSequence):
    """Base class for protein sequences.

    A `Protein` is a `Sequence` containing only characters
    used in the IUPAC protein lexicon.

    See Also
    --------
    Sequence

    Notes
    -----
    All uppercase and lowercase IUPAC protein characters are supported.

    """

    @classproperty
    @overrides(IUPACSequence)
    def nondegenerate_chars(cls):
        """Return the non-degenerate IUPAC protein characters.

        Returns
        -------
        set
            Non-degenerate IUPAC protein characters.

        """
        return set("ACDEFGHIKLMNPQRSTVWY")

    @classproperty
    @overrides(IUPACSequence)
    def degenerate_map(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate protein character to the set of
            non-degenerate IUPAC protein characters it represents.

        """
        return {
            "B": set("DN"), "Z": set("EQ"),
            "X": set("ACDEFGHIKLMNPQRSTVWY")
        }
