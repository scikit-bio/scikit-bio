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


class NucleotideSequence(with_metaclass(ABCMeta, IUPACSequence)):
    """Base class for nucleotide sequences.

    A `NucleotideSequence` is a `Sequence` with additional methods
    that are only applicable for nucleotide sequences, and containing only
    characters used in the IUPAC DNA or RNA lexicon.

    See Also
    --------
    Sequence

    Notes
    -----
    All uppercase and lowercase IUPAC DNA/RNA characters are supported.

    """

    @abstractmethod
    @classproperty
    def complement_map(cls):
        """Return the mapping of characters to their complements.

        Returns
        -------
        dict
            Mapping of characters to their complements.

        Notes
        -----
        Complements cannot be defined for a generic `NucleotideSequence`
        because the complement of 'A' is ambiguous.
        `NucleotideSequence.complement_map` will therefore be the empty dict.
        Thanks, nature...

        """
        pass

    @abstractmethod
    @classproperty
    def nondegenerate_chars(cls):
        """Return the non-degenerate IUPAC nucleotide characters.

        Returns
        -------
        set
            Non-degenerate IUPAC nucleotide characters.

        """
        pass

    @abstractmethod
    @classproperty
    def degenerate_map(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate nucleotide character to the set of
            non-degenerate IUPAC nucleotide characters it represents.

        """
        pass

    def _complement(self, reverse=False):
        """Returns `NucleotideSequence` that is (reverse) complement of `self`.

        Parameters
        ----------
        reverse : bool, optional
            If ``True``, reverse `self` before complementing.

        Returns
        -------
        NucelotideSequence
            The (reverse) complement of `self`. Specific type will be the same
            as ``type(self)``.

        Raises
        ------
        skbio.sequence.SequenceError
            If a character is present in the `NucleotideSequence` that is not
            in the complement map.

        Notes
        -----
        This private method centralizes the logic for `complement` and
        `reverse_complement`.

        """
        result = []
        complement_map = self.complement_map
        seq_iterator = reversed(self) if reverse else self
        for base in seq_iterator:
            # TODO fix me!
            base = str(base)
            try:
                result.append(complement_map[base])
            except KeyError:
                raise SequenceError(
                    "Don't know how to complement base %s. Is it in "
                    "%s.complement_map?" % (base, self.__class__.__name__))

        quality = self.quality
        if self._has_quality() and reverse:
            quality = self.quality[::-1]

        return self.to(sequence=''.join(result), quality=quality)

    def complement(self):
        """Return the complement of the `NucleotideSequence`

        Returns
        -------
        NucelotideSequence
            The complement of `self`. Specific type will be the same as
            ``type(self)``.

        Raises
        ------
        skbio.sequence.SequenceError
            If a character is present in the `NucleotideSequence` that is not
            in `self.complement_map`.

        See Also
        --------
        reverse_complement
        complement_map

        Notes
        -----
        The type, id, description, and quality scores of the result will be the
        same as `self`.

        """
        return self._complement()

    def is_reverse_complement(self, other):
        """Return True if `other` is the reverse complement of `self`

        Returns
        -------
        bool
            `True` if `other` is the reverse complement of `self` and `False`
            otherwise.

        Raises
        ------
        skbio.sequence.SequenceError
            If a character is present in `other` that is not in the
            `self.complement_map`.

        See Also
        --------
        reverse_complement

        """
        return self.equals(other.reverse_complement(),
                           ignore=['id', 'description', 'quality'])

    def reverse_complement(self):
        """Return the reverse complement of the `NucleotideSequence`

        Returns
        -------
        NucelotideSequence
            The reverse complement of `self`. Specific type will be the same as
            ``type(self)``.

        Raises
        ------
        skbio.sequence.SequenceError
            If a character is present in the `NucleotideSequence` that is not
            in `self.complement_map`.

        See Also
        --------
        complement
        complement_map
        is_reverse_complement

        Notes
        -----
        The type, id, and description of the result will be the same as `self`.
        If quality scores are present, they will be reversed and included in
        the resulting biological sequence.

        """
        return self._complement(reverse=True)
    rc = reverse_complement

    def find_features(self, feature_type, min_length=1, allow_gaps=False):
        """Search the sequence for features

        Parameters
        ----------
        feature_type : {'purine_run', 'pyrimidine_run'}
            The type of feature to find
        min_length : int, optional
            Defaults to 1. Only features at least as long as this will be
            returned
        allow_gaps : bool, optional
            Defaults to ``False``. If ``True``, then gaps will not be
            considered to disrupt a feature

        Returns
        -------
        generator
            Yields tuples of the start of the feature, the end of the feature,
            and the subsequence that composes the feature

        Examples
        --------
        >>> from skbio.sequence import NucleotideSequence
        >>> s = NucleotideSequence('G-AT.T')
        >>> list(s.find_features('purine_run'))
        [(0, 1, 'G'), (2, 3, 'A')]
        >>> list(s.find_features('purine_run', 2))
        []
        >>> list(s.find_features('purine_run', 2, allow_gaps=True))
        [(0, 3, 'G-A')]
        >>> list(s.find_features('pyrimidine_run', 2, allow_gaps=True))
        [(3, 6, 'T.T')]

        """
        gaps = re.escape(''.join(self.gap_chars))
        acceptable = gaps if allow_gaps else ''

        if feature_type == 'purine_run':
            pat_str = '([AGag%s]{%d,})' % (acceptable, min_length)
        elif feature_type == 'pyrimidine_run':
            pat_str = '([CTUctu%s]{%d,})' % (acceptable, min_length)
        else:
            raise ValueError("Unknown feature type: %s" % feature_type)

        pat = re.compile(pat_str)

        for hits in self.regex_iter(pat):
            if allow_gaps:
                degapped = hits[2]
                for gap_char in self.gap_chars:
                    degapped = degapped.replace(gap_char, '')
                if len(degapped) >= min_length:
                    yield hits
            else:
                yield hits
