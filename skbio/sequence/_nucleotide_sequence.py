# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import with_metaclass

import re
from abc import ABCMeta, abstractproperty

from skbio.sequence import SequenceError
from skbio.util import classproperty
from ._iupac_sequence import IUPACSequence


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

    @abstractproperty
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

    def complement(self, reverse=False):
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
        # TODO rewrite method for optimized performance
        complement_map = self.complement_map
        seq_iterator = reversed(self) if reverse else self
        result = [complement_map[str(base)] for base in seq_iterator]

        quality = self.quality
        if self._has_quality() and reverse:
            quality = self.quality[::-1]

        return self._to(sequence=''.join(result), quality=quality)

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
        return self.complement(reverse=True)

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
        # avoid computing the reverse complement if possible
        if len(self) != len(other):
            return False
        else:
            return other.reverse_complement()._string == self._string

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
        >>> from skbio import DNA
        >>> s = DNA('G-AT.T')
        >>> list(s.find_features('purine_run'))
        [slice(0, 1, None), slice(2, 3, None)]
        >>> list(s.find_features('purine_run', 2))
        []
        >>> list(s.find_features('purine_run', 2, allow_gaps=True))
        [slice(0, 3, None)]
        >>> list(s.find_features('pyrimidine_run', 2, allow_gaps=True))
        [slice(3, 6, None)]

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

        for hits in self.slices_from_regex(pat):
            if allow_gaps:
                degapped = self[hits].degap()
                if len(degapped) >= min_length:
                    yield hits
            else:
                yield hits
