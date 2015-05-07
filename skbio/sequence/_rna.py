# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util import classproperty, overrides
from ._nucleotide_sequence import NucleotideSequence
from ._iupac_sequence import IUPACSequence


class RNA(NucleotideSequence):
    """Base class for RNA sequences.

    An `RNA` is a `NucelotideSequence` that is restricted to only
    containing characters used in the IUPAC RNA lexicon.

    Notes
    -----
    All uppercase and lowercase IUPAC RNA characters are supported.

    """

    @classproperty
    @overrides(NucleotideSequence)
    def complement_map(cls):
        """Return the mapping of characters to their complements.

        Returns
        -------
        dict
            Mapping of characters to their complements.

        Notes
        -----
        The complement of a gap character is itself.

        """
        comp_map = {
            'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y',
            'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }

        comp_map.update({c: c for c in cls.gap_chars})
        return comp_map

    @classproperty
    @overrides(IUPACSequence)
    def nondegenerate_chars(cls):
        """Return the non-degenerate IUPAC RNA characters.

        Returns
        -------
        set
            Non-degenerate IUPAC RNA characters.

        """
        return set("ACGU")

    @classproperty
    @overrides(IUPACSequence)
    def degenerate_map(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate RNA character to the set of
            non-degenerate IUPAC RNA characters it represents.

        """
        return {
            "R": set("AG"), "Y": set("CU"), "M": set("AC"), "K": set("UG"),
            "W": set("AU"), "S": set("GC"), "B": set("CGU"), "D": set("AGU"),
            "H": set("ACU"), "V": set("ACG"), "N": set("ACGU")
        }
