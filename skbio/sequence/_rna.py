# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

from skbio.util import classproperty, overrides
from ._nucleotide_sequence import NucleotideSequence
from ._iupac_sequence import IUPACSequence
from ._protein import Protein
from ._genetic_code import GeneticCode


class RNA(NucleotideSequence):
    """Store RNA sequence data and optional associated metadata.

    Only characters in the IUPAC RNA character set [1]_ are supported.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
        Characters representing the RNA sequence itself.
    metadata : dict, optional
        Arbitrary metadata which applies to the entire sequence.
    positional_metadata : Pandas DataFrame consumable, optional
        Arbitrary per-character metadata. For example, quality data from
        sequencing reads. Must be able to be passed directly to the Pandas
        DataFrame constructor.
    validate : bool, optional
        If ``True``, validation will be performed to ensure that all sequence
        characters are in the IUPAC RNA character set. If ``False``, validation
        will not be performed. Turning off validation will improve runtime
        performance. If invalid characters are present, however, there is
        **no guarantee that operations performed on the resulting object will
        work or behave as expected.** Only turn off validation if you are
        certain that the sequence characters are valid. To store sequence data
        that is not IUPAC-compliant, use ``Sequence``.
    case_insenstive : bool, optional
        If ``True``, lowercase sequence characters will be converted to
        uppercase characters in order to be valid IUPAC RNA characters.

    Attributes
    ----------
    values
    metadata
    positional_metadata
    alphabet
    gap_chars
    nondegenerate_chars
    degenerate_chars
    degenerate_map
    complement_map

    See Also
    --------
    DNA

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    Examples
    --------
    >>> from skbio import RNA
    >>> RNA('ACCGAAU')
    RNA
    -----------------------------
    Stats:
        length: 7
        has gaps: False
        has degenerates: False
        has non-degenerates: True
    -----------------------------
    0 ACCGAAU

    Convert lowercase characters to uppercase:

    >>> RNA('AcCGaaU', lowercase=True)
    RNA
    -----------------------------
    Stats:
        length: 7
        has gaps: False
        has degenerates: False
        has non-degenerates: True
    -----------------------------
    0 ACCGAAU

    """

    @classproperty
    @overrides(NucleotideSequence)
    def complement_map(cls):
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
        return set("ACGU")

    @classproperty
    @overrides(IUPACSequence)
    def degenerate_map(cls):
        return {
            "R": set("AG"), "Y": set("CU"), "M": set("AC"), "K": set("UG"),
            "W": set("AU"), "S": set("GC"), "B": set("CGU"), "D": set("AGU"),
            "H": set("ACU"), "V": set("ACG"), "N": set("ACGU")
        }

    def translate(self, genetic_code=1, reading_frame=1):
        """
        Parameters
        ----------
        genetic_code : int, GeneticCode, optional
            The genetic code to use, will lookup the GeneticCode by int id or
            will use the GeneticCode provided.
        reading_frame : {1, 2, 3, -1, -2, -3}
            The reading frame to use. The number indicates the base to start
            translation on. If negative, will perform a reverse complement
            first.

        Returns
        -------
        Protein
            The translated protein.
        """
        if self.has_degenerates():
            raise NotImplementedError("scikit-bio does not currently support"
                                      " translation of degenerate sequences.")

        if reading_frame not in GeneticCode.reading_frames:
            raise ValueError("`reading_frame` must be one of %r not %r" %
                             (GeneticCode.reading_frames, reading_frame))

        if not isinstance(genetic_code, GeneticCode):
            genetic_code = GeneticCode.from_id(genetic_code)

        data = self
        offset = reading_frame - 1
        if reading_frame < 0:
            data = self.reverse_complement()
            offset = 1 - reading_frame

        data = data._bytes[offset:].copy()

        data = data[:data.size // 3 * 3].reshape((-1, 3))
        data[data == ord('U')] = 0
        data[data == ord('C')] = 1
        data[data == ord('A')] = 2
        data[data == ord('G')] = 3


        metadata = None
        if self.has_metadata():
            metadata = self.metadata

        index = (data * np.asarray([16, 4, 1], dtype=np.uint8)).sum(axis=1)
        return Protein(genetic_code.values[index], metadata=metadata)
