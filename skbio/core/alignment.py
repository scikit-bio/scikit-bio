#!/usr/bin/env python
r"""
Sequence collections and alignments (:mod:`skbio.core.alignment`)
================================================================

.. currentmodule:: skbio.core.alignment

This module provides functionality for working with biological sequence
collections and alignments. These can be composed of generic sequences, 
nucelotide sequences, DNA sequences, and RNA sequences. By default, input is
not validated, but all contructor methods take a validate option which checks
different features of the input based on ``SequenceCollection`` type.

Classes
-------

.. autosummary::
   :toctree: generated/

   SequenceCollection
   Alignment

Examples
--------
>>> from StringIO import StringIO
>>> from skbio.core.alignment import SequenceCollection, Alignment
>>> from skbio.core.sequence import DNA

>>> seqs = [DNA("ACC--G-GGTA.."), DNA("TCC--G-GGCA..")]
>>> a1 = Alignment(seqs)
>>> a1
<Alignment: n=2; mean +/- std length=13.00 +/- 0.00>

>>> seqs = [DNA("ACCGGG"), DNA("TCCGGGCA")]
>>> a1 = SequenceCollection(seqs)
>>> a1
<SequenceCollection: n=2; mean +/- std length=7.00 +/- 1.00>

>>> from skbio.parse.sequences import fasta_parse
>>> fasta_f = StringIO('>seq1\n'
...                    'CGATGTCGATCGATCGATCGATCAG\n'
...                    '>seq2\n'
...                    'CATCGATCGATCGATGCATGCATGCATG\n')
>>> s1 = SequenceCollection.from_fasta_records(fasta_parse(fasta_f), DNA)
>>> s1
<SequenceCollection: n=2; mean +/- std length=26.50 +/- 1.50>

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import Counter, defaultdict

import numpy as np

from scipy.stats import entropy

from skbio.core.exception import SequenceCollectionError
from skbio.core.distance import SymmetricDistanceMatrix

class SequenceCollection(object):
    """Base class for biological sequence collections.
    """

    @classmethod
    def from_fasta_records(cls, fasta_records, seq_constructor,
            validate=False):
        r"""Initialize a `SequenceCollection` object

        Parameters
        ----------
        fasta_records : generator of tuples
            The records to load into a new `SequenceCollection` object. These
            should be tuples of ``(sequence_identifier, sequence)``.
        seq_constructor : skbio.core.sequence.BiologicalSequence
        validate: bool, optional
            If True, runs the `is_valid` method after construction and raises
            `SequenceCollectionError` if ``is_valid == False``.

        Results
        -------
        SequenceCollection (or a derived class)
            The new `SequenceCollection` object.

        Raises
        ------
        skbio.core.exception.SequenceCollectionError
            If ``validate == True`` and ``is_valid == False``.

        See Also
        --------
        skbio.core.sequence.BiologicalSequence
        skbio.core.sequence.NucelotideSequence
        skbio.core.sequence.DNASequence
        skbio.core.sequence.RNASequence
        SequenceCollection
        Alignment
        skbio.parse.sequences
        skbio.parse.sequences.fasta_parse

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.parse.sequences import fasta_parse
        >>> from StringIO import StringIO
        >>> from skbio.core.sequence import DNA
        >>> fasta_f = StringIO('>seq1\nACCGT\n>seq2\nAACCGGT\n')
        >>> s1 = SequenceCollection.from_fasta_records(fasta_parse(fasta_f), DNA)
        >>> s1
        <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>
        
        >>> records = [('seq1', 'ACCGT'), ('seq2', 'AACCGGT')]
        >>> s1 = SequenceCollection.from_fasta_records(records, DNA)
        >>> s1
        <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

        """
        data = []
        for seq_id, seq in fasta_records:
            try:
                identifier, description = seq_id.split(None, 1)
            except ValueError:
                identifier = seq_id.strip()
                description = None
            data.append(seq_constructor(seq, identifier=identifier,
              description=description))
 
        return cls(data, validate=validate)

    def __init__(self, seqs, validate=False):
        r"""Initialize a `SequenceCollection` object

        Parameters
        ----------
        seqs : list of `BiologicalSequence` objects
            The `BiologicalSequence` objects to load into a new 
            `SequenceCollection` object.
        validate: bool, optional
            If True, runs the `is_valid` method after construction and raises
            `SequenceCollectionError` if ``is_valid == False``.

        Results
        -------
        SequenceCollection (or a derived class)
            The new `SequenceCollection` object.

        Raises
        ------
        skbio.core.exception.SequenceCollectionError
            If ``validate == True`` and ``is_valid == False``.

        See Also
        --------
        skbio.core.sequence.BiologicalSequence
        skbio.core.sequence.NucelotideSequence
        skbio.core.sequence.DNASequence
        skbio.core.sequence.RNASequence
        SequenceCollection
        Alignment
        skbio.parse.sequences
        skbio.parse.sequences.fasta_parse

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('ACCGT', identifier="seq1"), 
        ...              DNA('AACCGGT', identifier="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> s1
        <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

        """
        self._data = seqs
        self._identifier_to_index = dict([(seq.identifier, i) 
                for i, seq in enumerate(self._data)])
        if validate and not self.is_valid():
            raise SequenceCollectionError(
                "Something is wrong, and it's your fault.")

    def __eq__(self, other):
        """
        """
        if self.__class__ != other.__class__:
            return False
        elif len(self) != len(other):
            return False
        else:
            for self_seq, other_seq in zip(self,other):
                if self_seq != other_seq:
                    return False
        return True
    
    def __getitem__(self, index):
        """
        """
        if isinstance(index, str):
            return self.get_seq(index)
        else:
            return self._data[index]

    def __iter__(self):
        """
        """
        return iter(self._data)

    def __len__(self):
        """
        """
        return len(self._data)

    def __ne__(self, other):
        """
        """
        return not self.__eq__(other)

    def __repr__(self):
        """The repr method.

        Returns
        -------
        str
            Returns a string representation of the object.

        Notes
        -----

        Examples
        --------

        """
        cn = self.__class__.__name__
        count, center, spread = self.count_center_spread() 
        return "<%s: n=%d; mean +/- std length=%.2f +/- %.2f>" % (cn, count,
                center, spread)

    def count_center_spread(self, center_f=np.mean, spread_f=np.std):
        """
        """
        sequence_lengths = self.sequence_lengths()
        return (len(sequence_lengths), center_f(sequence_lengths),
                spread_f(sequence_lengths))

    def degap(self):
        """
        """
        result = []
        for seq in self:
            result.append(seq.degap())
        return SequenceCollection(result)

    def get_seq(self, identifier):
        """
        """
        return self[self._identifier_to_index[identifier]]

    def identifiers(self):
        """
        """
        return self._identifier_to_index.keys()

    def int_map(self, prefix=""):
        """
        """
        int_keys = []
        int_map = []
        for i, seq in enumerate(self):
            k = ("%s%d" % (prefix, i))
            int_map.append((k, seq))
            int_keys.append((k, seq.identifier))
        return dict(int_map), dict(int_keys) 

    def is_valid(self):
        """
        """
        return self._validate_character_set()

    def items(self):
        """
        """
        for seq in self:
            yield seq.identifier, seq

    def lower(self):
        """
        """
        result = []
        for seq in self:
            result.append(seq.lower())
        return self.__class__(result)

    def sequence_count(self):
        """
        """
        return len(self._data)

    def sequence_lengths(self):
        """
        """
        return [len(seq) for seq in self]

    def to_fasta(self):
        """
        """
        return ''.join([seq.to_fasta() for seq in self._data])
   
    def toFasta(self):
        print "Deprecation warning!"
        return self.to_fasta()
    
    def to_phylip(self):
        """
        """
        raise NotImplementedError

    def toPhylip(self):
        """
        """
        print "Deprecation warning!"
        return self.to_phylip()

    def upper(self):
        """
        """
        result = []
        for seq in self:
            result.append(seq.upper())
        return self.__class__(result)

    def _validate_character_set(self):
        """
        """
        for seq in self:
            if not seq.is_valid():
                return False
        return True


class Alignment(SequenceCollection):
    """
    """

    def distances(self):
        """
        """
        sequence_count = self.sequence_count()
        dm = np.zeros((sequence_count, sequence_count))
        identifiers = []
        for i in range(sequence_count):
            self_i = self[i]
            identifiers.append(self_i.identifier)
            for j in range(i):
                dm[i, j] = dm[j, i] = self_i.distance(self[j])
        return SymmetricDistanceMatrix(dm, identifiers)

    def subalignment(self, seqs_to_keep=None, positions_to_keep=None,
            invert_seqs_to_keep=False, invert_positions_to_keep=False):
        """
        """
        if seqs_to_keep is None:
            if invert_seqs_to_keep:
                return self.__class__([])
            else:
                def keep_seq(i, identifier):
                    return True
        else:
            seqs_to_keep = set(seqs_to_keep)
            if invert_seqs_to_keep:
                def keep_seq(i, identifier):
                    return not (identifier in seqs_to_keep or i in seqs_to_keep)
            else:
                def keep_seq(i, identifier):
                    return (identifier in seqs_to_keep or i in seqs_to_keep)

        if positions_to_keep is None:
            if invert_positions_to_keep:
                return self.__class__([])
            else:
                def keep_position(pos):
                    return True
        else:
            positions_to_keep = set(positions_to_keep)
            if invert_positions_to_keep:
                def keep_position(pos):
                    return not pos in positions_to_keep
            else:
                def keep_position(pos):
                    return pos in positions_to_keep

        result = []
        for i, seq in enumerate(self):
            if keep_seq(i, seq.identifier):
                new_seq = []
                for i,c in enumerate(seq):
                    if keep_position(i):
                        new_seq.append(c)
                # this is bad, we are calling join too much. this
                # should be addressed in issue #60
                result.append(seq.__class__(''.join(new_seq),
                    identifier = seq.identifier,
                    description=seq.description))
            else:
                continue
        return self.__class__(result)

    def is_valid(self):
        """
        """
        return super(Alignment, self).is_valid() and self._validate_lengths()

    def iter_positions(self, constructor=None):
        """
        """
        if constructor is None:
            def constructor(s):
                return s
        for i in range(self.sequence_length()):
            position = [constructor(seq[i]) for seq in self]
            yield position

    def majority_consensus(self, constructor=None):
        """
        """
        if constructor is None:
            constructor = self[0].__class__
        result = []
        for c in self.position_counters():
            # Counter.most_common returns an ordered list of the 
            # n most common (sequence, count) items in Counter. Here
            # we set n=1, and take only the character, not the count.
            result.append(c.most_common(1)[0][0])
        result = ''.join(result)
        return constructor(result)

    def omit_gap_positions(self, maximum_gap_frequency=0.0):
        """
        """
        position_frequencies = self.position_frequencies()
        gap_alphabet = self[0].gap_alphabet()
        
        positions_to_keep = []
        for i, f in enumerate(position_frequencies):
            gap_frequency = sum([f[c] for c in gap_alphabet])
            if gap_frequency <= maximum_gap_frequency:
                positions_to_keep.append(i)
        return self.subalignment(positions_to_keep=positions_to_keep)

    def omit_gap_sequences(self, maximum_gap_frequency=0.0):
        """
        """
        sequence_frequencies = self.sequence_frequencies()
        gap_alphabet = self[0].gap_alphabet()
        seqs_to_keep = []
        for seq, f in zip(self, sequence_frequencies):
            gap_frequency = sum([f[c] for c in gap_alphabet])
            if gap_frequency <= maximum_gap_frequency:
                seqs_to_keep.append(seq.identifier)
        return self.subalignment(seqs_to_keep=seqs_to_keep)

    def position_counters(self):
        """
        """
        result = []
        for p in self.iter_positions(constructor=str):
            result.append(Counter(p))
        return result

    def position_frequencies(self):
        """
        """
        result = []
        count = 1 / self.sequence_count()
        for p in self.iter_positions(constructor=str):
            current_freqs = defaultdict(int)
            for c in p:
                current_freqs[c] += count
            result.append(current_freqs)
        return result

    def position_entropies(self, base=None):
        """
        Parameters
        ----------
        base: float (optional)
            log base for entropy calculation (default: e)

        References
        ----------
        A Mathematical Theory of Communication, CE Shannon
        The Bell System Technical Journal (1948).
        """
        result = []
        for f in self.position_frequencies():
            result.append(entropy(f.values(), base=base))
        return result

    def sequence_frequencies(self):
        """
        """
        result = []
        count = 1 / self.sequence_length()
        for s in self:
            current_freqs = defaultdict(int)
            for c in s:
                current_freqs[c] += count
            result.append(current_freqs)
        return result

    def sequence_length(self):
        """
        """
        return len(self._data[0])

    def _validate_lengths(self):
        """
        """
        seq1_length = len(self[0])
        for seq in self[1:]:
            if seq1_length != len(seq):
                return False
        return True

