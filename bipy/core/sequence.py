#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division
from collections import Sequence
from itertools import izip

from bipy.core.exception import BiologicalSequenceError


class BiologicalSequence(Sequence):
    """ Base class for biological sequences 
        
        BiologicalSequence objects are immutable. Where applicable, 
         methods return a new object of the same class.
    """
    
    _alphabet = set()
    _gap_alphabet = set('-.')

    def __init__(self, sequence, identifier="", description="", 
                 validate=False):
        """ initialize a BiologicalSequence object

            sequence: the biological sequence as a python Sequence
             (e.g., a string, list, or tuple)
            identifier: the sequence identifier (e.g., an accession number)
            description: a description or comment about the sequence (e.g.,
             "green fluorescent protein")
            validate: if True, runs the is_valid method after construction
             and raises BiologicalSequenceError if is_valid == False
        """
        self._sequence = ''.join(sequence)
        self._identifier = identifier
        self._description = description

        if validate and not self.is_valid():
            unsupported_chars = self.unsupported_characters()
            raise BiologicalSequenceError(
                "Sequence contains unsupported characters: %s" 
                % (" ".join(unsupported_chars)))
 
    def __contains__(self, other):
        """ return True if other is contained in the BiologicalSequence
        """
        return other in self._sequence
   
    def __eq__(self, other):
        """ equality (==) operator
            
            BiologicalSequences are equal if their sequence is the same and
             they are the same type
        """
        if self.__class__ != other.__class__:
            return False
        elif self._sequence != other._sequence:
            return False
        else:
            return True

    def __getitem__(self, i):
        try:
            return self._sequence[i]
        except IndexError:
            raise IndexError(
                "Position %d is out of range for %r." % (i, self))
 
    def __hash__(self):
        return hash(self._sequence)
   
    def __iter__(self):
        return iter(self._sequence)

    def __len__(self):
        return len(self._sequence)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        first_ten = str(self)[:10]
        cn = self.__class__.__name__
        length = len(self)
        if length > 10:
            elipses = "..."
        else:
            elipses = ""
        return '<%s: %s%s (length: %d)>' % (cn, first_ten, elipses, length) 

    def __reversed__(self):
        return reversed(self._sequence)

    def __str__(self):
        return ''.join(self._sequence)

    def _hamming_distance(self, other):
        """ return the hamming distance to other based on the shorter sequence

            hamming distance is the number of substitutions to convert one
             sequence to the other
        """
        distance = 0
        for s, o in izip(self, other):
            if s != o:
                distance += 1
        return distance
    
    @property
    def alphabet(self):
        """ return the set of characters allowed in the BiologicalSequence
        """
        return self._alphabet

    @property
    def description(self):
        """ return the description of the sequence
        """
        return self._description

    @property
    def gap_alphabet(self):
        """ return the set of gap characters allowed in the BiologicalSequence
        """
        return self._gap_alphabet

    @property
    def identifier(self):
        """ return the identifier of the sequence
        """
        return self._identifier
    
    def count(self, subsequence):
        """ return the number of occurences of subsequence
        """
        return self._sequence.count(subsequence)
 
    def degap(self):
        """ return a new BiologicalSequence with gaps characters removed

            the type, identifier, and description of the result will be the 
             same as self
        """
        result = [e for e in self._sequence if e not in self._gap_alphabet]
        return self.__class__(result, identifier=self._identifier,
                              description=self._description)

    def distance(self, other, distance_fn=_hamming_distance):
        """ return the distance to other using an arbitrary distance function

            distance_fn must take two Sequence objects and is expected to
            return a number (integer or float). for example, see
            BiologicalSequence._hamming_distance.
        """
        return distance_fn(self, other)

    def fraction_diff(self, other):
        """ return fraction of positions that differ 
        
            based on self._hamming_distance between the sequences
        """
        min_edit_dist = self._hamming_distance(other)
        len_shorter = min(len(self), len(other))
        return min_edit_dist / len_shorter
    
    def fraction_same(self, other):
        """ return fraction of positions that are the same 
        
            based on self._hamming_distance between the sequences
        """
        return 1. - self.fraction_diff(other)

    def gap_maps(self):
        """ return tuples mapping positions bw gapped and ungapped seq

            two lists of integers are returned:
             the first is the length of the ungapped sequence, and each entry 
             is the position of that base in the gapped sequence. 
             the second is the length of the gapped sequence, and each entry is
             either None (if that position represents a gap) or the position of
             that base in the ungapped sequence.

            for example:
             BiologicalSequence('-ACCGA-TA-').gap_maps() ==
             ([1,2,3,4,5,7,8],[None,0,1,2,3,4,None,5,6,None])

             because:
             
              01234 56 (spaces are for visual aid)
              ACCGA TA
              ||||| || 
             -ACCGA-TA-
             0123456789

             so... 
             in the first list, position 0 maps to position 1, position 1
             maps to position 2, position 5 maps to position 7, ...
             and in the second list, position 0 doesn't map to anything (so
             it's None), position 1 maps to position 0, ...
        """
        degapped_to_gapped = []
        gapped_to_degapped = []
        non_gap_count = 0
        for i, e in enumerate(self):
            if self.is_gap(e):
                gapped_to_degapped.append(None)
            else:
                gapped_to_degapped.append(non_gap_count)
                degapped_to_gapped.append(i)
                non_gap_count += 1
        return degapped_to_gapped, gapped_to_degapped

    def gap_vector(self):
        """ return a list indicating positions containing gaps 

            for example:
             BiologicalSequence('..ACG--TT-').gap_vector() ==
             [True, True, False, False, False, True, True, False, False, True]
        """
        return map(self.is_gap, self._sequence)

    def unsupported_characters(self):
        """ return set of unsupported characters present in the sequence
        """
        return set(self) - self._alphabet - self._gap_alphabet

    def has_unsupported_characters(self):
        """ return True if unsupported characters are present

            unsupported characters are defined as any characters that are not
            in a BiologicalSequence's alphabet
        """
        all_supported = self._alphabet | self._gap_alphabet
        for e in self:
            if not e in all_supported:
                return True
        return False

    def index(self, subsequence):
        """ return the position where subsequence first occurs
        """
        try:
            return self._sequence.index(subsequence)
        except ValueError:
            raise ValueError(
                "%s is not present in %r." % (subsequence, self))
    
    @classmethod
    def is_gap(cls, char):
        """ return True if char is a gap character
        """
        return char in cls._gap_alphabet

    def is_gapped(self):
        """ return True if any gap characters are in the BiologicalSequence
        """
        for e in self:
            if self.is_gap(e):
                return True
        return False

    def is_valid(self):
        """ return True if the sequence is valid

            validity is defined as not containing any characters outside of
            alphabet and gap_alphabet
        """
        return not self.has_unsupported_characters()

    def to_fasta(self, field_delimiter=" ", terminal_character="\n"):
        """ return the sequence as a fasta-formatted string
          
            terminal_character: the last character to be included in the
             result (if you don't want a trailing newline or other character
             in the result, you can pass terminal_character="")
        """
        if self._description:
            header_line = '%s%s%s' % (self._identifier, field_delimiter,
                                      self._description)
        else:
            header_line = self._identifier

        return '>%s\n%s%s' % (
            header_line, str(self), terminal_character)


def complement(seq_iterator,klass):
    result = []
    for base in seq_iterator:
        try:
            # _complement_map would lose the leading _
            result.append(klass._complement_map[base])
        except KeyError:
            raise BiologicalSequenceError( 
                "Don't know how to complement base %s. Is it in "
                "%s.complement_map?" % (base, klass.__name__))
    return "".join(result)


def reverse_complement(seq_iterator,klass):
    return complement(reversed(seq_iterator),klass)
rc = reverse_complement


def dna_complement(seq_iterator):
    return complement(seq_iterator, DNA)

def dna_reverse_complement(seq_iterator):
    return reverse_complement(seq_iterator, DNA)
dna_rc = dna_reverse_complement


def rna_complment(seq_iterator):
    return complement(seq_iterator, RNA)


def rna_reverse_complement(seq_iterator):
    return reverse_complement(seq_iterator, RNA)
rna_rc = rna_reverse_complement


class NucleotideSequence(BiologicalSequence):
    """ class for representing nucleotide sequences
        
        all uppercase and lowercase IUPAC DNA/RNA characters are supported
    """

    iupac_standard_characters = set("ACGTU")
    iupac_degeneracies = {"R": set("AG"), "Y": set("CTU"), "M": set("AC"), 
                          "K": set("TUG"), "W": set("ATU"), "S": set("GC"), 
                          "B": set("CGTU"), "D": set("AGTU"), 
                          "H": set("ACTU"), "V": set("ACG"), "N": set("ACGTU")}
    iupac_degenerate_characters = set(iupac_degeneracies)
    iupac_characters = iupac_standard_characters | iupac_degenerate_characters
    # complement_map cannot be defined for a generic NucleotideSequence
    # as the complement of 'A' is ambiguous. thanks, nature...
    _complement_map = {}
    _alphabet = iupac_characters | set([c.lower() for c in iupac_characters])

    def _complement(self, seq_iterator):
        """ private method for complementing based on an iterator

            this centralizes the logic for complement and reverse_complement
        """
        result = complement(seq_iterator, self.__class__) 
        return self.__class__(result, self._identifier, self._description)

    @property
    def complement_map(self):
        """ return the mapping of characters to their complements
        """
        return self._complement_map

    def complement(self):
        """ return the complement of the sequence

            raises BiologicalSequence error if there is a character in the
             BiologicalSequence that's not in NucleotideSequence.complement_map
        """
        return self._complement(self)
    
    def is_reverse_complement(self, other):
        """ return True if NucleotideSequences are rev. comp. of one another
            
            raises BiologicalSequence error if there is a character in the
             BiologicalSequence that's not in NucleotideSequence.complement_map
        """
        return self == other.reverse_complement()

    def reverse_complement(self):
        """ return the reverse complement of the sequence

            raises BiologicalSequence error if there is a character in the
             BiologicalSequence that's not in NucleotideSequence.complement_map
        """
        return self._complement(reversed(self))
    rc = reverse_complement


class DNASequence(NucleotideSequence):
    """ class for representing DNA sequences
        
        all uppercase and lowercase IUPAC DNA characters are supported
    """

    iupac_standard_characters = set("ACGT")
    iupac_degeneracies = {"R": set("AG"), "Y": set("CT"), "M": set("AC"), 
                          "K": set("TG"), "W": set("AT"), "S": set("GC"), 
                          "B": set("CGT"), "D": set("AGT"), 
                          "H": set("ACT"), "V": set("ACG"), "N": set("ACGT")}
    iupac_degenerate_characters = set(iupac_degeneracies)
    iupac_characters = iupac_standard_characters | iupac_degenerate_characters
    _complement_map = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y', 'S': 'S',
        'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 
        'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'y': 'r', 'r': 'y', 
        's': 's', 'w': 'w', 'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd',
        'v': 'b', 'n': 'n'}
    _alphabet = iupac_characters | set([c.lower() for c in iupac_characters])

# class is accessible with alternative name for convenience  
DNA = DNASequence


class RNASequence(NucleotideSequence):
    """ class for representing RNA sequences
        
        all uppercase and lowercase IUPAC RNA characters are supported
    """
 
    iupac_standard_characters = set("ACGU")
    iupac_degeneracies = {"R": set("AG"), "Y": set("CU"), "M": set("AC"), 
                          "K": set("UG"), "W": set("AU"), "S": set("GC"), 
                          "B": set("CGU"), "D": set("AGU"), 
                          "H": set("ACU"), "V": set("ACG"), "N": set("ACGU")}
    iupac_degenerate_characters = set(iupac_degeneracies)
    iupac_characters = iupac_standard_characters | iupac_degenerate_characters
    _complement_map = {
        'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y', 'S': 'S',
        'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B',
        'N': 'N', 'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g', 'y': 'r', 'r': 'y',
        's': 's', 'w': 'w', 'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd',
        'v': 'b', 'n': 'n'}
    _alphabet = iupac_characters | set([c.lower() for c in iupac_characters])

# class is accessible with alternative name for convenience  
RNA = RNASequence
