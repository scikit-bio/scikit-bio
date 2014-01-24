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

LazyDeveloperError = NotImplementedError

class BiologicalSequenceError(Exception):
    pass

class BiologicalSequence(Sequence):
    """ Base class for biological sequences """
    
    _alphabet = set()
    _gap_alphabet = set('-.')

    def __init__(self, sequence, identifier="", description=""):

        self._sequence = ''.join(sequence)
        self._identifier = identifier
        self._description = description
 
    def __contains__(self, other):
        return other in self._sequence
   
    def __eq__(self, other):
        """ define equality
            
            BiologicalSequences are equal if their sequence is the same and
             they are the same type
        """
        if type(self) != type(other):
            return False
        elif self._sequence != other._sequence:
            return False
        else:
            return True

    def __getitem__(self, i):
        return self._sequence[i]
 
    def __hash__(self):
        return hash(self._sequence)
   
    def __iter__(self):
        return iter(self._sequence)

    def __len__(self):
        return len(self._sequence)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        first_ten = self._sequence[:10]
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
        return str(self._sequence)

    def _minimum_edit_distance(self,other):
        distance = 0
        for s, o in zip(self,other):
            if s != o:
                distance +=1
        return distance
    
    @property
    def Alphabet(self):
        return self._alphabet

    @property
    def Description(self):
        return self._description

    @property
    def GapAlphabet(self):
        return self._gap_alphabet

    @property
    def Identifier(self):
        return self._identifier
    
    def count(self, char):
        return self._sequence.count(char)
 
    def degap(self):
        result = [e for e in self._sequence if e not in self._gap_alphabet]
        return self.__class__(result, identifier=self._identifier,
                              description=self._description)

    def distance(self,other,distance_fn=_minimum_edit_distance):
        """
        """
        return distance_fn(self,other)

    def fractionDiff(self,other):
        min_edit_dist = self._minimum_edit_distance(other)
        len_shorter = min(len(self),len(other))
        return min_edit_dist / len_shorter
    
    def fractionSame(self,other):
        return 1. - self.fractionDiff(other)

    def gapMaps(self):
        """
        """
        degapped_to_gapped = []
        gapped_to_degapped = []
        non_gap_count = 0
        for i,e in enumerate(self):
            if self.isGap(e):
                gapped_to_degapped.append(None)
            else:
                gapped_to_degapped.append(non_gap_count)
                degapped_to_gapped.append(i)
                non_gap_count += 1
        return degapped_to_gapped, gapped_to_degapped

    def gapVector(self):
        return [e in self._gap_alphabet for e in self._sequence]

    def getUnsupportedCharacters(self):
        """ return set of unsupported characters present in the sequence
        """
        return set(self) - self._alphabet - self._gap_alphabet

    def hasUnsupportedCharacters(self):
        """ return True if unsupported characters are present

            unsupported characters are defined as any characters that are not
            in a BiologicalSequence's Alphabet
        """
        return len(self.getUnsupportedCharacters()) > 0

    def index(self, char):
        return self._sequence.index(char)

    def isGap(self, char):
        return char in self._gap_alphabet

    def isGapped(self):
        for e in self:
            if e in self._gap_alphabet:
                return True
        return False

    def toFasta(self, field_delimiter = " ", terminal_character="\n"):
        """ return the sequence as a fasta-formatted string
          
            terminal_character: the last character to be included in the
             string (default: \n (i.e., newline); if you don't want a trailing
             newline in the string, you can pass terminal_character="")
        """
        if self._identifier != "" and self._description != "":
            header_line = "%s%s%s" % (
             self._identifier, field_delimiter, self._description)
        elif self._identifier == "" and self._description == "":
            header_line = ""
        elif self._identifier:
            header_line = self._identifier
        elif self._description:
            header_line = "%s%s" % (field_delimiter, self._description)
        else:
            # we've exhausted the possibilities - it shouldn't be 
            # possible to get here, but just in case...
            raise BiologicalSequenceError(
             "Can't construct header line in BiologicalSequence.toFasta().")

        return '>%s\n%s%s' % (
         header_line, self._sequence, terminal_character)


class NucleotideSequence(BiologicalSequence):
    """ Base class for nucleotide sequences """

    _complement_map = {}
    _alphabet = set('ACGTURYMKWSBDHVNacgturymkwsbdhvn')

    def _complement(self, seq_iterator):
        result = []
        for base in seq_iterator:
            try:
                result.append(self._complement_map[base])
            except KeyError:
                raise BiologicalSequenceError( 
                 "Don't know how to complement base %s. "
                 "Is it a known base in your nucleotide alphabet?" % base)
        return self.__class__(result, self._identifier, self._description)

    def complement(self):
        return self._complement(self)

    def reverse_complement(self):
        return self._complement(reversed(self))
    
    def isReverseComplement(self,other):
        return self == other.reverse_complement()

class DNASequence(NucleotideSequence):
 
    _complement_map = {
     'A':'T', 'T':'A', 'G':'C', 'C':'G', 'Y':'R', 'R':'Y', 'S':'S',
     'W':'W', 'K':'M', 'M':'K', 'B':'V', 'D':'H', 'H':'D', 'V':'B', 'N':'N',
     'a':'t', 't':'a', 'g':'c', 'c':'g', 'y':'r', 'r':'y', 's':'s',
     'w':'w', 'k':'m', 'm':'k', 'b':'v', 'd':'h', 'h':'d', 'v':'b', 'n':'n'}
    _alphabet = set('ACGTRYMKWSBDHVNacgtrymkwsbdhvn')

# class is accessible with alternative capitalization scheme for convenience  
DnaSequence = DNASequence

class RNASequence(NucleotideSequence):
  
    _complement_map = {
     'A':'U', 'U':'A', 'G':'C', 'C':'G', 'Y':'R', 'R':'Y', 'S':'S',
     'W':'W', 'K':'M', 'M':'K', 'B':'V', 'D':'H', 'H':'D', 'V':'B', 'N':'N',
     'a':'u', 'u':'a', 'g':'c', 'c':'g', 'y':'r', 'r':'y', 's':'s',
     'w':'w', 'k':'m', 'm':'k', 'b':'v', 'd':'h', 'h':'d', 'v':'b', 'n':'n'}
    _alphabet = set('ACGURYMKWSBDHVNacgurymkwsbdhvn')

# class is accessible with alternative capitalization scheme for convenience  
RnaSequence = RNASequence
