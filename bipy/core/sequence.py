#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import Sequence

class BiologicalSequenceError(Exception):
    pass

class BiologicalSequence(Sequence):
    """ Base class for biological sequences """
    
    def __init__(self, sequence, identifier="", description=""):
        self._sequence = sequence
        self._identifier = identifier
        self._description = description
    
    def __getitem__(self, i):
        return self._sequence[i]
    
    def __len__(self):
        return len(self._sequence)

    def __str__(self):
        return str(self._sequence)

    def __repr__(self):
        return str(self)

    def __iter__(self):
        return iter(self._sequence)

    def __reversed__(self):
        return reversed(self._sequence)

    def toFasta(self, terminal_character="\n"):
        """ return the sequence as a fasta-formatted string
          
            terminal_character: the last character to be included in the
             string (default: \n (i.e., newline); if you don't want a trailing
             newline in the string, you can pass terminal_character="")
        """
        return '>%s %s\n%s%s' % (self._identifier, self._description, 
                                 self._sequence, terminal_character)


class NucleotideSequence(BiologicalSequence):
    """ Base class for nucleotide sequences """

    _complement_map = {}

    def _complement(self, seq_iterator):
        result = []
        for base in seq_iterator:
            try:
                result.append(self._complement_map[base])
            except KeyError:
                raise BiologicalSequenceError( 
                 "Don't know how to complement base %s. "
                 "Is it a known base in your nucleotide alphabet?" % base)
        return NucleotideSequence(result, self._identifier, self._description)

    def complement(self):
        return self._complement(self)

    def reverse_complement(self):
        return self._complement(reversed(self))

class DNASequence(NucleotideSequence):
    
    _complement_map = {
     'A':'T', 'T':'A', 'G':'C', 'C':'G', 'a':'t', 't':'a', 'g':'c', 'c':'g'} 

# class is accessible with alternative capitalization scheme for convenience  
DnaSequence = DNASequence

class RNASequence(NucleotideSequence):
    
    _complement_map = {
     'A':'U', 'U':'A', 'G':'C', 'C':'G', 'a':'u', 'u':'a', 'g':'c', 'c':'g'} 

# class is accessible with alternative capitalization scheme for convenience  
RnaSequence = RNASequence
