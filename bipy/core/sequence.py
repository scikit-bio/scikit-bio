#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import Sequence

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

    def toFasta(self, terminal_character="\n"):
        """ return the sequence as a fasta-formatted string
          
            terminal_character: the last character to be included in the
             string (default: \n (i.e., newline); if you don't want a trailing
             newline in the string, you can pass terminal_character="")
        """
        return '>%s %s\n%s%s' % (self._identifier, self._description, 
                                 self._sequence, terminal_character)


