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
    
    def __new__(list,sequence):
        self = sequence
    
    def __getitem__(self,i):
        return self[i]
    
    def __len__(self):
        return len(self)