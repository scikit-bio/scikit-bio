# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function


class BiologicalSequenceError(Exception):
    """General error for biological sequence validation failures."""
    pass


class GeneticCodeError(Exception):
    """Base class exception used by the GeneticCode class"""
    pass


class GeneticCodeInitError(ValueError, GeneticCodeError):
    """Exception raised by the GeneticCode class upon a bad initialization"""
    pass


class InvalidCodonError(KeyError, GeneticCodeError):
    """Exception raised by the GeneticCode class if __getitem__ fails"""
    pass
