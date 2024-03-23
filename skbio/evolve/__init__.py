r"""Evolutionary statistics (:mod:`skbio.evolve`)
=============================================

.. currentmodule:: skbio.evolve

This package contains statistics pertaining to phylogenies and evolution.

Cophylogenetic methods
----------------------

These functions test for correlation between phylogenies or representations of
evolutionary distance (for example, genetic distance matrices).

Functions
^^^^^^^^^

.. autosummary::
   :toctree: generated/

   hommola_cospeciation

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._hommola import hommola_cospeciation

__all__ = ["hommola_cospeciation"]
