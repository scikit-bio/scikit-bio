"""
Warnings (:mod:`skbio.core.warning`)
====================================

.. currentmodule:: skbio.core.warning

This module defines custom warning classes used throughout the core scikit-bio
codebase.

Warnings
--------

.. autosummary::
   :toctree: generated/

   EfficiencyWarning

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function


class EfficiencyWarning(Warning):
    pass
