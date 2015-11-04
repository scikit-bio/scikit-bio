r"""
Gene ontology enrichment analysis (:mod:`skbio.goenrich`)
==============================================================

.. currentmodule:: skbio.goenrich

Utility functions for programmatic gene ontology (GO) enrichment
analysis.

Enrichment
----------

.. currentmodule:: skbio.goenrich.enrich
.. autosummary::
    
   analyze 

Export
------

.. autosummary::

    export

Tools
-----

.. autosummary::

    tools

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util import TestRunner
test = TestRunner(__file__).test

import skbio.goenrich.enrich
import skbio.goenrich.export
import skbio.goenrich.tools
