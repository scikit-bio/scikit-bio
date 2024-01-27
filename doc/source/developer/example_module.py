r"""
Numbers (:mod:`skbio.numbers`)
==============================

.. currentmodule:: skbio.numbers

Numbers holds a sequence of numbers, and defines several statistical
operations (mean, stdev, etc.) FrequencyDistribution holds a mapping from
items (not necessarily numbers) to counts, and defines operations such as
Shannon entropy and frequency normalization.


Classes
-------

.. autosummary::
    :toctree: generated/

    Numbers

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from random import choice, random

import numpy as np
from utils import indices


class Numbers(list):
    pass


class FrequencyDistribution(dict):
    pass
