# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.standard_library import hooks

from re import compile as re_compile
from collections import Counter, defaultdict, Hashable
from unittest import TestCase, main
from itertools import product, chain

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import euclidean

from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.sequence import SequenceError
from skbio.util._testing import IDValidationTests

with hooks():
    from itertools import zip_longest


if __name__ == "__main__":
    main()
