# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from skbio.util import TestRunner

test = TestRunner(__file__).test

if __name__ == '__main__':
    if test():
        sys.exit(0)
    else:
        sys.exit(1)
