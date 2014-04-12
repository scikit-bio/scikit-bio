# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import inspect

def get_data_path(fn):
    """Return path to filename ``fn`` in the data folder.

    Useful to load files in tests/data/* when testing."""
    # getouterframes returns a list of tuples: the second tuple
    # contains info about the caller, and the second element is its
    # filname
    callers_filename = inspect.getouterframes(inspect.currentframe())[1][1]
    path = os.path.dirname(os.path.abspath(callers_filename))
    data_path = os.path.join(path, 'data', fn)
    return data_path
