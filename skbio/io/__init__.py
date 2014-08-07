r"""
"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


from ._registry import (write, read, guess_format, get_writer, get_reader,
                        get_identifier, list_write_formats, list_read_formats,
                        register_writer, register_reader, register_identifier)

__all__ = ['write', 'read', 'guess_format', 'get_writer', 'get_reader',
           'get_identifier', 'list_write_formats', 'list_read_formats',
           'register_writer', 'register_reader', 'register_identifier']

from numpy.testing import Tester
test = Tester().test
