# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from six import string_types as base_string

import os
import unittest

import numpy.testing as npt

from skbio.util import get_data_path
from skbio.util.testing import IDValidationTests


def test_get_data_path():
    fn = 'parrot'
    path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(path, 'data', fn)
    data_path_2 = get_data_path(fn)
    npt.assert_string_equal(data_path_2, data_path)


class TestIDValidationTestsWithKwargs(IDValidationTests,unittest.TestCase):

    class ExampleObjectKwargs(object):
        def __init__(self, arg1=3, arg2='b', id=''):
            self.arg1 = arg1
            self.arg2 = arg2
            self._id = self.validate_id(id)

        @property
        def id(self):
            return self._id

        def validate_id(self, id):
            if isinstance(id, base_string):
                return id
            raise ValueError('ID is not valid.')

    def setUp(self):
        self.id_cls = self.ExampleObjectKwargs
        self.id_kwargs = {'arg1': 1, 'arg2': 'a'}

    def test_create_instance(self):
        instance = self.create_instance()
        self.assertEqual(instance.arg1, 1)
        self.assertEqual(instance.arg2, 'a')


if __name__ == '__main__':
    import nose
    nose.runmodule()
