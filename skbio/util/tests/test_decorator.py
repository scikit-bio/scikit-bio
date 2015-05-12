# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
import unittest

from skbio.util import classproperty, overrides
from skbio.util._exception import OverrideError


class TestOverrides(unittest.TestCase):
    def test_raises_when_missing(self):
        class A(object):
            pass

        with self.assertRaises(OverrideError):
            class B(A):
                @overrides(A)
                def test(self):
                    pass

    def test_doc_inherited(self):
        class A(object):
            def test(self):
                """Docstring"""
                pass

        class B(A):
            @overrides(A)
            def test(self):
                pass

        self.assertEqual(B.test.__doc__, "Docstring")

    def test_doc_not_inherited(self):
        class A(object):
            def test(self):
                """Docstring"""
                pass

        class B(A):
            @overrides(A)
            def test(self):
                """Different"""
                pass

        self.assertEqual(B.test.__doc__, "Different")


class TestClassProperty(unittest.TestCase):
    def test_getter_only(self):
        class Foo(object):
            _foo = 42

            @classproperty
            def foo(cls):
                return cls._foo

        # class-level getter
        self.assertEqual(Foo.foo, 42)

        # instance-level getter
        f = Foo()
        self.assertEqual(f.foo, 42)

        with self.assertRaises(AttributeError):
            f.foo = 4242


if __name__ == '__main__':
    unittest.main()
