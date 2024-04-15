# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import inspect
import warnings

from skbio.util import classproperty
from skbio.util._decorator import overrides, classonlymethod
from skbio.util._exception import OverrideError


class TestClassOnlyMethod(unittest.TestCase):
    def test_works_on_class(self):
        class A:
            @classonlymethod
            def example(cls):
                return cls

        self.assertEqual(A.example(), A)

    def test_fails_on_instance(self):
        class A:
            @classonlymethod
            def example(cls):
                pass

        with self.assertRaises(TypeError) as e:
            A().example()

        self.assertIn('A.example', str(e.exception))
        self.assertIn('instance', str(e.exception))

    def test_matches_classmethod(self):
        class A:
            pass

        def example(cls, thing):
            """doc"""

        A.example1 = classmethod(example)
        A.example2 = classonlymethod(example)

        self.assertEqual(A.__dict__['example1'].__func__, example)
        self.assertEqual(A.__dict__['example2'].__func__, example)

        self.assertEqual(A.example1.__doc__, example.__doc__)
        self.assertEqual(A.example2.__doc__, example.__doc__)

        self.assertEqual(A.example1.__name__, example.__name__)
        self.assertEqual(A.example2.__name__, example.__name__)

    def test_passes_args_kwargs(self):
        self.ran_test = False

        class A:
            @classonlymethod
            def example(cls, arg1, arg2, kwarg1=None, kwarg2=None,
                        default=5):
                self.assertEqual(arg1, 1)
                self.assertEqual(arg2, 2)
                self.assertEqual(kwarg1, '1')
                self.assertEqual(kwarg2, '2')
                self.assertEqual(default, 5)
                self.ran_test = True

        A.example(1, *[2], kwarg2='2', **{'kwarg1': '1'})
        self.assertTrue(self.ran_test)


class TestOverrides(unittest.TestCase):
    def test_raises_when_missing(self):
        class A:
            pass

        with self.assertRaises(OverrideError):
            class B(A):
                @overrides(A)
                def test(self):
                    pass

    def test_doc_inherited(self):
        class A:
            def test(self):
                """Docstring"""
                pass

        class B(A):
            @overrides(A)
            def test(self):
                pass

        self.assertEqual(B.test.__doc__, "Docstring")

    def test_doc_not_inherited(self):
        class A:
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
        class Foo:
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
