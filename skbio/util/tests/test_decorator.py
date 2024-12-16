# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
import unittest
import inspect
import warnings

from skbio.util import classproperty
from skbio.util._decorator import (
    overrides,
    classproperty,
    classonlymethod,
    deprecated,
    aliased,
    register_aliases,
    params_aliased,
)
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


class TestDeprecated(unittest.TestCase):
    def test_deprecated(self):

        @deprecated("1.0", msg="It will be removed soon.")
        def foo(param):
            """I am a function.

            Parameters
            ----------
            param : str
                I am a parameter.

            """
            pass

        msg = ".. deprecated:: 1.0 It will be removed soon."
        self.assertIn(msg, foo.__doc__)
        with self.assertWarns(DeprecationWarning) as ctx:
            foo("bar")
        msg = "`foo` has been deprecated since 1.0. It will be removed soon."
        self.assertEqual(str(ctx.warning), msg)


class TestAliases(unittest.TestCase):
    def test_aliased(self):

        # apply to a function
        @aliased("bar")
        def foo(param):
            return param

        self.assertEqual(foo._alias[0], "bar")
        self.assertTrue(callable(foo._alias[1]))
        self.assertTrue(foo._alias[1]._skipdoc)
        self.assertEqual(foo.__doc__, f"Alias: ``bar``\n")
        self.assertEqual(foo(42), 42)
        self.assertEqual(foo._alias[1](42), 42)

        # apply to a class method
        class Foo:

            @aliased("bar")
            def foo(cls, param):
                return param

        self.assertEqual(Foo.foo._alias[0], "bar")

    def test_register_aliases(self):

        @register_aliases
        class Foo:

            @aliased("bar", ver="1.0", warn=True)
            def foo(self, param):
                return param

        # class level
        self.assertTrue(hasattr(Foo, "bar"))
        self.assertTrue(Foo.bar._skipdoc)

        # instance level
        f = Foo()
        self.assertTrue(hasattr(f, "bar"))

        # check warning
        self.assertFalse(hasattr(f.foo, "_warned"))
        msg = "`bar` was renamed to `foo` in 1.0."
        with self.assertWarnsRegex(DeprecationWarning, msg):
            f.bar(42)
        self.assertTrue(hasattr(f.foo, "_warned"))

    def test_params_aliased(self):

        @params_aliased([
            ("param1", "alias1", None, False),
            ("param2", "alias2", "1.0", True),
        ])
        def foo(param1, param2):
            """I am a function.

            Parameters
            ----------
            param1 : int
                The first parameter.
            param2 : int
                The second parameter.
            
            Returns
            -------
            int
                Sum of the two parameters.

            """
            return param1 + param2

        self.assertIn("Renamed from ``alias2``.", foo.__doc__)
        self.assertEqual(foo(1, 2), 3)
        self.assertEqual(foo(param1=1, param2=2), 3)
        self.assertEqual(foo(alias1=1, param2=2), 3)
        msg = "`foo`'s parameter `alias2` was renamed to `param2` in 1.0."
        with self.assertWarnsRegex(DeprecationWarning, msg):
            obs = foo(param1=1, alias2=2)
        self.assertEqual(obs, 3)
        with self.assertRaises(TypeError):
            foo(1, 2, alias1=2)


if __name__ == '__main__':
    unittest.main()
