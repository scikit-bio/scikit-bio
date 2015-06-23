# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
import unittest
import inspect
import warnings

from skbio.util import classproperty, overrides
from skbio.util._decorator import stable, experimental, deprecated
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


class TestStable(unittest.TestCase):

    def _get_f(self, as_of):
        @stable(as_of=as_of)
        def f(x, y=42):
            """ Add 42, or something else, to x.

                Parameters
                ----------
                x : int, x
                y : int, optional

            """
            return x + y
        return f

    def test_function_output(self):
        f = self._get_f('0.1.0')
        self.assertEqual(f(1), 43)

    def test_function_docstring(self):
        f = self._get_f('0.1.0')
        self.assertTrue('State: Stable as of 0.1.0.' in f.__doc__)
        f = self._get_f('0.1.1')
        self.assertTrue('State: Stable as of 0.1.1.' in f.__doc__)

    def test_function_signature(self):
        f = self._get_f('0.1.0')
        expected = inspect.ArgSpec(
            args=['x', 'y'], varargs=None, keywords=None, defaults=(42,))
        self.assertEqual(inspect.getargspec(f), expected)
        self.assertEqual(f.__name__, 'f')

    def test_missing_kwarg(self):
        self.assertRaises(ValueError, stable)
        self.assertRaises(ValueError, stable, '0.1.0')


class TestExperimental(unittest.TestCase):

    def _get_f(self, as_of):
        @experimental(as_of=as_of)
        def f(x, y=42):
            """ Add 42, or something else, to x.

                Parameters
                ----------
                x : int, x
                y : int, optional

            """
            return x + y
        return f

    def test_function_output(self):
        f = self._get_f('0.1.0')
        self.assertEqual(f(1), 43)

    def test_function_docstring(self):
        f = self._get_f('0.1.0')
        self.assertTrue(
            'State: Experimental as of 0.1.0.' in f.__doc__)
        f = self._get_f('0.1.1')
        self.assertTrue(
            'State: Experimental as of 0.1.1.' in f.__doc__)

    def test_function_signature(self):
        f = self._get_f('0.1.0')
        expected = inspect.ArgSpec(
            args=['x', 'y'], varargs=None, keywords=None, defaults=(42,))
        self.assertEqual(inspect.getargspec(f), expected)
        self.assertEqual(f.__name__, 'f')

    def test_missing_kwarg(self):
        self.assertRaises(ValueError, experimental)
        self.assertRaises(ValueError, experimental, '0.1.0')


class TestDeprecated(unittest.TestCase):

    def _get_f(self, as_of, until, reason):
        @deprecated(as_of=as_of, until=until, reason=reason)
        def f(x, y=42):
            """ Add 42, or something else, to x.

                Parameters
                ----------
                x : int, x
                y : int, optional

            """
            return x + y
        return f

    def test_function_output(self):
        f = self._get_f('0.1.0', until='0.1.4',
                        reason='You should now use skbio.g().')
        self.assertEqual(f(1), 43)

    def test_deprecation_warning(self):
        f = self._get_f('0.1.0', until='0.1.4',
                        reason='You should now use skbio.g().')
        # adapted from SO example here: http://stackoverflow.com/a/3892301
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            f(1)
            self.assertTrue(issubclass(w[-1].category, DeprecationWarning))
            expected_str = "is deprecated as of scikit-bio version 0.1.0"
            self.assertTrue(expected_str in str(w[-1].message))

    def test_function_docstring(self):
        f = self._get_f('0.1.0', until='0.1.4',
                        reason='You should now use skbio.g().')
        e1 = ("        .. note:: Deprecated as of 0.1.0 for "
              "removal in 0.1.4. You should now")
        e2 = "                  use skbio.g()."
        self.assertTrue(e1 in f.__doc__)
        self.assertTrue(e2 in f.__doc__)

        f = self._get_f('0.1.1', until='0.1.5',
                        reason='You should now use skbio.h().')
        e1 = ("        .. note:: Deprecated as of 0.1.1 for "
              "removal in 0.1.5. You should now")
        e2 = "                  use skbio.h()."
        self.assertTrue(e1 in f.__doc__)
        self.assertTrue(e2 in f.__doc__)

    def test_function_signature(self):
        f = self._get_f('0.1.0', until='0.1.4',
                        reason='You should now use skbio.g().')
        expected = inspect.ArgSpec(
            args=['x', 'y'], varargs=None, keywords=None, defaults=(42,))
        self.assertEqual(inspect.getargspec(f), expected)
        self.assertEqual(f.__name__, 'f')

    def test_missing_kwarg(self):
        self.assertRaises(ValueError, deprecated)
        self.assertRaises(ValueError, deprecated, '0.1.0')
        self.assertRaises(ValueError, deprecated, as_of='0.1.0')
        self.assertRaises(ValueError, deprecated, as_of='0.1.0', until='0.1.4')

if __name__ == '__main__':
    unittest.main()
