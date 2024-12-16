# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import warnings

from skbio.util._warning import _warn_once, _warn_deprecated, _warn_param_deprecated


class TestWarning(unittest.TestCase):
    def test_warn_once(self):

        def foo(param):
            pass

        # function warning
        self.assertFalse(hasattr(foo, "_warned"))
        wtype = FutureWarning
        msg = "`foo` will become `bar` in 2.0."
        with self.assertWarns(wtype) as ctx:
            _warn_once(foo, wtype, msg)
        self.assertEqual(str(ctx.warning), msg)
        self.assertTrue(hasattr(foo, "_warned"))
        with self.assertRaises(AssertionError):
            self.assertWarns(wtype, _warn_once, foo, wtype, msg)

        # parameter warning
        self.assertFalse(hasattr(foo, "_warned_params"))
        wtype = DeprecationWarning
        msg = "`param` is deprecated as of 3.0."
        with self.assertWarns(wtype) as ctx:
            _warn_once(foo, wtype, msg, "param")
        self.assertEqual(str(ctx.warning), msg)
        self.assertIn("param", foo._warned_params)
        with self.assertRaises(AssertionError):
            self.assertWarns(wtype, _warn_once, foo, wtype, msg, "param")


    def test_warn_deprecated(self):

        def foo():
            pass

        self.assertFalse(hasattr(foo, "_warned"))
        with self.assertWarns(DeprecationWarning) as ctx:
            _warn_deprecated(foo)
        exp = "`foo` is deprecated."
        self.assertEqual(str(ctx.warning), exp)
        self.assertTrue(hasattr(foo, "_warned"))
        with self.assertRaises(AssertionError):
            self.assertWarns(DeprecationWarning, _warn_deprecated, foo)

        def foo():
            pass

        with self.assertWarns(DeprecationWarning) as ctx:
            _warn_deprecated(foo, ver="1.0")
        exp = "`foo` has been deprecated since 1.0."
        self.assertEqual(str(ctx.warning), exp)

        def foo():
            pass

        with self.assertWarns(DeprecationWarning) as ctx:
            _warn_deprecated(foo, ver="1.0", msg="Use `bar` instead.")
        exp = "`foo` has been deprecated since 1.0. Use `bar` instead."
        self.assertEqual(str(ctx.warning), exp)

        def foo():
            pass

        msg="Use `bar` instead of `foo`."
        with self.assertWarns(DeprecationWarning) as ctx:
            _warn_deprecated(foo, ver="1.0", msg=msg, append=False)
        self.assertEqual(str(ctx.warning), msg)

    def test_warn_param_deprecated(self):

        def foo(param1, param2):
            pass

        self.assertFalse(hasattr(foo, "_warned_params"))

        # parameter 1
        with self.assertWarns(DeprecationWarning) as ctx:
            _warn_param_deprecated(foo, "param1")
        exp = "`foo`'s parameter `param1` is deprecated."
        self.assertEqual(str(ctx.warning), exp)
        self.assertIn("param1", foo._warned_params)
        with self.assertRaises(AssertionError):
            self.assertWarns(DeprecationWarning, _warn_param_deprecated, foo, "param1")

        # parameter 2
        with self.assertWarns(DeprecationWarning) as ctx:
            _warn_param_deprecated(foo, "param2", ver="1.0")
        exp = "`foo`'s parameter `param2` has been deprecated since 1.0."
        self.assertEqual(str(ctx.warning), exp)
        self.assertIn("param2", foo._warned_params)
        with self.assertRaises(AssertionError):
            self.assertWarns(DeprecationWarning, _warn_param_deprecated, foo, "param2")


if __name__ == '__main__':
    unittest.main()
