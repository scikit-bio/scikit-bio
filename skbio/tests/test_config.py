# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose

from skbio import get_config, set_config, reset_config
from skbio.stats.ordination import pcoa, center_distance_matrix
from skbio._config import _resolve_engine
from skbio.util import numba_code


class TestOptions(TestCase):
    def setUp(self):
        self.original = get_config()
        reset_config()

    def tearDown(self):
        for option, value in self.original.items():
            set_config(option, value)

    def test_set_config_bad_option(self):
        with self.assertRaisesRegex(KeyError, "Unknown option: 'nonsense'."):
            set_config("nonsense", "asdf")

    def test_set_config_bad_value(self):
        with self.assertRaisesRegex(
            ValueError, "Unsupported value 'asdf' for 'table_output'."
        ):
            set_config("table_output", "asdf")

    def test_get_config_bad_option(self):
        with self.assertRaisesRegex(KeyError, "Unknown option: 'frontend'."):
            get_config("frontend")

    def test_engine_default_is_cython(self):
        self.assertEqual(get_config("compute_engine"), "cython")

    def test_set_engine_valid(self):
        set_config("compute_engine", "numba")
        self.assertEqual(get_config("compute_engine"), "numba")

    def test_set_engine_bad_value(self):
        with self.assertRaisesRegex(
            ValueError, "Unsupported value 'julia' for 'compute_engine'."
        ):
            set_config("compute_engine", "julia")

    def test_set_engine_fast(self):
        set_config("compute_engine", "fast")
        self.assertEqual(get_config("compute_engine"), "fast")

    def test_get_all_is_copy(self):
        set_config("compute_engine", "fast")
        options = get_config()
        self.assertEqual(options, {"compute_engine": "fast", "table_output": "pandas"})
        options["compute_engine"] = "numba"
        options["unknown"] = "value"
        self.assertEqual(get_config("compute_engine"), "fast")
        self.assertNotIn("unknown", get_config())

    def test_reset_one(self):
        set_config("compute_engine", "fast")
        set_config("table_output", "numpy")
        reset_config("compute_engine")
        self.assertEqual(get_config("compute_engine"), "cython")
        self.assertEqual(get_config("table_output"), "numpy")
        reset_config("table_output")
        self.assertEqual(get_config("table_output"), "pandas")

    def test_reset_all(self):
        set_config("compute_engine", "numba")
        set_config("table_output", "polars")
        reset_config()
        self.assertEqual(get_config(), {
            "compute_engine": "cython", "table_output": "pandas"})
        reset_config()
        self.assertEqual(get_config("compute_engine"), "cython")

    def test_invalid_option_does_not_change_settings(self):
        for option in ("unknown", "engine"):
            before = get_config()
            for action in (get_config, reset_config):
                with self.assertRaisesRegex(KeyError, "Unknown option"):
                    action(option)
            with self.assertRaisesRegex(KeyError, "Unknown option"):
                set_config(option, "numba")
            self.assertEqual(get_config(), before)

    def test_invalid_value_does_not_change_settings(self):
        before = get_config()
        with self.assertRaises(ValueError):
            set_config("compute_engine", "invalid")
        self.assertEqual(get_config(), before)


class TestResolveEngine(TestCase):
    def setUp(self):
        self._original = get_config("compute_engine")

    def tearDown(self):
        set_config("compute_engine", self._original)

    def test_none_uses_global_default(self):
        set_config("compute_engine", "cython")
        self.assertEqual(_resolve_engine(None, ("cython", "numba")), "cython")

    def test_explicit_cython(self):
        self.assertEqual(_resolve_engine("cython", ("cython", "numba")), "cython")

    def test_unsupported_value_raises(self):
        with self.assertRaisesRegex(ValueError, "engine='julia' is not supported"):
            _resolve_engine("julia", ("cython", "numba"))

    def test_engine_not_in_supported_raises(self):
        with self.assertRaisesRegex(ValueError, "engine='numba' is not supported"):
            _resolve_engine("numba", ("cython",))

    def test_numba_requested_but_absent_raises(self):
        # Simulate Numba not being installed.
        import builtins

        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "numba":
                raise ImportError("no numba")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=fake_import):
            with self.assertRaisesRegex(ImportError, "requires the optional numba"):
                _resolve_engine("numba", ("cython", "numba"))

    def test_fast_resolves_to_what_the_caller_names(self):
        # A target the resolver could not have arrived at on its own shows that
        # the caller's value is what gets used, and needs no optional
        # dependency to check. The name is deliberately nonsense so that it
        # cannot be read as an engine scikit-bio might one day support.
        with self.assertRaisesRegex(
            ValueError, "engine='SantaGoesSkiing' is not supported"
        ):
            _resolve_engine("fast", ("cython", "numba"), fast="SantaGoesSkiing")
        self.assertEqual(
            _resolve_engine("fast", ("cython", "numba"), fast="cython"), "cython"
        )

    @numba_code
    def test_fast_resolves_to_numba_when_that_is_the_target(self):
        # The production case: every wired call site passes "numba" when numba
        # imports. Marked, since resolving to it imports numba.
        self.assertEqual(
            _resolve_engine("fast", ("cython", "numba"), fast="numba"), "numba"
        )

    def test_fast_without_a_target_falls_back_to_the_default(self):
        # A function with nothing faster to offer does not pass fast=, and
        # engine="fast" then has to be a no-op rather than an error.
        set_config("compute_engine", "cython")
        self.assertEqual(_resolve_engine("fast", ("cython", "numba")), "cython")

    def test_fast_is_resolved_after_the_global_default(self):
        set_config("compute_engine", "fast")
        self.assertEqual(
            _resolve_engine(None, ("cython", "numba"), fast="cython"), "cython"
        )

    def test_fast_stays_a_no_op_when_the_default_is_itself_fast(self):
        set_config("compute_engine", "fast")
        self.assertEqual(_resolve_engine(None, ("cython", "numba")), "cython")
        self.assertEqual(_resolve_engine("fast", ("cython", "numba")), "cython")

    def test_explicit_engine_overrides_global(self):
        set_config("compute_engine", "numba")
        self.assertEqual(_resolve_engine("cython", ("cython", "numba")), "cython")

    @numba_code
    def test_global_fast_selects_numba(self):
        set_config("compute_engine", "fast")
        self.assertEqual(
            _resolve_engine(None, ("cython", "numba"), fast="numba"), "numba"
        )

    @numba_code
    def test_pcoa_global_numba_bypasses_binaries(self):
        data = np.array([[0., 1., 2.], [1., 0., 1.], [2., 1., 0.]])
        set_config("compute_engine", "numba")
        module = "skbio.stats.ordination._principal_coordinate_analysis"
        with patch(module + "._skbb_pcoa_fsvd_available") as available:
            result = pcoa(data, method="fsvd", dimensions=1, seed=0)
            available.assert_not_called()
        assert_allclose(result.eigvals, [2.], atol=1e-12)

    def test_center_global_fast(self):
        data = np.array([[0., 2.], [2., 0.]])
        set_config("compute_engine", "fast")
        result = center_distance_matrix(data)
        assert_allclose(result, [[1., -1.], [-1., 1.]])
        assert_allclose(data, [[0., 2.], [2., 0.]])

    def test_fast_target_still_checked_against_supported(self):
        with self.assertRaisesRegex(ValueError, "engine='numba' is not supported"):
            _resolve_engine("fast", ("cython",), fast="numba")


if __name__ == "__main__":
    main()
