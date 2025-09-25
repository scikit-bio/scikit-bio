# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import importlib

from unittest import mock

from skbio.util import get_package


class TestGetPackage(TestCase):
    def test_get_numpy(self):
        """Import NumPy, a mandatory dependency."""
        obs = get_package("numpy")
        self.assertIsNotNone(obs)
        self.assertIsInstance(obs.ndarray, type)

    def test_get_numpy_testing(self):
        """Import NumPy's testing module."""
        obs = get_package("numpy.testing")
        self.assertIsNotNone(obs)
        self.assertTrue(callable(obs.assert_array_equal))

    def test_get_polars_import_error(self):
        """Import Polars, pretending it doesn't exist."""
        msg = (
            'Optional dependency "polars" is not found. '
            "Install it to use relevant functionalities."
        )
        with mock.patch.object(importlib, "import_module", side_effect=ImportError):
            with self.assertRaises(ImportError) as cm:
                get_package("polars")
            self.assertEqual(str(cm.exception), msg)


if __name__ == "__main__":
    main()
