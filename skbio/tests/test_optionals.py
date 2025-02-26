# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import importlib

from unittest import mock

from skbio._optionals import _get_package


class TestGetPackage(unittest.TestCase):
    def test_get_polars_import_error(self):
        with mock.patch.object(importlib, "import_module", side_effect=ImportError):
            with self.assertRaises(ImportError) as c:
                _get_package("polars")
            self.assertIn(
                "Using the polars backend requires the polars package to be installed.",
                str(c.exception),
            )


if __name__ == "__main__":
    unittest.main()
