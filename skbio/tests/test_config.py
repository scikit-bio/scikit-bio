# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio._config import get_option, set_option


class TestOptions(TestCase):
    def test_set_option_bad_option(self):
        with self.assertRaisesRegex(ValueError, "Unknown option: 'nonsense'"):
            set_option("nonsense", "asdf")

    def test_set_option_bad_value(self):
        with self.assertRaisesRegex(
            ValueError, "Unsupported value 'asdf' for 'tabular_backend'"
        ):
            set_option("tabular_backend", "asdf")

    def test_get_option_bad_option(self):
        with self.assertRaisesRegex(ValueError, "Unknown option: 'frontend'"):
            get_option("frontend")


if __name__ == "__main__":
    main()
