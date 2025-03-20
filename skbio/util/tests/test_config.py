# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.util.config import get_config, set_config


class TestOptions(TestCase):
    def test_set_config_bad_option(self):
        with self.assertRaisesRegex(ValueError, "Unknown option: 'nonsense'"):
            set_config("nonsense", "asdf")

    def test_set_config_bad_value(self):
        with self.assertRaisesRegex(
            ValueError, "Unsupported value 'asdf' for 'output'"
        ):
            set_config("output", "asdf")

    def test_get_config_bad_option(self):
        with self.assertRaisesRegex(ValueError, "Unknown option: 'frontend'"):
            get_config("frontend")


if __name__ == "__main__":
    main()
