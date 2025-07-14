# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np

from skbio.binaries._util import py_to_bin_random_seed

class UtilTests(TestCase):

    def test_py_to_bin_random_seed(self):

        # returns default value (-1) if seed or generator is None
        obs = py_to_bin_random_seed()
        self.assertEqual(obs, -1)

        # returns given value if seed is integer
        obs = py_to_bin_random_seed(seed_or_generator=42)
        self.assertEqual(obs, 42)

        # raise error if seed is negative
        with self.assertRaisesRegex(ValueError, "seed must be a non-negative number"):
            py_to_bin_random_seed(-1)

        # return an int between 0 and maxint if Generator is provided
        maxint = np.iinfo(np.int32).max + 1
        rng = np.random.default_rng(42)
        obs = py_to_bin_random_seed(rng)
        self.assertTrue(obs >= 0 and obs <= maxint)

        # return an int between 0 and maxint if RandomState is provided
        state = np.random.RandomState(42)
        obs = py_to_bin_random_seed(state)
        self.assertTrue(obs >= 0 and obs <= maxint)

        # raise error if invalid seed
        with self.assertRaisesRegex(ValueError, "Invalid seed. It must be an integer or a random generator instance."):
            py_to_bin_random_seed(0.5)


if __name__ == '__main__':
    main()