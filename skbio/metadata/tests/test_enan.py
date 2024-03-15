# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from skbio.metadata._enan import get_payload_from_nan, make_nan_with_payload


class TestNanPayloads(unittest.TestCase):
    def test_normal_nan(self):
        normal_nan = float('nan')

        payload, namespace = get_payload_from_nan(normal_nan)

        self.assertIs(payload, None)
        self.assertIs(namespace, None)

    def test_roundtrip_payload(self):
        for namespace in range(0, 256):
            for payload in range(-50, 500):
                nan = make_nan_with_payload(payload, namespace)
                new_payload, new_namespace = get_payload_from_nan(nan)

                self.assertEqual(namespace, new_namespace)
                self.assertEqual(payload, new_payload)

                self.assertNotEqual(nan, nan)

    def test_user_namespace_default(self):
        nan = make_nan_with_payload(42)

        payload, namespace = get_payload_from_nan(nan)

        self.assertEqual(42, payload)
        self.assertEqual(255, namespace)

        self.assertNotEqual(nan, nan)
