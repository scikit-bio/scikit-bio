# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import numpy as np
import numpy.testing as npt

from skbio.sequence.tpm import jc69, k2p, f81, hky85, tn93


class TestJC69(TestCase):
    def test_instance(self):
        check = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        npt.assert_array_equal(jc69(0).data, check)

        check = np.array(
            [
                [0.6350628, 0.1216457, 0.1216457, 0.1216457],
                [0.1216457, 0.6350628, 0.1216457, 0.1216457],
                [0.1216457, 0.1216457, 0.6350628, 0.1216457],
                [0.1216457, 0.1216457, 0.1216457, 0.6350628],
            ]
        )
        npt.assert_array_equal(jc69(0.5).data.round(7), check)

        check = np.array(
            [
                [0.4476979, 0.1841007, 0.1841007, 0.1841007],
                [0.1841007, 0.4476979, 0.1841007, 0.1841007],
                [0.1841007, 0.1841007, 0.4476979, 0.1841007],
                [0.1841007, 0.1841007, 0.1841007, 0.4476979],
            ]
        )
        npt.assert_array_equal(jc69(1).data.round(7), check)

        check = np.array(
            [
                [0.2500012, 0.2499996, 0.2499996, 0.2499996],
                [0.2499996, 0.2500012, 0.2499996, 0.2499996],
                [0.2499996, 0.2499996, 0.2500012, 0.2499996],
                [0.2499996, 0.2499996, 0.2499996, 0.2500012],
            ]
        )
        npt.assert_array_equal(jc69(10).data.round(7), check)

        check = np.array(
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
            ]
        )
        npt.assert_array_equal(jc69(14).data.round(7), check)

        with self.assertRaises(ValueError):
            jc69(-1)

        with self.assertRaises(TypeError):
            jc69(1, seqtype="Protein")


class TestK2P(TestCase):
    def test_instance(self):
        check = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        npt.assert_array_equal(k2p(0, kappa=0.5).data, check)
        npt.assert_array_equal(k2p(0, kappa=0.5).data, k2p(0, kappa=1).data)

        check = np.array(
            [
                [0.6816196, 0.1216457, 0.0750889, 0.1216457],
                [0.1216457, 0.6816196, 0.1216457, 0.0750889],
                [0.0750889, 0.1216457, 0.6816196, 0.1216457],
                [0.1216457, 0.0750889, 0.1216457, 0.6816196],
            ]
        )
        npt.assert_array_equal(k2p(0.5, kappa=0.5).data.round(7), check)

        check = np.array(
            [
                [0.6350628, 0.1216457, 0.1216457, 0.1216457],
                [0.1216457, 0.6350628, 0.1216457, 0.1216457],
                [0.1216457, 0.1216457, 0.6350628, 0.1216457],
                [0.1216457, 0.1216457, 0.1216457, 0.6350628],
            ]
        )
        npt.assert_array_equal(k2p(0.5, kappa=1).data.round(7), check)
        npt.assert_allclose(k2p(0.5, kappa=1).data, jc69(0.5).data)

        check = np.array(
            [
                [0.499839, 0.1841007, 0.1319596, 0.1841007],
                [0.1841007, 0.499839, 0.1841007, 0.1319596],
                [0.1319596, 0.1841007, 0.499839, 0.1841007],
                [0.1841007, 0.1319596, 0.1841007, 0.499839],
            ]
        )
        npt.assert_array_equal(k2p(1.0, kappa=0.5).data.round(7), check)

        check = np.array(
            [
                [0.2500231, 0.2499996, 0.2499777, 0.2499996],
                [0.2499996, 0.2500231, 0.2499996, 0.2499777],
                [0.2499777, 0.2499996, 0.2500231, 0.2499996],
                [0.2499996, 0.2499777, 0.2499996, 0.2500231],
            ]
        )
        npt.assert_array_equal(k2p(10, kappa=0.5).data.round(7), check)

        check = np.array(
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
            ]
        )
        npt.assert_array_equal(k2p(17, kappa=0.5).data.round(7), check)

        with self.assertRaises(ValueError):
            k2p(-1, kappa=0.5)

        with self.assertRaises(ValueError):
            k2p(1, kappa=0.0)

        with self.assertRaises(ValueError):
            k2p(1, kappa=-1)

        with self.assertRaises(ValueError):
            k2p(1, kappa=2)


class TestF81(TestCase):
    def test_instance(self):
        check = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        npt.assert_array_equal(f81(0, freqs=(0.25, 0.25, 0.25, 0.25)).data, check)
        npt.assert_allclose(
            f81(0, freqs=(0.25, 0.25, 0.25, 0.25)).data,
            f81(0, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )

        check = np.array(
            [
                [0.4464394, 0.1107121, 0.1660682, 0.2767803],
                [0.0, 0.5571515, 0.1660682, 0.2767803],
                [0.0, 0.1107121, 0.6125076, 0.2767803],
                [0.0, 0.1107121, 0.1660682, 0.7232197],
            ]
        )
        npt.assert_array_equal(
            f81(0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(7), check
        )

        check = np.array(
            [
                [0.6350628, 0.1216457, 0.1216457, 0.1216457],
                [0.1216457, 0.6350628, 0.1216457, 0.1216457],
                [0.1216457, 0.1216457, 0.6350628, 0.1216457],
                [0.1216457, 0.1216457, 0.1216457, 0.6350628],
            ]
        )
        npt.assert_array_equal(
            f81(0.5, freqs=(0.25, 0.25, 0.25, 0.25)).data.round(7), check
        )
        npt.assert_allclose(
            f81(0.5, freqs=(0.25, 0.25, 0.25, 0.25)).data,
            jc69(0.5).data,
        )

        check = np.array(
            [
                [0.1993081, 0.1601384, 0.2402076, 0.4003459],
                [0.0, 0.3594465, 0.2402076, 0.4003459],
                [0.0, 0.1601384, 0.4395157, 0.4003459],
                [0.0, 0.1601384, 0.2402076, 0.5996541],
            ]
        )
        npt.assert_array_equal(
            f81(1.0, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(7), check
        )

        check = np.array(
            [
                [0.0003145, 0.1999371, 0.2999056, 0.4998427],
                [0.0, 0.2002516, 0.2999056, 0.4998427],
                [0.0, 0.1999371, 0.3002202, 0.4998427],
                [0.0, 0.1999371, 0.2999056, 0.5001573],
            ]
        )
        npt.assert_array_equal(f81(5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(7), check)

        check = np.array(
            [
                [0.1245145, 0.3501942, 0.0, 0.5252913],
                [0.0, 0.4747087, 0.0, 0.5252913],
                [0.0, 0.3501942, 0.1245145, 0.5252913],
                [0.0, 0.3501942, 0.0, 0.6498058],
            ]
        )
        npt.assert_array_equal(f81(1, freqs=(0.0, 0.4, 0.0, 0.6)).data.round(7), check)

        check = np.array(
            [
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
            ]
        )
        npt.assert_array_equal(f81(11, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(7), check)

        with self.assertRaises(ValueError):
            f81(-1, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            f81(1, freqs=(0.5, 0.5, 0.5, 0.5))

        with self.assertRaises(ValueError):
            f81(1, freqs=(0.0, 0.2, 0.3, -0.5))


class TestHKY85(TestCase):
    def test_instance(self):
        check = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        npt.assert_array_equal(
            hky85(0, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data, check
        )
        npt.assert_array_equal(
            hky85(0, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            hky85(0, kappa=1, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )
        npt.assert_array_equal(
            hky85(0, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            hky85(0, kappa=0.5, freqs=(0.1, 0.4, 0.4, 0.1)).data,
        )

        check = np.array(
            [
                [0.5461625, 0.1107121, 0.0663451, 0.2767803],
                [0.0, 0.6283823, 0.1660682, 0.2055495],
                [0.0, 0.1107121, 0.6125076, 0.2767803],
                [0.0, 0.0822198, 0.1660682, 0.751712],
            ]
        )
        npt.assert_array_equal(
            hky85(0.5, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(7), check
        )
        check = np.array(
            [
                [0.6001238, 0.2124794, 0.1342769, 0.0531198],
                [0.0531198, 0.7008315, 0.2124794, 0.0335692],
                [0.0335692, 0.2124794, 0.7008315, 0.0531198],
                [0.0531198, 0.1342769, 0.2124794, 0.6001238],
            ]
        )
        npt.assert_array_equal(
            hky85(0.5, kappa=0.5, freqs=(0.1, 0.4, 0.4, 0.1)).data.round(7), check
        )

        # HKY85 should always produce results same results as simpler models for
        # certain parameter values.
        npt.assert_allclose(
            hky85(0.5, kappa=1.0, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            f81(0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )
        npt.assert_allclose(
            hky85(0.5, kappa=0.5, freqs=(0.25, 0.25, 0.25, 0.25)).data,
            k2p(0.5, kappa=0.5).data,
        )
        npt.assert_allclose(
            hky85(0.5, kappa=1, freqs=(0.25, 0.25, 0.25, 0.25)).data, jc69(0.5).data
        )

        check = np.array(
            [
                [0.2982935, 0.1601384, 0.1412222, 0.4003459],
                [0.0, 0.4301503, 0.2402076, 0.3296421],
                [0.0, 0.1601384, 0.4395157, 0.4003459],
                [0.0, 0.1318568, 0.2402076, 0.6279356],
            ]
        )
        npt.assert_array_equal(
            hky85(1, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(7), check
        )

        check = np.array(
            [
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
            ]
        )
        npt.assert_array_equal(
            hky85(14, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(7), check
        )

        with self.assertRaises(ValueError):
            hky85(-1, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            hky85(1, kappa=0.0, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            hky85(1, kappa=-1, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            hky85(1, kappa=2, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            hky85(1, kappa=0.5, freqs=(0.5, 0.5, 0.5, 0.5))

        with self.assertRaises(ValueError):
            hky85(1, kappa=0.5, freqs=(0.0, 0.2, 0.3, -0.5))


class TestTN93(TestCase):
    def test_instance(self):
        check = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        npt.assert_array_equal(
            tn93(0, kappa_r=0.4, kappa_y=0.6, freqs=(0.0, 0.2, 0.3, 0.5)).data, check
        )
        npt.assert_allclose(
            tn93(0, kappa_r=0.6, kappa_y=0.4, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            tn93(0, kappa_r=1, kappa_y=1, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )
        npt.assert_allclose(
            tn93(0, kappa_r=0.5, kappa_y=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            tn93(0, kappa_r=0.4, kappa_y=0.6, freqs=(0.1, 0.4, 0.4, 0.1)).data,
        )

        check = np.array(
            [
                [0.5686352, 0.1107121, 0.0438724, 0.2767803],
                [0.0, 0.6129648, 0.1660682, 0.2209671],
                [0.0, 0.1107121, 0.6125076, 0.2767803],
                [0.0, 0.0883868, 0.1660682, 0.745545],
            ]
        )
        npt.assert_array_equal(
            tn93(0.5, kappa_r=0.4, kappa_y=0.6, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(
                7
            ),
            check,
        )

        check = np.array(
            [
                [0.5245779, 0.1107121, 0.0879297, 0.2767803],
                [0.0, 0.6444343, 0.1660682, 0.1894976],
                [0.0, 0.1107121, 0.6125076, 0.2767803],
                [0.0, 0.075799, 0.1660682, 0.7581328],
            ]
        )
        npt.assert_array_equal(
            tn93(0.5, kappa_r=0.6, kappa_y=0.4, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(
                7
            ),
            check,
        )

        check = np.array(
            [
                [0.6176215, 0.2124794, 0.1167793, 0.0531198],
                [0.0531198, 0.6966197, 0.2124794, 0.037781],
                [0.0291948, 0.2124794, 0.7052059, 0.0531198],
                [0.0531198, 0.1511242, 0.2124794, 0.5832766],
            ]
        )
        npt.assert_array_equal(
            tn93(0.5, kappa_r=0.4, kappa_y=0.6, freqs=(0.1, 0.4, 0.4, 0.1)).data.round(
                7
            ),
            check,
        )

        # TN93 should always produce results same results as simpler models for
        # certain parameter values.
        npt.assert_allclose(
            tn93(0.5, kappa_r=0.5, kappa_y=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            hky85(0.5, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )
        npt.assert_allclose(
            tn93(0.5, kappa_r=1, kappa_y=1, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            f81(0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )
        npt.assert_allclose(
            tn93(0.5, kappa_r=0.5, kappa_y=0.5, freqs=(0.25, 0.25, 0.25, 0.25)).data,
            k2p(0.5, kappa=0.5).data,
        )
        npt.assert_allclose(
            tn93(0.5, kappa_r=1, kappa_y=1, freqs=(0.25, 0.25, 0.25, 0.25)).data,
            jc69(0.5).data,
        )

        check = np.array(
            [
                [0.323346, 0.1601384, 0.1161697, 0.4003459],
                [0.0, 0.4136421, 0.2402076, 0.3461503],
                [0.0, 0.1601384, 0.4395157, 0.4003459],
                [0.0, 0.1384601, 0.2402076, 0.6213323],
            ]
        )
        npt.assert_array_equal(
            tn93(1, kappa_r=0.4, kappa_y=0.6, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(7),
            check,
        )

        check = np.array(
            [
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
            ]
        )
        npt.assert_array_equal(
            tn93(15, kappa_r=0.4, kappa_y=0.6, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(
                7
            ),
            check,
        )

        with self.assertRaises(ValueError):
            tn93(-1, kappa_r=0.5, kappa_y=0.5, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            tn93(1, kappa_r=0.0, kappa_y=0.5, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            tn93(1, kappa_r=-1, kappa_y=0.5, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            tn93(1, kappa_r=2, kappa_y=0.5, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            tn93(1, kappa_r=0.5, kappa_y=0.0, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            tn93(1, kappa_r=0.5, kappa_y=-1, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            tn93(1, kappa_r=0.5, kappa_y=2, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            tn93(1, kappa_r=0.5, kappa_y=0.5, freqs=(0.5, 0.5, 0.5, 0.5))

        with self.assertRaises(ValueError):
            tn93(1, kappa_r=0.5, kappa_y=0.5, freqs=(0.0, 0.2, 0.3, -0.5))
