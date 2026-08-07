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

from skbio.sequence.transition import jc69, k2p, f81, hky85, tn93


class TestJC69(TestCase):
    def test_instance(self):
        obs = jc69(0).data
        exp = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = jc69(0.5).data.round(5)
        exp = np.array(
            [
                [0.63506, 0.12165, 0.12165, 0.12165],
                [0.12165, 0.63506, 0.12165, 0.12165],
                [0.12165, 0.12165, 0.63506, 0.12165],
                [0.12165, 0.12165, 0.12165, 0.63506],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = jc69(1).data.round(5)
        exp = np.array(
            [
                [0.44770, 0.18410, 0.18410, 0.18410],
                [0.18410, 0.44770, 0.18410, 0.18410],
                [0.18410, 0.18410, 0.44770, 0.18410],
                [0.18410, 0.18410, 0.18410, 0.44770],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = jc69(10).data.round(5)
        exp = np.array(
            [
                [0.25000, 0.25000, 0.25000, 0.25000],
                [0.25000, 0.25000, 0.25000, 0.25000],
                [0.25000, 0.25000, 0.25000, 0.25000],
                [0.25000, 0.25000, 0.25000, 0.25000],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = jc69(14).data.round(5)
        exp = np.array(
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        with self.assertRaises(ValueError):
            jc69(-1)

        with self.assertRaises(TypeError):
            jc69(1, seqtype="Protein")


class TestK2P(TestCase):
    def test_instance(self):
        obs = k2p(0, kappa=0.5).data
        exp = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)
        npt.assert_array_equal(k2p(0, kappa=0.5).data, k2p(0, kappa=1).data)

        obs = k2p(0.5, kappa=0.5).data.round(5)
        exp = np.array(
            [
                [0.68162, 0.12165, 0.07509, 0.12165],
                [0.12165, 0.68162, 0.12165, 0.07509],
                [0.07509, 0.12165, 0.68162, 0.12165],
                [0.12165, 0.07509, 0.12165, 0.68162],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = k2p(0.5, kappa=1).data.round(5)
        exp = np.array(
            [
                [0.63506, 0.12165, 0.12165, 0.12165],
                [0.12165, 0.63506, 0.12165, 0.12165],
                [0.12165, 0.12165, 0.63506, 0.12165],
                [0.12165, 0.12165, 0.12165, 0.63506],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)
        npt.assert_allclose(k2p(0.5, kappa=1).data, jc69(0.5).data)

        obs = k2p(1.0, kappa=0.5).data.round(5)
        exp = np.array(
            [
                [0.49984, 0.18410, 0.13196, 0.18410],
                [0.18410, 0.49984, 0.18410, 0.13196],
                [0.13196, 0.18410, 0.49984, 0.18410],
                [0.18410, 0.13196, 0.18410, 0.49984],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = k2p(10, kappa=0.5).data.round(5)
        exp = np.array(
            [
                [0.25002, 0.25000, 0.24998, 0.25000],
                [0.25000, 0.25002, 0.25000, 0.24998],
                [0.24998, 0.25000, 0.25002, 0.25000],
                [0.25000, 0.24998, 0.25000, 0.25002],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = k2p(17, kappa=0.5).data.round(5)
        exp = np.array(
            [
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
                [0.25, 0.25, 0.25, 0.25],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

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
        obs = f81(0, freqs=(0.25, 0.25, 0.25, 0.25)).data
        exp = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)
        npt.assert_allclose(
            f81(0, freqs=(0.25, 0.25, 0.25, 0.25)).data,
            f81(0, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )

        obs = f81(0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(5)
        exp = np.array(
            [
                [0.44644, 0.11071, 0.16607, 0.27678],
                [0.0, 0.55715, 0.16607, 0.27678],
                [0.0, 0.11071, 0.61251, 0.27678],
                [0.0, 0.11071, 0.16607, 0.72322],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = f81(0.5, freqs=(0.25, 0.25, 0.25, 0.25)).data.round(5)
        exp = np.array(
            [
                [0.63506, 0.12165, 0.12165, 0.12165],
                [0.12165, 0.63506, 0.12165, 0.12165],
                [0.12165, 0.12165, 0.63506, 0.12165],
                [0.12165, 0.12165, 0.12165, 0.63506],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)
        npt.assert_allclose(
            f81(0.5, freqs=(0.25, 0.25, 0.25, 0.25)).data,
            jc69(0.5).data,
        )

        obs = f81(1.0, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(5)
        exp = np.array(
            [
                [0.19931, 0.16014, 0.24021, 0.40035],
                [0.0, 0.35945, 0.24021, 0.40035],
                [0.0, 0.16014, 0.43952, 0.40035],
                [0.0, 0.16014, 0.24021, 0.59965],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = f81(5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(5)
        exp = np.array(
            [
                [0.00031, 0.19994, 0.29991, 0.49984],
                [0.0, 0.20025, 0.29991, 0.49984],
                [0.0, 0.19994, 0.30022, 0.49984],
                [0.0, 0.19994, 0.29991, 0.50016],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = f81(1, freqs=(0.0, 0.4, 0.0, 0.6)).data.round(5)
        exp = np.array(
            [
                [0.12451, 0.35019, 0.0, 0.52529],
                [0.0, 0.47471, 0.0, 0.52529],
                [0.0, 0.35019, 0.12451, 0.52529],
                [0.0, 0.35019, 0.0, 0.64981],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = f81(11, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(5)
        exp = np.array(
            [
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        with self.assertRaises(ValueError):
            f81(-1, freqs=(0.0, 0.2, 0.3, 0.5))

        with self.assertRaises(ValueError):
            f81(1, freqs=(0.5, 0.5, 0.5, 0.5))

        with self.assertRaises(ValueError):
            f81(1, freqs=(0.0, 0.2, 0.3, -0.5))


class TestHKY85(TestCase):
    def test_instance(self):
        obs = hky85(0, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data
        exp = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)
        npt.assert_array_equal(
            hky85(0, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            hky85(0, kappa=1, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )
        npt.assert_array_equal(
            hky85(0, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            hky85(0, kappa=0.5, freqs=(0.1, 0.4, 0.4, 0.1)).data,
        )

        obs = hky85(0.5, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(5)
        exp = np.array(
            [
                [0.54616, 0.11071, 0.06635, 0.27678],
                [0.0, 0.62838, 0.16607, 0.20555],
                [0.0, 0.11071, 0.61251, 0.27678],
                [0.0, 0.08222, 0.16607, 0.75171],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)
        obs = hky85(0.5, kappa=0.5, freqs=(0.1, 0.4, 0.4, 0.1)).data.round(5)
        exp = np.array(
            [
                [0.60012, 0.21248, 0.13428, 0.05312],
                [0.05312, 0.70083, 0.21248, 0.03357],
                [0.03357, 0.21248, 0.70083, 0.05312],
                [0.05312, 0.13428, 0.21248, 0.60012],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

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

        obs = hky85(1, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(5)
        exp = np.array(
            [
                [0.29829, 0.16014, 0.14122, 0.40035],
                [0.0, 0.43015, 0.24021, 0.32964],
                [0.0, 0.16014, 0.43952, 0.40035],
                [0.0, 0.13186, 0.24021, 0.62794],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = hky85(14, kappa=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(5)
        exp = np.array(
            [
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

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
        obs = tn93(0, kappa_r=0.4, kappa_y=0.6, freqs=(0.0, 0.2, 0.3, 0.5)).data
        exp = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)
        npt.assert_allclose(
            tn93(0, kappa_r=0.6, kappa_y=0.4, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            tn93(0, kappa_r=1, kappa_y=1, freqs=(0.0, 0.2, 0.3, 0.5)).data,
        )
        npt.assert_allclose(
            tn93(0, kappa_r=0.5, kappa_y=0.5, freqs=(0.0, 0.2, 0.3, 0.5)).data,
            tn93(0, kappa_r=0.4, kappa_y=0.6, freqs=(0.1, 0.4, 0.4, 0.1)).data,
        )

        obs = tn93(0.5, kappa_r=0.4, kappa_y=0.6, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(
                5
            )
        exp = np.array(
            [
                [0.56864, 0.11071, 0.04387, 0.27678],
                [0.0, 0.61296, 0.16607, 0.22097],
                [0.0, 0.11071, 0.61251, 0.27678],
                [0.0, 0.08839, 0.16607, 0.74554],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = tn93(0.5, kappa_r=0.6, kappa_y=0.4, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(
                5
            )
        exp = np.array(
            [
                [0.52458, 0.11071, 0.08793, 0.27678],
                [0.0, 0.64443, 0.16607, 0.18950],
                [0.0, 0.11071, 0.61251, 0.27678],
                [0.0, 0.07580, 0.16607, 0.75813],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = tn93(0.5, kappa_r=0.4, kappa_y=0.6, freqs=(0.1, 0.4, 0.4, 0.1)).data.round(
                5
            )
        exp = np.array(
            [
                [0.61762, 0.21248, 0.11678, 0.05312],
                [0.05312, 0.69662, 0.21248, 0.03778],
                [0.02919, 0.21248, 0.70521, 0.05312],
                [0.05312, 0.15112, 0.21248, 0.58328],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

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

        obs = tn93(1, kappa_r=0.4, kappa_y=0.6, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(5)
        exp = np.array(
            [
                [0.32335, 0.16014, 0.11617, 0.40035],
                [0.0, 0.41364, 0.24021, 0.34615],
                [0.0, 0.16014, 0.43952, 0.40035],
                [0.0, 0.13846, 0.24021, 0.62133],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

        obs = tn93(15, kappa_r=0.4, kappa_y=0.6, freqs=(0.0, 0.2, 0.3, 0.5)).data.round(
                5
            )
        exp = np.array(
            [
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
                [0.0, 0.2, 0.3, 0.5],
            ], dtype=np.float32
        )
        npt.assert_array_equal(obs, exp)

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


if __name__ == "__main__":
    main()
