# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Array API backend tests for composition functions.

These tests verify that composition functions (closure, clr, clr_inv, ilr,
ilr_inv, alr, alr_inv) produce correct results across all supported array
backends (NumPy, JAX, PyTorch, CuPy) and devices (CPU, GPU).

By default only NumPy is tested. Set environment variables to test others::

    SKBIO_ARRAY_API_BACKEND=jax   python -m skbio.test
    SKBIO_ARRAY_API_BACKEND=all   python -m skbio.test
    SKBIO_ARRAY_API_BACKEND=torch SKBIO_DEVICE=cuda python -m skbio.test

"""

from unittest import TestCase, main

import numpy as np
from numpy.random import rand, randint

from skbio.stats.composition import (
    closure, clr, clr_inv, ilr, ilr_inv, alr, alr_inv,
)
from skbio.stats.composition._base import _gram_schmidt_basis
from skbio.util._testing import ArrayAPITestMixin, backends


# Tolerance for functions that lose precision on some backends (ilr, alr, etc.)
_RELAXED_RTOL = 1e-4
_STRICT_RTOL = 1e-7


class TestClosureArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for closure()."""

    def setUp(self):
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-4
        self.axis = 1

    @backends("numpy", "jax", "torch", "cupy")
    def test_closure_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        expected_np = self.data_int / np.sum(
            self.data_int, keepdims=True, axis=self.axis
        )
        result = closure(data, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assert_close(result, expected_np)

    @backends("numpy", "jax", "torch", "cupy")
    def test_closure_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        expected_np = self.data_real / np.sum(
            self.data_real, keepdims=True, axis=self.axis
        )
        result = closure(data, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assert_close(result, expected_np)


class TestCLRArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for clr()."""

    def setUp(self):
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-4
        self.axis = 1

        lmat = np.log(self.data_int)
        self.expected_int = lmat - np.mean(lmat, axis=self.axis, keepdims=True)

        lmat = np.log(self.data_real)
        self.expected_real = lmat - np.mean(lmat, axis=self.axis, keepdims=True)

    @backends("numpy", "jax", "torch", "cupy")
    def test_clr_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        result = clr(data, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(result.shape, data.shape)
        self.assert_close(result, self.expected_int, rtol=_RELAXED_RTOL)

    @backends("numpy", "jax", "torch", "cupy")
    def test_clr_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        result = clr(data, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(result.shape, data.shape)
        self.assert_close(result, self.expected_real, rtol=_RELAXED_RTOL)


class TestCLRInvArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for clr_inv()."""

    def setUp(self):
        self.axis = 1
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-4

        tmp = np.exp(
            self.data_int - np.max(self.data_int, axis=self.axis, keepdims=True)
        )
        self.expected_int = tmp / np.sum(tmp, axis=self.axis, keepdims=True)

        tmp = np.exp(
            self.data_real - np.max(self.data_real, axis=self.axis, keepdims=True)
        )
        self.expected_real = tmp / np.sum(tmp, axis=self.axis, keepdims=True)

    @backends("numpy", "jax", "torch", "cupy")
    def test_clr_inv_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        result = clr_inv(data, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(result.shape, data.shape)
        self.assert_close(result, self.expected_int, rtol=_RELAXED_RTOL)

    @backends("numpy", "jax", "torch", "cupy")
    def test_clr_inv_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        result = clr_inv(data, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(result.shape, data.shape)
        self.assert_close(result, self.expected_real, rtol=_RELAXED_RTOL)


class TestALRArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for alr()."""

    def setUp(self):
        self.ref_idx = 2
        self.axis = 1
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-6

        denom = self.data_int[:, [self.ref_idx], :]
        numerators = np.delete(self.data_int, self.ref_idx, axis=self.axis)
        self.expected_int = np.log(numerators / denom)

        denom = self.data_real[:, [self.ref_idx], :]
        numerators = np.delete(self.data_real, self.ref_idx, axis=self.axis)
        self.expected_real = np.log(numerators / denom)

        self.expected_shape_int = [
            d if i != self.axis else d - 1
            for i, d in enumerate(self.data_int.shape)
        ]
        self.expected_shape_real = [
            d if i != self.axis else d - 1
            for i, d in enumerate(self.data_real.shape)
        ]

    @backends("numpy", "jax", "torch", "cupy")
    def test_alr_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        result = alr(data, self.ref_idx, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_int)
        self.assert_close(result, self.expected_int, rtol=_RELAXED_RTOL)

    @backends("numpy", "jax", "torch", "cupy")
    def test_alr_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        result = alr(data, self.ref_idx, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_real)
        self.assert_close(result, self.expected_real, rtol=_RELAXED_RTOL)


class TestALRInvArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for alr_inv()."""

    def setUp(self):
        self.ref_idx = 2
        self.axis = 1
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-6

        exp_mat = np.exp(self.data_int)
        ones = np.ones_like(exp_mat[:, 0, :])
        exp_mat = np.insert(exp_mat, self.ref_idx, ones, axis=self.axis)
        self.expected_int = exp_mat / exp_mat.sum(axis=self.axis, keepdims=True)

        exp_mat = np.exp(self.data_real)
        ones = np.ones_like(exp_mat[:, 0, :])
        exp_mat = np.insert(exp_mat, self.ref_idx, ones, axis=self.axis)
        self.expected_real = exp_mat / exp_mat.sum(axis=self.axis, keepdims=True)

        self.expected_shape_int = [
            d if i != self.axis else d + 1
            for i, d in enumerate(self.data_int.shape)
        ]
        self.expected_shape_real = [
            d if i != self.axis else d + 1
            for i, d in enumerate(self.data_real.shape)
        ]

    @backends("numpy", "jax", "torch", "cupy")
    def test_alr_inv_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        result = alr_inv(data, self.ref_idx, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_int)
        self.assert_close(result, self.expected_int, rtol=_RELAXED_RTOL)

    @backends("numpy", "jax", "torch", "cupy")
    def test_alr_inv_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        result = alr_inv(data, self.ref_idx, self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_real)
        self.assert_close(result, self.expected_real, rtol=_RELAXED_RTOL)


class TestILRDefaultBasisArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for ilr() with default (Gram-Schmidt) basis."""

    def setUp(self):
        self.axis = 1
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-6

        V = _gram_schmidt_basis(6)
        x = clr(self.data_int, self.axis)
        x = np.moveaxis(x, self.axis, -1)
        self.expected_int = np.moveaxis(x @ V.T, -1, self.axis)

        x = clr(self.data_real, self.axis)
        x = np.moveaxis(x, self.axis, -1)
        self.expected_real = np.moveaxis(x @ V.T, -1, self.axis)

        self.expected_shape_int = [
            d if i != self.axis else d - 1
            for i, d in enumerate(self.data_int.shape)
        ]
        self.expected_shape_real = [
            d if i != self.axis else d - 1
            for i, d in enumerate(self.data_real.shape)
        ]

    @backends("numpy", "jax", "torch", "cupy")
    def test_ilr_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        result = ilr(data, axis=self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_int)
        self.assert_close(result, self.expected_int, rtol=_RELAXED_RTOL)

    @backends("numpy", "jax", "torch", "cupy")
    def test_ilr_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        result = ilr(data, axis=self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_real)
        self.assert_close(result, self.expected_real, rtol=_RELAXED_RTOL)


class TestILRInvDefaultBasisArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for ilr_inv() with default (Gram-Schmidt) basis."""

    def setUp(self):
        self.axis = 1
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-6

        V = _gram_schmidt_basis(7)
        y = np.moveaxis(self.data_int, self.axis, -1)
        y = np.exp(y @ V)
        y = y / np.sum(y, axis=-1, keepdims=True)
        self.expected_int = np.moveaxis(y, -1, self.axis)

        y = np.moveaxis(self.data_real, self.axis, -1)
        y = np.exp(y @ V)
        y = y / np.sum(y, axis=-1, keepdims=True)
        self.expected_real = np.moveaxis(y, -1, self.axis)

        self.expected_shape_int = [
            d if i != self.axis else d + 1
            for i, d in enumerate(self.data_int.shape)
        ]
        self.expected_shape_real = [
            d if i != self.axis else d + 1
            for i, d in enumerate(self.data_real.shape)
        ]

    @backends("numpy", "jax", "torch", "cupy")
    def test_ilr_inv_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        result = ilr_inv(data, axis=self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_int)
        self.assert_close(result, self.expected_int, rtol=_RELAXED_RTOL)

    @backends("numpy", "jax", "torch", "cupy")
    def test_ilr_inv_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        result = ilr_inv(data, axis=self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_real)
        self.assert_close(result, self.expected_real, rtol=_RELAXED_RTOL)


class TestILRHelmertBasisArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for ilr() with Helmert sub-composition basis."""

    def setUp(self):
        self.axis = 1
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-6

        D = 6
        V = np.zeros((D, D - 1))
        for i in range(1, D):
            V[:i, i - 1] = 1 / i
            V[i, i - 1] = -1
        V = V / np.linalg.norm(V, axis=0)
        V = V.T  # (n_basis, n_component)

        x = clr(self.data_int, self.axis)
        x = np.moveaxis(x, self.axis, -1)
        self.expected_int = np.moveaxis(x @ V.T, -1, self.axis)

        x = clr(self.data_real, self.axis)
        x = np.moveaxis(x, self.axis, -1)
        self.expected_real = np.moveaxis(x @ V.T, -1, self.axis)

        self.expected_shape_int = [
            d if i != self.axis else d - 1
            for i, d in enumerate(self.data_int.shape)
        ]
        self.expected_shape_real = [
            d if i != self.axis else d - 1
            for i, d in enumerate(self.data_real.shape)
        ]

    @backends("numpy", "jax", "torch", "cupy")
    def test_ilr_helmert_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        result = ilr(data, axis=self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_int)
        self.assert_close(result, self.expected_int, rtol=_RELAXED_RTOL)

    @backends("numpy", "jax", "torch", "cupy")
    def test_ilr_helmert_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        result = ilr(data, axis=self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_real)
        self.assert_close(result, self.expected_real, rtol=_RELAXED_RTOL)


class TestILRInvHelmertBasisArrayAPI(TestCase, ArrayAPITestMixin):
    """Array API tests for ilr_inv() with Helmert sub-composition basis."""

    def setUp(self):
        self.axis = 1
        self.data_int = randint(1, 100, (4, 6, 3))
        self.data_real = rand(4, 6, 3) + 1e-6

        D = 7
        V = np.zeros((D, D - 1))
        for i in range(1, D):
            V[:i, i - 1] = 1 / i
            V[i, i - 1] = -1
        V = V / np.linalg.norm(V, axis=0)
        V = V.T  # (n_basis, n_component)

        y = np.moveaxis(self.data_int, self.axis, -1)
        y = np.exp(y @ V)
        y = y / np.sum(y, axis=-1, keepdims=True)
        self.expected_int = np.moveaxis(y, -1, self.axis)

        y = np.moveaxis(self.data_real, self.axis, -1)
        y = np.exp(y @ V)
        y = y / np.sum(y, axis=-1, keepdims=True)
        self.expected_real = np.moveaxis(y, -1, self.axis)

        self.expected_shape_int = [
            d if i != self.axis else d + 1
            for i, d in enumerate(self.data_int.shape)
        ]
        self.expected_shape_real = [
            d if i != self.axis else d + 1
            for i, d in enumerate(self.data_real.shape)
        ]

    @backends("numpy", "jax", "torch", "cupy")
    def test_ilr_inv_helmert_int(self, xp, device):
        data = self.make_array(xp, device, self.data_int)
        result = ilr_inv(data, axis=self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_int)
        self.assert_close(result, self.expected_int, rtol=_RELAXED_RTOL)

    @backends("numpy", "jax", "torch", "cupy")
    def test_ilr_inv_helmert_real(self, xp, device):
        data = self.make_array(xp, device, self.data_real)
        result = ilr_inv(data, axis=self.axis)

        self.assert_type_preserved(result, xp, device)
        self.assertEqual(list(result.shape), self.expected_shape_real)
        self.assert_close(result, self.expected_real, rtol=_RELAXED_RTOL)


if __name__ == "__main__":
    main()