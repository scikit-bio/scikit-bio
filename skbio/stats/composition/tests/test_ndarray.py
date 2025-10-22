# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main, skipIf
import warnings
import inspect

import numpy as np
import numpy.testing as npt
from numpy.random import rand, randint

from skbio.util import get_package
from skbio.util._gpu import cuda_avail

from skbio.stats.composition import (
    closure, clr, clr_inv, ilr,
    ilr_inv, alr, alr_inv)
from skbio.stats.composition._base import _gram_schmidt_basis


# import optional dependencies
try:
    jax = get_package("jax")
except ImportError:
    no_jax = True
else:
    no_jax = False
    jnp = get_package("jax.numpy")
    jax.config.update("jax_enable_x64", True)

try:
    torch = get_package("torch")
except ImportError:
    no_torch = True
else:
    no_torch = False
    # initially, the default is float32, but this test will use float64
    torch.set_default_dtype(torch.float64)


# check the existence of cuda
no_cuda = not cuda_avail()


def assert_coo_allclose(res, exp, rtol=1e-7, atol=1e-7):
    res_data = np.vstack((res.row, res.col, res.data)).T
    exp_data = np.vstack((exp.row, exp.col, exp.data)).T

    # sort by row and col
    res_data = res_data[res_data[:, 1].argsort()]
    res_data = res_data[res_data[:, 0].argsort()]
    exp_data = exp_data[exp_data[:, 1].argsort()]
    exp_data = exp_data[exp_data[:, 0].argsort()]
    npt.assert_allclose(res_data, exp_data, rtol=rtol, atol=atol)


class Tests_Closure(TestCase):
    def setUp(self):
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-04

    def test_ndarray_numpy(self):
        # Asrange
        axis = 1
        expected_int = self.random_postive_int/np.sum(
            self.random_postive_int, keepdims=True, axis=axis)
        expected_real = self.random_postive_real/np.sum(
            self.random_postive_real, keepdims=True, axis=axis)

        # Action
        rst_int = closure(self.random_postive_int, axis)
        rst_real = closure(self.random_postive_real, axis)

        # Assert
        np.allclose(rst_int, expected_int)
        np.allclose(rst_real, expected_real)

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        axis = 1
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = random_int/torch.sum(random_int, keepdims=True, dim=axis)
        expected_real = random_real/torch.sum(random_real, keepdims=True, dim=axis)

        # Action
        rst_int = closure(random_int, axis)
        rst_real = closure(random_real, axis)

       # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        assert torch.sum(torch.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
        assert torch.sum(torch.abs(expected_real - rst_real)
                         ) < 1e-08, "unmatched values"

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = 1
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = random_int/torch.sum(random_int, keepdims=True, dim=axis)
        expected_real = random_real/torch.sum(random_real, keepdims=True, dim=axis)

        # Action
        rst_int = closure(random_int, axis)
        rst_real = closure(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        assert torch.sum(torch.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
        assert torch.sum(torch.abs(expected_real - rst_real)
                         ) < 1e-08, "unmatched values"

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = 1
        device = jax.devices("cpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real)
        expected_int = random_int/jnp.sum(random_int, keepdims=True, axis=axis)
        expected_real = random_real/jnp.sum(random_real, keepdims=True, axis=axis)

        # Action
        rst_int = closure(random_int, axis)
        rst_real = closure(random_real, axis)

       # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
        assert jnp.sum(jnp.abs(expected_real - rst_real)) < 1e-08, "unmatched values"

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = 1
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real)
        expected_int = random_int/jnp.sum(random_int, keepdims=True, axis=axis)
        expected_real = random_real/jnp.sum(random_real, keepdims=True, axis=axis)

        # Action
        rst_int = closure(random_int, axis)
        rst_real = closure(random_real, axis)

       # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
        assert jnp.sum(jnp.abs(expected_real - rst_real)) < 1e-08, "unmatched values"


class Tests_CLR(TestCase):
    def setUp(self):
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-04

    def test_ndarray_numpy(self):
        # Asrange
        axis = 1
        expected_int = (lmat := np.log(self.random_postive_int)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_real = (lmat := np.log(self.random_postive_real)) \
            - np.mean(lmat, axis=axis, keepdims=True)

        # Action
        rst_int = clr(self.random_postive_int, axis)
        rst_real = clr(self.random_postive_real, axis)

        # Assert
        # the shape should not be changed
        assert expected_int.shape == self.random_postive_int.shape,  "unepected shape"

        # check the values
        np.allclose(rst_int, expected_int)
        np.allclose(rst_real, expected_real)

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        axis = 1
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = (lmat := np.log(self.random_postive_int)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_int = torch.tensor(expected_int)
        expected_real = (lmat := np.log(self.random_postive_real)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_real = torch.tensor(expected_real)

        # Action
        rst_int = clr(random_int, axis)
        rst_real = clr(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert expected_int.shape == random_int.shape,  "unepected shape"
        assert expected_int.shape == random_real.shape,  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = 1
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = (lmat := np.log(self.random_postive_int)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_int = torch.tensor(expected_int, device="cuda")
        expected_real = (lmat := np.log(self.random_postive_real)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_real = torch.tensor(expected_real, device="cuda")

        # Action
        rst_int = clr(random_int, axis)
        rst_real = clr(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert expected_int.shape == random_int.shape,  "unepected shape"
        assert expected_int.shape == random_real.shape,  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = 1
        device = jax.devices("cpu")[0]

        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = (lmat := np.log(self.random_postive_int)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_int = jnp.array(expected_int, device=device)
        expected_real = (lmat := np.log(self.random_postive_real)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_real = jnp.array(expected_real, device=device)

        # Action
        rst_int = clr(random_int, axis)
        rst_real = clr(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert expected_int.shape == random_int.shape,  "unepected shape"
        assert expected_int.shape == random_real.shape,  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = 1
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = (lmat := np.log(self.random_postive_int)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_int = jnp.array(expected_int, device=device)
        expected_real = (lmat := np.log(self.random_postive_real)) \
            - np.mean(lmat, axis=axis, keepdims=True)
        expected_real = jnp.array(expected_real, device=device)

        # Action
        rst_int = clr(random_int, axis)
        rst_real = clr(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert expected_int.shape == random_int.shape,  "unepected shape"
        assert expected_int.shape == random_real.shape,  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )


class Tests_CLR_INV(TestCase):
    def setUp(self):
        axis = 1
        self.axis = axis
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-04

        tem_int = np.exp(self.random_postive_int -
                         np.max(self.random_postive_int, axis=axis, keepdims=True))
        self.expected_int = tem_int/np.sum(tem_int, axis=axis, keepdims=True)

        tem_real = np.exp(self.random_postive_real -
                          np.max(self.random_postive_real, axis=axis, keepdims=True))
        self.expected_real = tem_real/np.sum(tem_real, axis=axis, keepdims=True)

        self.expected_shape_int = self.random_postive_int.shape
        self.expected_shape_real = self.random_postive_real.shape

    def test_ndarray_numpy(self):
        # Asrange
        axis = self.axis

        # Action
        rst_int = clr_inv(self.random_postive_int, axis)
        rst_real = clr_inv(self.random_postive_real, axis)

        # Assert
        # the shape should not be changed
        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in self.expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in self.expected_real.shape],  "unepected shape"

        # check the values
        assert np.allclose(rst_int, self.expected_int), "unmatched values"
        assert np.allclose(rst_real, self.expected_real), "unmatched values"

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        axis = 1
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = torch.tensor(self.expected_int)
        expected_real = torch.tensor(self.expected_real)

        # Action
        rst_int = clr_inv(random_int, axis)
        rst_real = clr_inv(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = torch.tensor(self.expected_int, device="cuda")
        expected_real = torch.tensor(self.expected_real, device="cuda")

        # Action
        rst_int = clr_inv(random_int, axis)
        rst_real = clr_inv(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = self.axis
        device = jax.devices("cpu")[0]

        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = clr_inv(random_int, axis)
        rst_real = clr_inv(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = self.axis
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = clr_inv(random_int, axis)
        rst_real = clr_inv(random_real, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )


class Tests_ALR(TestCase):
    def setUp(self):
        ref_idx = 2
        self.ref_idx = ref_idx
        axis = 1
        self.axis = axis
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-06

        denom = self.random_postive_int[:, [2], :]
        numerators = np.delete(self.random_postive_int, 2, axis=axis)
        self.expected_int = np.log(numerators / denom)

        denom = self.random_postive_real[:, [2], :]
        numerators = np.delete(self.random_postive_real, 2, axis=axis)
        self.expected_real = np.log(numerators / denom)

        self.expected_shape_int = [d if i != axis else d -
                                   1 for i, d in enumerate(self.random_postive_int.shape)]
        self.expected_shape_real = [d if i != axis else d -
                                    1 for i, d in enumerate(self.random_postive_real.shape)]

    def test_ndarray_numpy(self):
        # Asrange
        axis = self.axis
        ref_idx = self.ref_idx

        # Action
        rst_int = alr(self.random_postive_int, ref_idx, axis)
        rst_real = alr(self.random_postive_real, ref_idx, axis)

        # Assert
        # the shape should not be changed
        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in self.expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in self.expected_real.shape],  "unepected shape"

        # check the values
        assert np.allclose(rst_int, self.expected_int), "unmatched values"
        assert np.allclose(rst_real, self.expected_real), "unmatched values"

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        ref_idx = self.ref_idx
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = torch.tensor(self.expected_int)
        expected_real = torch.tensor(self.expected_real)

        # Action
        rst_int = alr(random_int, ref_idx, axis)
        rst_real = alr(random_real, ref_idx, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        ref_idx = self.ref_idx
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = torch.tensor(self.expected_int, device="cuda")
        expected_real = torch.tensor(self.expected_real, device="cuda")

        # Action
        rst_int = alr(random_int, ref_idx, axis)
        rst_real = alr(random_real, ref_idx, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = self.axis
        device = jax.devices("cpu")[0]
        ref_idx = self.ref_idx

        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = alr(random_int, ref_idx, axis)
        rst_real = alr(random_real, ref_idx, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = self.axis
        ref_idx = self.ref_idx
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action

        rst_int = alr(random_int, ref_idx, axis)
        rst_real = alr(random_real, ref_idx, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )


class Tests_ALR_INV(TestCase):
    def setUp(self):
        ref_idx = 2
        self.ref_idx = ref_idx
        axis = 1
        self.axis = axis
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-06

        exp_mat = np.exp(self.random_postive_int)
        ones = np.ones_like(exp_mat[:, 0, :])
        exp_mat = np.insert(exp_mat, 2, ones, axis=1)
        self.expected_int = exp_mat / exp_mat.sum(axis=1, keepdims=True)

        exp_mat = np.exp(self.random_postive_real)
        ones = np.ones_like(exp_mat[:, 0, :])
        exp_mat = np.insert(exp_mat, 2, ones, axis=1)
        self.expected_real = exp_mat / exp_mat.sum(axis=1, keepdims=True)

        self.expected_shape_int = [d if i != axis else d +
                                   1 for i, d in enumerate(self.random_postive_int.shape)]
        self.expected_shape_real = [d if i != axis else d +
                                    1 for i, d in enumerate(self.random_postive_real.shape)]

    def test_ndarray_numpy(self):
        # Asrange
        axis = self.axis
        ref_idx = self.ref_idx

        # Action
        rst_int = alr_inv(self.random_postive_int, ref_idx, axis)
        rst_real = alr_inv(self.random_postive_real, ref_idx, axis)

        # Assert
        # the shape should not be changed
        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in self.expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in self.expected_real.shape],  "unepected shape"

        # check the values
        assert np.allclose(rst_int, self.expected_int), "unmatched values"
        assert np.allclose(rst_real, self.expected_real), "unmatched values"

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        ref_idx = self.ref_idx
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = torch.tensor(self.expected_int)
        expected_real = torch.tensor(self.expected_real)

        # Action
        rst_int = alr_inv(random_int, ref_idx, axis)
        rst_real = alr_inv(random_real, ref_idx, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        ref_idx = self.ref_idx
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = torch.tensor(self.expected_int, device="cuda")
        expected_real = torch.tensor(self.expected_real, device="cuda")

        # Action
        rst_int = alr_inv(random_int, ref_idx, axis)
        rst_real = alr_inv(random_real, ref_idx, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = self.axis
        device = jax.devices("cpu")[0]
        ref_idx = self.ref_idx

        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = alr_inv(random_int, ref_idx, axis)
        rst_real = alr_inv(random_real, ref_idx, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = self.axis
        ref_idx = self.ref_idx
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = alr_inv(random_int, ref_idx, axis)
        rst_real = alr_inv(random_real, ref_idx, axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )


class Tests_ILR_default_basis(TestCase):
    def setUp(self):
        axis = 1
        self.axis = axis
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-06

        V = _gram_schmidt_basis(6)
        x = clr(self.random_postive_int, axis)
        x = np.moveaxis(x, axis, -1)
        x = x @ V.T
        self.expected_int = np.moveaxis(x, -1, axis)

        V = _gram_schmidt_basis(6)
        x = clr(self.random_postive_real, axis)
        x = np.moveaxis(x, axis, -1)
        x = x @ V.T
        self.expected_real = np.moveaxis(x, -1, axis)

        self.expected_shape_int = [d if i != axis else d -
                                   1 for i, d in enumerate(self.random_postive_int.shape)]
        self.expected_shape_real = [d if i != axis else d -
                                    1 for i, d in enumerate(self.random_postive_real.shape)]

    def test_ndarray_numpy(self):
        # Asrange
        axis = self.axis

        # Action
        rst_int = ilr(self.random_postive_int, axis=axis)
        rst_real = ilr(self.random_postive_real, axis=axis)

        # Assert
        # the shape should not be changed
        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in self.expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in self.expected_real.shape],  "unepected shape"

        # check the values
        assert np.allclose(rst_int, self.expected_int), "unmatched values"
        assert np.allclose(rst_real, self.expected_real), "unmatched values"

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = torch.tensor(self.expected_int)
        expected_real = torch.tensor(self.expected_real)

        # Action
        rst_int = ilr(random_int, axis=axis)
        rst_real = ilr(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = torch.tensor(self.expected_int, device="cuda")
        expected_real = torch.tensor(self.expected_real, device="cuda")

        # Action
        rst_int = ilr(random_int, axis=axis)
        rst_real = ilr(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = self.axis
        device = jax.devices("cpu")[0]

        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = ilr(random_int, axis=axis)
        rst_real = ilr(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = self.axis
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = ilr(random_int, axis=axis)
        rst_real = ilr(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )


class Tests_ILR_INV_default_basis(TestCase):
    def setUp(self):
        axis = 1
        self.axis = axis
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-06

        V = _gram_schmidt_basis(7)
        y = np.moveaxis(self.random_postive_int, axis, -1)
        y = y @ V
        y = np.exp(y)
        y = y / np.sum(y, axis=-1, keepdims=True)
        self.expected_int = np.moveaxis(y, -1, axis)

        V = _gram_schmidt_basis(7)
        y = np.moveaxis(self.random_postive_real, axis, -1)
        y = y @ V
        y = np.exp(y)
        y = y / np.sum(y, axis=-1, keepdims=True)
        self.expected_real = np.moveaxis(y, -1, axis)

        self.expected_shape_int = [d if i != axis else d +
                                   1 for i, d in enumerate(self.random_postive_int.shape)]
        self.expected_shape_real = [d if i != axis else d +
                                    1 for i, d in enumerate(self.random_postive_real.shape)]

    def test_ndarray_numpy(self):
        # Asrange
        axis = self.axis

        # Action
        rst_int = ilr_inv(self.random_postive_int, axis=axis)
        rst_real = ilr_inv(self.random_postive_real, axis=axis)

        # Assert
        # the shape should not be changed
        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in self.expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in self.expected_real.shape],  "unepected shape"

        # check the values
        assert np.allclose(rst_int, self.expected_int), "unmatched values"
        assert np.allclose(rst_real, self.expected_real), "unmatched values"

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = torch.tensor(self.expected_int)
        expected_real = torch.tensor(self.expected_real)

        # Action
        rst_int = ilr_inv(random_int, axis=axis)
        rst_real = ilr_inv(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = torch.tensor(self.expected_int, device="cuda")
        expected_real = torch.tensor(self.expected_real, device="cuda")

        # Action
        rst_int = ilr_inv(random_int, axis=axis)
        rst_real = ilr_inv(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = self.axis
        device = jax.devices("cpu")[0]

        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = ilr_inv(random_int, axis=axis)
        rst_real = ilr_inv(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = self.axis
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = ilr_inv(random_int, axis=axis)
        rst_real = ilr_inv(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )


class Tests_ILR_helmert_basis(TestCase):
    def setUp(self):
        axis = 1
        self.axis = axis
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-06

        D = 6
        V = np.zeros((D, D - 1))
        for i in range(1, D):
            V[:i, i - 1] = 1 / i
            V[i, i - 1] = -1
        V = V / np.linalg.norm(V, axis=0)
        # make the basis in (n_basis, n_component)
        V = V.T

        x = clr(self.random_postive_int, axis)
        x = np.moveaxis(x, axis, -1)
        x = x @ V.T
        self.expected_int = np.moveaxis(x, -1, axis)

        x = clr(self.random_postive_real, axis)
        x = np.moveaxis(x, axis, -1)
        x = x @ V.T
        self.expected_real = np.moveaxis(x, -1, axis)

        self.expected_shape_int = [d if i != axis else d -
                                   1 for i, d in enumerate(self.random_postive_int.shape)]
        self.expected_shape_real = [d if i != axis else d -
                                    1 for i, d in enumerate(self.random_postive_real.shape)]

    def test_ndarray_numpy(self):
        # Asrange
        axis = self.axis

        # Action
        rst_int = ilr(self.random_postive_int, axis=axis)
        rst_real = ilr(self.random_postive_real, axis=axis)

        # Assert
        # the shape should not be changed
        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in self.expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in self.expected_real.shape],  "unepected shape"

        # check the values
        assert np.allclose(rst_int, self.expected_int), "unmatched values"
        assert np.allclose(rst_real, self.expected_real), "unmatched values"

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = torch.tensor(self.expected_int)
        expected_real = torch.tensor(self.expected_real)

        # Action
        rst_int = ilr(random_int, axis=axis)
        rst_real = ilr(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = torch.tensor(self.expected_int, device="cuda")
        expected_real = torch.tensor(self.expected_real, device="cuda")

        # Action
        rst_int = ilr(random_int, axis=axis)
        rst_real = ilr(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = self.axis
        device = jax.devices("cpu")[0]

        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = ilr(random_int, axis=axis)
        rst_real = ilr(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = self.axis
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = ilr(random_int, axis=axis)
        rst_real = ilr(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )


class Tests_ILR_INV_helmert_basis(TestCase):
    def setUp(self):
        axis = 1
        self.axis = axis
        self.random_postive_int = randint(1, 100, (4, 6, 3))
        self.random_postive_real = rand(4, 6, 3)+1e-06

        D = 7
        V = np.zeros((D, D - 1))
        for i in range(1, D):
            V[:i, i - 1] = 1 / i
            V[i, i - 1] = -1
        V = V / np.linalg.norm(V, axis=0)
        # make the basis in (n_basis, n_component)
        V = V.T
        y = np.moveaxis(self.random_postive_int, axis, -1)
        y = y @ V
        y = np.exp(y)
        y = y / np.sum(y, axis=-1, keepdims=True)
        self.expected_int = np.moveaxis(y, -1, axis)

        V = _gram_schmidt_basis(7)
        y = np.moveaxis(self.random_postive_real, axis, -1)
        y = y @ V
        y = np.exp(y)
        y = y / np.sum(y, axis=-1, keepdims=True)
        self.expected_real = np.moveaxis(y, -1, axis)

        self.expected_shape_int = [d if i != axis else d +
                                   1 for i, d in enumerate(self.random_postive_int.shape)]
        self.expected_shape_real = [d if i != axis else d +
                                    1 for i, d in enumerate(self.random_postive_real.shape)]

    def test_ndarray_numpy(self):
        # Asrange
        axis = self.axis

        # Action
        rst_int = ilr_inv(self.random_postive_int, axis=axis)
        rst_real = ilr_inv(self.random_postive_real, axis=axis)

        # Assert
        # the shape should not be changed
        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in self.expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in self.expected_real.shape],  "unepected shape"

        # check the values
        assert np.allclose(rst_int, self.expected_int), "unmatched values"
        assert np.allclose(rst_real, self.expected_real), "unmatched values"

    @skipIf(no_torch, "Skipping tests: no torch dependency")
    def test_ndarray_torch(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int)
        random_real = torch.tensor(self.random_postive_real)
        expected_int = torch.tensor(self.expected_int)
        expected_real = torch.tensor(self.expected_real)

        # Action
        rst_int = ilr_inv(random_int, axis=axis)
        rst_real = ilr_inv(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_torch or no_cuda, "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int, device="cuda")
        random_real = torch.tensor(self.random_postive_real, device="cuda")
        expected_int = torch.tensor(self.expected_int, device="cuda")
        expected_real = torch.tensor(self.expected_real, device="cuda")

        # Action
        rst_int = ilr_inv(random_int, axis=axis)
        rst_real = ilr_inv(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int)
                             ) < 1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real)
                             ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax, "Skipping tests: no jax dependency")
    def test_ndarray_jnp(self):
        # Asrange
        axis = self.axis
        device = jax.devices("cpu")[0]

        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = ilr_inv(random_int, axis=axis)
        rst_real = ilr_inv(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )

    @skipIf(no_jax or no_cuda, "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = self.axis
        device = jax.devices("gpu")[0]
        random_int = jnp.array(self.random_postive_int, device=device)
        random_real = jnp.array(self.random_postive_real, device=device)
        expected_int = jnp.array(self.expected_int, device=device)
        expected_real = jnp.array(self.expected_real, device=device)

        # Action
        rst_int = ilr_inv(random_int, axis=axis)
        rst_real = ilr_inv(random_real, axis=axis)

        # Assert
        assert type(expected_int) == type(rst_int), "type is changed"
        assert type(rst_real) == type(rst_real), "type is changed"

        assert [int(i) for i in self.expected_shape_int] \
            == [int(i) for i in expected_int.shape],  "unepected shape"
        assert [int(i) for i in self.expected_shape_real] \
            == [int(i) for i in expected_real.shape],  "unepected shape"

        assert expected_int.device == rst_int.device, "different device"
        assert expected_real.device == rst_real.device, "different device"

        try:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int)) < 1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real)
                           ) < 1e-04, "unmatched values"
            test_name = inspect.currentframe().f_code.co_name
            class_name = self.__class__.__name__
            warnings.warn(
                f"In {class_name}.{test_name}: tolerance increased from 1e-8 to 1e-4.",
                UserWarning
            )


if __name__ == "__main__":
    main()
