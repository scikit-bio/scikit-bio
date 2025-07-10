# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase,  main, skipIf
from copy import deepcopy
# switched to skbio warning
import warnings

import numpy as np
import torch
import numpy.testing as npt
from numpy.exceptions import AxisError
from numpy.random import normal, rand, randint
import pandas as pd
import pandas.testing as pdt
from scipy.sparse import coo_matrix
from scipy.stats import f_oneway, ConstantInputWarning

from skbio import TreeNode
from skbio.util import assert_data_frame_almost_equal
from skbio.stats.distance import DistanceMatrixError
from skbio.stats.composition import (
    _check_composition, _check_orthogonality,
    closure, multi_replace, perturb, perturb_inv, power, inner, clr, clr_inv, ilr,
    ilr_inv, alr, alr_inv, sbp_basis, _gram_schmidt_basis, centralize, _check_sig_test,
    _check_p_adjust, ancom, vlr, pairwise_vlr, tree_basis, dirmult_ttest, dirmult_lme)

try:
    import jax
    import jax.numpy as jnp
    no_jax = False
except ImportError:
    no_jax = True
    
try:
    import torch
    no_torch = False
except ImportError:
    no_torch = True

import subprocess
def no_cuda_available():
    try:
        result = subprocess.run(
            ["nvidia-smi"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return not result.returncode == 0
    except FileNotFoundError:
        return True

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
        self.random_postive_int = randint(1,100, (4,6,3))
        self.random_postive_real = rand(4,6,3)+1e-04
    
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
        
        assert torch.sum(torch.abs(expected_int - rst_int))<1e-08, "unmatched values"
        assert torch.sum(torch.abs(expected_real - rst_real))<1e-08, "unmatched values"
    
    @skipIf(no_torch or no_cuda_available(), "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = 1
        random_int = torch.tensor(self.random_postive_int, device='cuda')
        random_real = torch.tensor(self.random_postive_real, device='cuda')
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
        
        assert torch.sum(torch.abs(expected_int - rst_int))<1e-08, "unmatched values"
        assert torch.sum(torch.abs(expected_real - rst_real))<1e-08, "unmatched values"

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
        
        assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-08, "unmatched values"
        assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-08, "unmatched values"
        
    @skipIf(no_jax or no_cuda_available(), "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = 1
        device = jax.devices('gpu')[0]
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
        
        assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-08, "unmatched values"
        assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-08, "unmatched values"

        
class Tests_CLR(TestCase):
    def setUp(self):
        self.random_postive_int = randint(1,100, (4,6,3))
        self.random_postive_real = rand(4,6,3)+1e-04
    
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
            assert torch.sum(torch.abs(expected_int - rst_int))<1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real))<1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int))<1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real))<1e-04, "unmatched values"
            warnings.warn("tolerance is increased from 1e-8 to 1e-4.")
    
    @skipIf(no_torch or no_cuda_available(), "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = 1
        random_int = torch.tensor(self.random_postive_int, device='cuda')
        random_real = torch.tensor(self.random_postive_real, device='cuda')
        expected_int = (lmat := np.log(self.random_postive_int)) \
                                    - np.mean(lmat, axis=axis, keepdims=True)
        expected_int = torch.tensor(expected_int, device='cuda')
        expected_real = (lmat := np.log(self.random_postive_real)) \
                                    - np.mean(lmat, axis=axis, keepdims=True)
        expected_real = torch.tensor(expected_real, device='cuda')
        
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
            assert torch.sum(torch.abs(expected_int - rst_int))<1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real))<1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int))<1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real))<1e-04, "unmatched values"
            warnings.warn("tolerance is increased from 1e-8 to 1e-4.")

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
            assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-04, "unmatched values"
            warnings.warn("tolerance is increased from 1e-8 to 1e-4.")
        
        
    @skipIf(no_jax or no_cuda_available(), "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = 1
        device = jax.devices('gpu')[0]
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
            assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-04, "unmatched values"
            warnings.warn("tolerance is increased from 1e-8 to 1e-4.")


class Tests_CLR_INV(TestCase):
    def setUp(self):
        axis = 1
        self.axis = axis
        self.random_postive_int = randint(1,100, (4,6,3))
        self.random_postive_real = rand(4,6,3)+1e-04
        
        tem_int = np.exp(self.random_postive_int- \
                            np.max(self.random_postive_int, axis=axis, keepdims=True))
        self.expected_int = tem_int/np.sum(tem_int, axis=axis, keepdims=True)
        
        tem_real = np.exp(self.random_postive_real- \
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
        np.allclose(rst_int, self.expected_int)
        np.allclose(rst_real, self.expected_real)
    
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
            assert torch.sum(torch.abs(expected_int - rst_int))<1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real))<1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int))<1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real))<1e-04, "unmatched values"
            warnings.warn("tolerance is increased from 1e-8 to 1e-4.")
        
    
    @skipIf(no_torch or no_cuda_available(), "Skipping tests: no torch dependency or no cuda")
    def test_ndarray_torch_cuda(self):
        # Asrange
        axis = self.axis
        random_int = torch.tensor(self.random_postive_int, device='cuda')
        random_real = torch.tensor(self.random_postive_real, device='cuda')
        expected_int = torch.tensor(self.expected_int, device='cuda')
        expected_real = torch.tensor(self.expected_real, device='cuda')
        
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
            assert torch.sum(torch.abs(expected_int - rst_int))<1e-08, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real))<1e-08, "unmatched values"
        except:
            assert torch.sum(torch.abs(expected_int - rst_int))<1e-04, "unmatched values"
            assert torch.sum(torch.abs(expected_real - rst_real))<1e-04, "unmatched values"
            warnings.warn("tolerance is increased from 1e-8 to 1e-4.")

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
            assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-04, "unmatched values"
            warnings.warn("tolerance is increased from 1e-8 to 1e-4.")
        

    @skipIf(no_jax or no_cuda_available(), "Skipping tests: no jax dependency or no cuda")
    def test_ndarray_jnp_gpu(self):
        # Asrange
        axis = self.axis
        device = jax.devices('gpu')[0]
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
            assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-08, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-08, "unmatched values"
        except:
            assert jnp.sum(jnp.abs(expected_int - rst_int))<1e-04, "unmatched values"
            assert jnp.sum(jnp.abs(expected_real - rst_real))<1e-04, "unmatched values"
            warnings.warn("tolerance is increased from 1e-8 to 1e-4.")


if __name__ == '__main__':
    main()