import unittest
from unittest import TestCase, skipIf

import numpy as np
import numpy.testing as npt

from skbio.util import get_package
from skbio.util._gpu import cuda_avail

from skbio.diversity.alpha import shannon, simpson, dominance, simpson_e

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
    torch.set_default_dtype(torch.float64)

no_cuda = not cuda_avail()


class Tests_Alpha_ArrayAPI(TestCase):
    def setUp(self):
        np.random.seed(0)
        self.counts = np.random.randint(0, 50, (6, 8))
        self.counts[0, [1, 5]] = 0
        self.counts[2, [0, 7]] = 0
        self.counts[5, :] = 0

        self.exp_shannon = np.array([shannon(r) for r in self.counts])
        self.exp_simpson = np.array([simpson(r) for r in self.counts])
        self.exp_dominance = np.array([dominance(r) for r in self.counts])
        self.exp_simpson_e = np.array([simpson_e(r) for r in self.counts])

    def test_numpy_2d(self):
        npt.assert_allclose(shannon(self.counts, axis=1), self.exp_shannon)
        npt.assert_allclose(simpson(self.counts, axis=1), self.exp_simpson)
        npt.assert_allclose(dominance(self.counts, axis=1), self.exp_dominance)
        npt.assert_allclose(simpson_e(self.counts, axis=1), self.exp_simpson_e)

    @skipIf(no_torch, "no torch")
    def test_torch(self):
        mat = torch.tensor(self.counts)
        res = shannon(mat, axis=1)
        exp = torch.tensor(self.exp_shannon)
        assert type(res) == type(exp)
        assert res.device == exp.device
        torch.testing.assert_close(res, exp, equal_nan=True)

    @skipIf(no_torch or no_cuda, "no torch or cuda")
    def test_torch_cuda(self):
        mat = torch.tensor(self.counts, device="cuda")
        res = shannon(mat, axis=1)
        exp = torch.tensor(self.exp_shannon, device="cuda")
        assert type(res) == type(exp)
        assert res.device == exp.device
        torch.testing.assert_close(res, exp, equal_nan=True)

    @skipIf(no_jax, "no jax")
    def test_jnp(self):
        device = jax.devices("cpu")[0]
        mat = jnp.array(self.counts, device=device)
        res = shannon(mat, axis=1)
        exp = jnp.array(self.exp_shannon, device=device)
        assert type(res) == type(exp)
        assert res.device == exp.device
        npt.assert_allclose(np.array(res), np.array(exp))

    @skipIf(no_jax or no_cuda, "no jax or cuda")
    def test_jnp_gpu(self):
        device = jax.devices("gpu")[0]
        mat = jnp.array(self.counts, device=device)
        res = shannon(mat, axis=1)
        exp = jnp.array(self.exp_shannon, device=device)
        assert type(res) == type(exp)
        assert res.device == exp.device
        npt.assert_allclose(np.array(res), np.array(exp))


if __name__ == "__main__":
    unittest.main()
