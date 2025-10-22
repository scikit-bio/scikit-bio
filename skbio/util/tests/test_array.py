# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main, skipIf

import numpy as np
import pandas as pd
import numpy.testing as npt
import array_api_compat as aac

from skbio.util import get_package
from skbio.util._array import ingest_array


# import optional dependencies
cp = get_package("cupy", raise_error=False)
jax = get_package("jax", raise_error=False)
torch = get_package("torch", raise_error=False)
xr = get_package("xarray", raise_error=False)
dask = get_package("dask", raise_error=False)
if dask is not None:
    da = get_package("dask.array")


class TestIngestArray(TestCase):

    def test_ingest_array_numpy(self):
        # NumPy arrays are kept as-is.
        arr = np.arange(5)
        xp, obs = ingest_array(arr)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        self.assertIs(obs, arr)
        npt.assert_array_equal(obs, np.arange(5))

        xp, obs = ingest_array(arr, to_numpy=True)
        self.assertIs(xp, aac.numpy)
        self.assertIs(obs, arr)

        # 2-D array
        arr = np.arange(9).reshape(3, 3)
        xp, obs = ingest_array(arr)
        self.assertIs(xp, aac.numpy)
        self.assertIs(obs, arr)
        npt.assert_array_equal(obs, np.arange(9).reshape(3, 3))

    def test_ingest_array_sequence(self):
        # Plain Python lists are converted into NumPy arrays.
        lst = list(range(5))
        xp, obs = ingest_array(lst)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.asarray(lst))

        tup = tuple(range(5))
        xp, obs = ingest_array(tup)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.asarray(tup))

        ran = range(5)
        xp, obs = ingest_array(ran)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.asarray(ran))

    def test_ingest_array_pyarray(self):
        # Python's built-in array is not an API-compatible array library.
        import array as pyarray
        arr = pyarray.array('i', range(5))
        xp, obs = ingest_array(arr)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.asarray(arr))

    @skipIf(cp is None, "CuPy is not available for unit tests.")
    def test_ingest_array_cupy(self):
        # CuPy arrays are kept as-is in default mode.
        arr = cp.arange(5)
        xp, obs = ingest_array(arr)
        self.assertTrue(aac.is_cupy_namespace(xp))
        self.assertIsInstance(obs, cp.ndarray)
        self.assertIs(obs, arr)
        npt.assert_array_equal(cp.asnumpy(obs), np.arange(5))

        # Convert CuPy arrays into NumPy arrays.
        xp, obs = ingest_array(arr, to_numpy=True)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(torch is None, "PyTorch is not available for unit tests.")
    def test_ingest_array_torch(self):
        # PyTorch tensors are kept as-is in default mode.
        ten = torch.arange(5)
        xp, obs = ingest_array(ten)
        self.assertTrue(aac.is_torch_namespace(xp))
        self.assertIsInstance(obs, torch.Tensor)
        self.assertIs(obs, ten)
        npt.assert_array_equal(obs.numpy(), np.arange(5))

        # Convert PyTorch tensors into NumPy arrays.
        xp, obs = ingest_array(ten, to_numpy=True)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(jax is None, "JAX is not available for unit tests.")
    def test_ingest_array_jax(self):
        # JAX arrays are kept as-is in default mode.
        arr = jax.numpy.arange(5)
        xp, obs = ingest_array(arr)
        self.assertTrue(aac.is_jax_namespace(xp))
        self.assertIsInstance(obs, jax.Array)
        self.assertIs(obs, arr)
        npt.assert_array_equal(np.asarray(obs), np.arange(5))

        # Convert JAX arrays into NumPy arrays.
        xp, obs = ingest_array(arr, to_numpy=True)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(dask is None, "Dask is not available for unit tests.")
    def test_ingest_array_dask(self):
        # Dask arrays are kept as-is in default mode.
        arr = da.arange(20, chunks=4)
        xp, obs = ingest_array(arr)
        self.assertTrue(aac.is_dask_namespace(xp))
        self.assertIsInstance(obs, da.Array)
        self.assertIs(obs, arr)
        npt.assert_array_equal(obs.compute(), np.arange(20))

        # Convert Dask arrays into NumPy arrays.
        xp, obs = ingest_array(arr, to_numpy=True)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(20))

    def test_ingest_array_scalar(self):
        # Scalars are converted into 0-D NumPy arrays.
        xp, obs = ingest_array(42)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.array(42))
        self.assertEqual(obs.ndim, 0)

        # Strings are also scalars (not treated as sequences of characters).
        xp, obs = ingest_array("hello")
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.array("hello"))
        self.assertEqual(obs.ndim, 0)

        xp, obs = ingest_array(None)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.array(None))
        self.assertEqual(obs.ndim, 0)

    def test_ingest_array_pandas(self):
        # Pandas Series' are not API-compatible arrays, but they can be casted into
        # NumPy arrays.
        s = pd.Series(range(5), index=list("abcde"))
        xp, obs = ingest_array(s)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(s.to_numpy(), np.arange(5))

        # So are Pandas DataFrames.
        df = s.to_frame()
        xp, obs = ingest_array(df)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(df.to_numpy(), np.arange(5).reshape(-1, 1))

    @skipIf(xr is None, "Xarray is not available for unit tests.")
    def test_ingest_array_xarray(self):
        # Xarray DataArrays are not API-compatible arrays, but they can be casted into
        # Numpy arrays.
        xarr = xr.DataArray(range(5))
        xp, obs = ingest_array(xarr)
        self.assertIs(xp, aac.numpy)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(xarr.to_numpy(), np.arange(5))


if __name__ == "__main__":
    main()
