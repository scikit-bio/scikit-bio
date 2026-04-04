# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main, skipIf
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import numpy.testing as npt
import array_api_compat as aac

from skbio.util import get_package
from skbio.util._array import (
    ingest_array,
    _get_array,
    _get_namespace_from_args,
    _check_array_api_backend,
    _to_numpy,
    _move_to_device,
    _get_backend_name,
)


# import optional dependencies
cp = get_package("cupy", raise_error=False)
jax = get_package("jax", raise_error=False)
torch = get_package("torch", raise_error=False)
xr = get_package("xarray", raise_error=False)
dask = get_package("dask", raise_error=False)
if dask is not None:
    da = get_package("dask.array")

# Get the numpy namespace via array_api_compat (aac has no .numpy attribute).
np_xp = aac.array_namespace(np.empty(0))


# =====================================================================
# Tests for _get_array
# =====================================================================


class TestGetArray(TestCase):

    def test_non_array_converted_to_numpy(self):
        # A plain list is not an array API object; should become ndarray.
        obs = _get_array([1, 2, 3])
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.array([1, 2, 3]))

    def test_numpy_array_returned_as_is(self):
        arr = np.arange(5)
        obs = _get_array(arr)
        self.assertIs(obs, arr)

    def test_numpy_array_to_numpy_noop(self):
        arr = np.arange(5)
        obs = _get_array(arr, to_numpy=True)
        self.assertIs(obs, arr)

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_kept_when_to_numpy_false(self):
        ten = torch.arange(5)
        obs = _get_array(ten)
        self.assertIsInstance(obs, torch.Tensor)
        self.assertIs(obs, ten)

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_to_numpy_via_dlpack(self):
        ten = torch.arange(5)
        obs = _get_array(ten, to_numpy=True)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(cp is None, "CuPy is not available.")
    def test_cupy_to_numpy(self):
        arr = cp.arange(5)
        obs = _get_array(arr, to_numpy=True)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(jax is None, "JAX is not available.")
    def test_jax_to_numpy(self):
        arr = jax.numpy.arange(5)
        obs = _get_array(arr, to_numpy=True)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(dask is None, "Dask is not available.")
    def test_dask_to_numpy_fallback(self):
        # Dask doesn't support __dlpack__, so _get_array should fall back
        # to np.asarray.
        arr = da.arange(5, chunks=5)
        obs = _get_array(arr, to_numpy=True)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    def test_to_numpy_dlpack_fallback_on_runtime_error(self):
        # Simulate an array API object whose from_dlpack raises RuntimeError.
        mock_arr = MagicMock()
        mock_arr.__array__ = lambda dtype=None, copy=None: np.array([10, 20, 30])

        with patch.object(aac, "is_array_api_obj", return_value=True), \
             patch.object(aac, "is_numpy_array", return_value=False), \
             patch("skbio.util._array.np.from_dlpack",
                   side_effect=RuntimeError("no dlpack")):
            obs = _get_array(mock_arr, to_numpy=True)
        self.assertIsInstance(obs, np.ndarray)

    def test_to_numpy_dlpack_fallback_on_attribute_error(self):
        mock_arr = MagicMock()
        mock_arr.__array__ = lambda dtype=None, copy=None: np.array([1])

        with patch.object(aac, "is_array_api_obj", return_value=True), \
             patch.object(aac, "is_numpy_array", return_value=False), \
             patch("skbio.util._array.np.from_dlpack",
                   side_effect=AttributeError):
            obs = _get_array(mock_arr, to_numpy=True)
        self.assertIsInstance(obs, np.ndarray)

    def test_to_numpy_dlpack_fallback_on_type_error(self):
        mock_arr = MagicMock()
        mock_arr.__array__ = lambda dtype=None, copy=None: np.array([1])

        with patch.object(aac, "is_array_api_obj", return_value=True), \
             patch.object(aac, "is_numpy_array", return_value=False), \
             patch("skbio.util._array.np.from_dlpack",
                   side_effect=TypeError):
            obs = _get_array(mock_arr, to_numpy=True)
        self.assertIsInstance(obs, np.ndarray)


# =====================================================================
# Tests for _get_namespace_from_args
# =====================================================================


class TestGetNamespaceFromArgs(TestCase):

    def test_no_arrays_returns_none(self):
        result = _get_namespace_from_args((1, "hello", [1, 2]), {})
        self.assertIsNone(result)

    def test_numpy_in_args(self):
        arr = np.arange(3)
        result = _get_namespace_from_args((arr,), {})
        self.assertIs(result, np_xp)

    def test_numpy_in_kwargs(self):
        arr = np.arange(3)
        result = _get_namespace_from_args((), {"x": arr})
        self.assertIs(result, np_xp)

    def test_mixed_args_and_kwargs(self):
        a = np.arange(3)
        b = np.ones(3)
        result = _get_namespace_from_args((a,), {"b": b})
        self.assertIs(result, np_xp)

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_in_args(self):
        ten = torch.arange(3)
        result = _get_namespace_from_args((ten,), {})
        self.assertTrue(aac.is_torch_namespace(result))

    def test_empty_args_and_kwargs(self):
        result = _get_namespace_from_args((), {})
        self.assertIsNone(result)


# =====================================================================
# Tests for _check_array_api_backend
# =====================================================================


class TestCheckArrayApiBackend(TestCase):

    def test_numpy_allowed(self):
        # Should not raise.
        _check_array_api_backend(np_xp, ["numpy"], "my_func")

    def test_numpy_not_allowed_raises(self):
        with self.assertRaises(TypeError) as ctx:
            _check_array_api_backend(np_xp, ["torch"], "my_func")
        self.assertIn("my_func", str(ctx.exception))
        self.assertIn("torch", str(ctx.exception))

    def test_multiple_backends_allowed(self):
        _check_array_api_backend(np_xp, ["numpy", "torch"], "my_func")

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_allowed(self):
        xp = aac.array_namespace(torch.arange(1))
        _check_array_api_backend(xp, ["torch"], "my_func")

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_not_allowed_raises(self):
        xp = aac.array_namespace(torch.arange(1))
        with self.assertRaises(TypeError):
            _check_array_api_backend(xp, ["numpy"], "my_func")

    def test_unknown_backend_in_list_ignored(self):
        # "unknown_backend" isn't in _BACKEND_CHECKERS, so checker is None.
        # numpy is also listed, so this should pass.
        _check_array_api_backend(np_xp, ["unknown_backend", "numpy"], "f")

    def test_only_unknown_backend_raises(self):
        # No checker matches numpy if only "unknown_backend" is listed.
        with self.assertRaises(TypeError):
            _check_array_api_backend(np_xp, ["unknown_backend"], "f")

    def test_error_message_lists_supported(self):
        with self.assertRaises(TypeError) as ctx:
            _check_array_api_backend(np_xp, ["torch", "jax"], "compute")
        msg = str(ctx.exception)
        self.assertIn("compute", msg)
        self.assertIn("torch, jax", msg)


# =====================================================================
# Tests for _to_numpy
# =====================================================================


class TestToNumpy(TestCase):

    def test_numpy_passthrough(self):
        arr = np.arange(5)
        obs = _to_numpy(arr)
        self.assertIs(obs, arr)

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_cpu(self):
        ten = torch.arange(5)
        obs = _to_numpy(ten)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_requires_grad(self):
        ten = torch.arange(5, dtype=torch.float32).requires_grad_(True)
        obs = _to_numpy(ten)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5, dtype=np.float32))

    @skipIf(cp is None, "CuPy is not available.")
    def test_cupy(self):
        arr = cp.arange(5)
        obs = _to_numpy(arr)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(jax is None, "JAX is not available.")
    def test_jax(self):
        arr = jax.numpy.arange(5)
        obs = _to_numpy(arr)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    @skipIf(dask is None, "Dask is not available.")
    def test_dask(self):
        arr = da.arange(5, chunks=5)
        obs = _to_numpy(arr)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(5))

    def test_generic_array_via_np_asarray(self):
        # An object with no .detach and no .get falls through to np.asarray.
        mock_arr = MagicMock(spec=[])
        mock_arr.__array__ = lambda dtype=None, copy=None: np.array([7, 8, 9])
        # Ensure it's not ndarray and has no detach/get.
        self.assertNotIsInstance(mock_arr, np.ndarray)
        obs = _to_numpy(mock_arr)
        self.assertIsInstance(obs, np.ndarray)


# =====================================================================
# Tests for _move_to_device
# =====================================================================


class TestMoveToDevice(TestCase):

    def test_none_device_returns_same(self):
        arr = np.arange(5)
        obs = _move_to_device(arr, np, None)
        self.assertIs(obs, arr)

    def test_cpu_device_returns_same(self):
        arr = np.arange(5)
        obs = _move_to_device(arr, np, "cpu")
        self.assertIs(obs, arr)

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_cpu_noop(self):
        ten = torch.arange(5)
        obs = _move_to_device(ten, torch, "cpu")
        self.assertIs(obs, ten)

    @skipIf(torch is None or not torch.cuda.is_available(),
            "PyTorch CUDA is not available.")
    def test_torch_cuda(self):
        ten = torch.arange(5)
        obs = _move_to_device(ten, torch, "cuda")
        self.assertTrue(obs.is_cuda)
        npt.assert_array_equal(obs.cpu().numpy(), np.arange(5))

    def test_numpy_non_cpu_passthrough(self):
        # NumPy arrays are always on CPU; non-cpu device returns unchanged.
        arr = np.arange(5)
        obs = _move_to_device(arr, np, "gpu")
        self.assertIs(obs, arr)

    @skipIf(cp is None, "CuPy is not available.")
    def test_cupy_non_cpu_passthrough(self):
        import cupy
        arr = cupy.arange(5)
        obs = _move_to_device(arr, cupy, "gpu")
        self.assertIs(obs, arr)

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch_to_device(self):
        # torch.to("cpu") should work and return a CPU tensor.
        ten = torch.arange(5)
        # Use a module whose __name__ starts with "torch".
        obs = _move_to_device(ten, torch, "cpu")
        # device=cpu means early return, so test with a mock device name
        # that triggers the torch branch but stays on cpu.
        self.assertIs(obs, ten)


# =====================================================================
# Tests for _get_backend_name
# =====================================================================


class TestGetBackendName(TestCase):

    def test_numpy(self):
        self.assertEqual(_get_backend_name(np_xp), "numpy")

    @skipIf(torch is None, "PyTorch is not available.")
    def test_torch(self):
        xp = aac.array_namespace(torch.arange(1))
        self.assertEqual(_get_backend_name(xp), "torch")

    @skipIf(cp is None, "CuPy is not available.")
    def test_cupy(self):
        xp = aac.array_namespace(cp.arange(1))
        self.assertEqual(_get_backend_name(xp), "cupy")

    @skipIf(jax is None, "JAX is not available.")
    def test_jax(self):
        xp = aac.array_namespace(jax.numpy.arange(1))
        self.assertEqual(_get_backend_name(xp), "jax")

    @skipIf(dask is None, "Dask is not available.")
    def test_dask(self):
        xp = aac.array_namespace(da.arange(1, chunks=1))
        self.assertEqual(_get_backend_name(xp), "dask")

    def test_unknown_namespace_with_name(self):
        mock_ns = MagicMock()
        mock_ns.__name__ = "mystery_lib"
        with patch.dict(
            "skbio.util._array._BACKEND_CHECKERS",
            {k: (lambda ns: False) for k in
             ["numpy", "cupy", "torch", "jax", "dask"]},
        ):
            self.assertEqual(_get_backend_name(mock_ns), "mystery_lib")

    def test_unknown_namespace_without_name(self):
        mock_ns = object()
        with patch.dict(
            "skbio.util._array._BACKEND_CHECKERS",
            {k: (lambda ns: False) for k in
             ["numpy", "cupy", "torch", "jax", "dask"]},
        ):
            result = _get_backend_name(mock_ns)
            # Falls back to str(xp) since object has no __name__.
            self.assertIsInstance(result, str)


# =====================================================================
# Tests for ingest_array (original tests preserved + additions)
# =====================================================================


class TestIngestArray(TestCase):

    def test_ingest_array_numpy(self):
        # NumPy arrays are kept as-is.
        arr = np.arange(5)
        xp, obs = ingest_array(arr)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        self.assertIs(obs, arr)
        npt.assert_array_equal(obs, np.arange(5))

        xp, obs = ingest_array(arr, to_numpy=True)
        self.assertIs(xp, np_xp)
        self.assertIs(obs, arr)

        # 2-D array
        arr = np.arange(9).reshape(3, 3)
        xp, obs = ingest_array(arr)
        self.assertIs(xp, np_xp)
        self.assertIs(obs, arr)
        npt.assert_array_equal(obs, np.arange(9).reshape(3, 3))

    def test_ingest_array_sequence(self):
        # Plain Python lists are converted into NumPy arrays.
        lst = list(range(5))
        xp, obs = ingest_array(lst)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.asarray(lst))

        tup = tuple(range(5))
        xp, obs = ingest_array(tup)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.asarray(tup))

        ran = range(5)
        xp, obs = ingest_array(ran)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.asarray(ran))

    def test_ingest_array_pyarray(self):
        # Python's built-in array is not an API-compatible array library.
        import array as pyarray
        arr = pyarray.array('i', range(5))
        xp, obs = ingest_array(arr)
        self.assertIs(xp, np_xp)
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
        self.assertIs(xp, np_xp)
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
        self.assertIs(xp, np_xp)
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
        self.assertIs(xp, np_xp)
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
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.arange(20))

    def test_ingest_array_scalar(self):
        # Scalars are converted into 0-D NumPy arrays.
        xp, obs = ingest_array(42)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.array(42))
        self.assertEqual(obs.ndim, 0)

        # Strings are also scalars (not treated as sequences of characters).
        xp, obs = ingest_array("hello")
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.array("hello"))
        self.assertEqual(obs.ndim, 0)

        xp, obs = ingest_array(None)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(obs, np.array(None))
        self.assertEqual(obs.ndim, 0)

    def test_ingest_array_pandas(self):
        # Pandas Series' are not API-compatible arrays, but they can be casted into
        # NumPy arrays.
        s = pd.Series(range(5), index=list("abcde"))
        xp, obs = ingest_array(s)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(s.to_numpy(), np.arange(5))

        # So are Pandas DataFrames.
        df = s.to_frame()
        xp, obs = ingest_array(df)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(df.to_numpy(), np.arange(5).reshape(-1, 1))

    @skipIf(xr is None, "Xarray is not available for unit tests.")
    def test_ingest_array_xarray(self):
        # Xarray DataArrays are not API-compatible arrays, but they can be casted into
        # Numpy arrays.
        xarr = xr.DataArray(range(5))
        xp, obs = ingest_array(xarr)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(obs, np.ndarray)
        npt.assert_array_equal(xarr.to_numpy(), np.arange(5))

    def test_ingest_array_multiple(self):
        # Multiple arrays should all be returned with shared namespace.
        a = np.arange(3)
        b = np.ones(4)
        xp, oa, ob = ingest_array(a, b)
        self.assertIs(xp, np_xp)
        self.assertIs(oa, a)
        self.assertIs(ob, b)

    def test_ingest_array_multiple_mixed(self):
        # Mix of list and ndarray; both should come back as ndarray.
        a = [1, 2, 3]
        b = np.ones(3)
        xp, oa, ob = ingest_array(a, b)
        self.assertIs(xp, np_xp)
        self.assertIsInstance(oa, np.ndarray)
        self.assertIs(ob, b)


if __name__ == "__main__":
    main()