# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Array API test infrastructure for unittest.TestCase.

This module provides tools to write a single test method that automatically
runs against multiple array API backends (NumPy, JAX, PyTorch, CuPy, etc.)
without duplicating test code.

Key components:

- ``@backends()`` decorator — runs a test method once per backend/device combo,
  using ``unittest.subTest()`` for per-backend reporting.
- ``ArrayAPITestMixin`` — mixin providing ``make_array()``, ``assert_close()``,
  and ``assert_type_preserved()`` helpers.

Environment variables:

- ``SKBIO_ARRAY_API_BACKEND`` — which backend(s) to test:
    - unset or empty: numpy only (fast default for development)
    - ``"jax"``, ``"torch"``, ``"cupy"``, etc.: that specific backend
    - ``"all"``: every installed backend
- ``SKBIO_DEVICE`` — device filter (``"cpu"``, ``"cuda"``, ``"gpu"``).
  Only relevant for backends with multi-device support.

Example usage::

    from unittest import TestCase
    from skbio.util._testing import ArrayAPITestMixin, backends

    class TestClosure(TestCase, ArrayAPITestMixin):

        @backends("numpy", "jax", "torch", "cupy")
        def test_basic(self, xp, device):
            data = self.make_array(xp, device, [[2.0, 2.0, 6.0]])
            result = closure(data)
            expected = self.make_array(xp, device, [[0.2, 0.2, 0.6]])
            self.assert_type_preserved(result, xp, device)
            self.assert_close(result, expected)

"""

import os
import functools
import unittest

import numpy as np
import numpy.testing as npt

from skbio.util._array import _to_numpy, _move_to_device, _get_backend_name


# -------------------------------------------------------------------------------------
# Environment configuration
# -------------------------------------------------------------------------------------


def _read_env():
    """Read backend/device environment variables.

    Returns
    -------
    tuple of (str, str or None)
        ``(backend, device)`` where backend is ``""`` if unset.

    """
    backend = os.environ.get("SKBIO_ARRAY_API_BACKEND", "").strip()
    device = os.environ.get("SKBIO_DEVICE", "").strip() or None
    return backend, device


# -------------------------------------------------------------------------------------
# Backend discovery
# -------------------------------------------------------------------------------------


def _get_available_backends():
    """Return dict of ``{name: (namespace, [devices])}`` for installed backends.

    Called once at module import time. Only backends that are actually installed
    will be included.

    """
    backends = {"numpy": (np, ["cpu"])}

    try:
        import jax
        import jax.numpy as jnp

        jax.config.update("jax_enable_x64", True)
        devices = ["cpu"]
        try:
            if jax.devices("gpu"):
                devices.append("gpu")
        except RuntimeError:
            pass
        backends["jax"] = (jnp, devices)
    except ImportError:
        pass

    try:
        import torch

        torch.set_default_dtype(torch.float64)
        devices = ["cpu"]
        if torch.cuda.is_available():
            devices.append("cuda")
        backends["torch"] = (torch, devices)
    except ImportError:
        pass

    try:
        import cupy

        backends["cupy"] = (cupy, ["cuda"])
    except ImportError:
        pass

    return backends


# Populated once at import time.
_BACKENDS = _get_available_backends()


def _should_run(backend_name, device):
    """Check whether a backend/device combo should execute.

    Decision is based on the ``SKBIO_ARRAY_API_BACKEND`` and ``SKBIO_DEVICE``
    environment variables.

    """
    env_backend, env_device = _read_env()

    # No env var set → numpy/cpu only (fast default).
    if not env_backend:
        return backend_name == "numpy" and device == "cpu"

    # "all" → run everything installed, optionally filtered by device.
    if env_backend == "all":
        if env_device and env_device != device:
            return False
        return True

    # Specific backend requested.
    if env_backend != backend_name:
        return False

    # Optionally filter by device.
    if env_device and env_device != device:
        return False

    return True


# -------------------------------------------------------------------------------------
# Backend-agnostic assertions
# -------------------------------------------------------------------------------------


def xp_assert_close(actual, desired, rtol=1e-7, atol=0):
    """Assert two arrays are element-wise equal within a tolerance.

    Converts both arrays to NumPy via :func:`~skbio.util._array._to_numpy`
    (handles GPU arrays correctly), then delegates to
    ``numpy.testing.assert_allclose``.

    Parameters
    ----------
    actual : array
        Array obtained.
    desired : array
        Array desired.
    rtol : float, optional
        Relative tolerance. Default is 1e-7.
    atol : float, optional
        Absolute tolerance. Default is 0.

    Raises
    ------
    AssertionError
        If arrays are not equal within the tolerance.

    """
    npt.assert_allclose(_to_numpy(actual), _to_numpy(desired), rtol=rtol, atol=atol)


def xp_assert_equal(actual, desired):
    """Assert two arrays are element-wise exactly equal.

    Converts both arrays to NumPy via :func:`~skbio.util._array._to_numpy`,
    then delegates to ``numpy.testing.assert_array_equal``.

    Parameters
    ----------
    actual : array
        Array obtained.
    desired : array
        Array desired.

    Raises
    ------
    AssertionError
        If arrays are not exactly equal.

    """
    npt.assert_array_equal(_to_numpy(actual), _to_numpy(desired))


# -------------------------------------------------------------------------------------
# @backends() decorator
# -------------------------------------------------------------------------------------


def backends(*backend_names, cpu_only=False):
    """Decorator: run a test method once per supported backend/device.

    Uses ``unittest.subTest()`` for per-backend reporting. Works with
    ``unittest.TestCase`` — no pytest migration needed.

    Parameters
    ----------
    *backend_names : str
        Which backends to run on (e.g. ``"numpy"``, ``"jax"``, ``"torch"``,
        ``"cupy"``). If empty, defaults to all four.
    cpu_only : bool, optional
        If True, skip GPU devices even if available.

    Raises
    ------
    RuntimeError
        If ``SKBIO_ARRAY_API_BACKEND`` names a specific backend that is not
        installed. This is intentionally a hard error so that CI jobs do not
        silently pass when a backend is missing.

    Notes
    -----
    The decorated method receives two extra arguments: ``xp`` (the array
    namespace) and ``device`` (a device string like ``"cpu"`` or ``"cuda"``).

    Example::

        class TestShannon(unittest.TestCase, ArrayAPITestMixin):

            @backends("numpy", "jax", "torch")
            def test_basic(self, xp, device):
                counts = self.make_array(xp, device, [1, 2, 3, 0, 5],
                                         dtype=xp.float64)
                result = shannon(counts)
                self.assert_type_preserved(result, xp, device)

    """
    if not backend_names:
        backend_names = ("numpy", "jax", "torch", "cupy")

    def decorator(test_func):
        @functools.wraps(test_func)
        def wrapper(self):
            env_backend, _ = _read_env()

            # Hard error: a specific backend was requested but isn't installed.
            if env_backend and env_backend != "all" and env_backend not in _BACKENDS:
                raise RuntimeError(
                    f"SKBIO_ARRAY_API_BACKEND={env_backend!r} but "
                    f"{env_backend!r} is not installed. "
                    f"Installed backends: {sorted(_BACKENDS.keys())}"
                )

            ran_any = False
            for name in backend_names:
                if name not in _BACKENDS:
                    continue
                xp, devices = _BACKENDS[name]
                for device in devices:
                    if cpu_only and device != "cpu":
                        continue
                    if not _should_run(name, device):
                        continue

                    ran_any = True
                    with self.subTest(backend=name, device=device):
                        test_func(self, xp, device)

            if not ran_any:
                raise unittest.SkipTest(f"No matching backends for {backend_names}")

        return wrapper

    return decorator


# -------------------------------------------------------------------------------------
# ArrayAPITestMixin
# -------------------------------------------------------------------------------------


class ArrayAPITestMixin:
    """Mixin providing array-API test helpers for ``unittest.TestCase``.

    Mix into your test class alongside ``TestCase``::

        class TestMyFunc(unittest.TestCase, ArrayAPITestMixin):
            ...

    Provides:

    - ``make_array(xp, device, data)`` — create an array on the right
      backend and device.
    - ``assert_close(actual, expected)`` — backend-agnostic numerical
      comparison (converts to NumPy internally).
    - ``assert_type_preserved(result, xp, device)`` — verify the result
      is from the correct backend and on the correct device.

    """

    def make_array(self, xp, device, data, dtype=None):
        """Create an array on the appropriate backend and device.

        Parameters
        ----------
        xp : module
            Array namespace (``numpy``, ``jax.numpy``, ``torch``, ``cupy``).
        device : str or None
            Device string (``"cpu"``, ``"cuda"``, ``"gpu"``).
        data : array_like
            Input data (list, NumPy array, scalar, etc.).
        dtype : dtype, optional
            Desired dtype. If ``None``, uses the backend's default.

        Returns
        -------
        array
            Array on the specified backend and device.

        """
        if dtype is not None:
            arr = xp.asarray(data, dtype=dtype)
        else:
            arr = xp.asarray(data)

        return _move_to_device(arr, xp, device)

    def assert_close(self, actual, expected, rtol=1e-7, atol=1e-7):
        """Assert two arrays are numerically close.

        Both arrays are converted to NumPy before comparison, so this works
        for any backend including GPU arrays.

        """
        npt.assert_allclose(
            _to_numpy(actual),
            _to_numpy(expected),
            rtol=rtol,
            atol=atol,
        )

    def assert_type_preserved(self, result, xp, device):
        """Assert result is from the correct backend and on the right device.

        Checks two things:

        1. The result array's module matches the expected namespace.
        2. If a non-CPU device was requested, the result is on that device.

        """
        result_name = _get_backend_name(
            __import__("array_api_compat").array_namespace(result)
        )
        expected_name = _get_backend_name(xp)
        self.assertEqual(
            result_name,
            expected_name,
            f"Backend not preserved: result is {result_name}, expected {expected_name}",
        )

        if device not in ("cpu", None) and hasattr(result, "device"):
            self.assertIn(
                device,
                str(result.device),
                f"Device not preserved: result on {result.device}, expected {device}",
            )
