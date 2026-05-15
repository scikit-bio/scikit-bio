Array API and GPU support
=========================

Scikit-bio is gradually adopting the `Python array API standard
<https://data-apis.org/array-api/latest/>`_ to allow functions to work
across multiple array libraries (NumPy, PyTorch, JAX, CuPy, Dask) and run on
GPUs where supported. This page covers what contributors need to know
to write, test, and run array-API-compatible code.

Overview
--------

The array API standard defines a common interface that several array
libraries implement. Functions written against this standard accept
arrays from any compliant library and dispatch to that library's
implementation — meaning the same function can run on CPU via NumPy,
on GPU via PyTorch or CuPy, or on TPU via JAX.

Currently supported backends:

- **NumPy** (CPU only) — the default, always available
- **Torch** — CPU or CUDA
- **Jax** — CPU or GPU
- **CuPy** — CUDA only

Only a small fraction of scikit-bio currently has array-API support.
Most functions still operate on NumPy arrays exclusively. This is being
expanded over time.


Writing array API compatible code
-----------------------------------------

Use ``ingest_array`` at function entry points to obtain the array
namespace (``xp``) and a converted array:

.. code-block:: python

   from skbio.util._array import ingest_array

   def my_function(arr):
       xp, arr = ingest_array(arr)
       result = xp.sum(arr, axis=0)
       return result

The returned ``xp`` is the namespace corresponding to the input array's
backend. Use ``xp.func()`` for all array operations rather than
``np.func()`` — this is what makes the function backend-agnostic.

.. code-block:: python

   def my_function(arr):
       xp, arr = ingest_array(arr)
       return xp.sum(arr)
       ...


**Things to avoid:**

- Hardcoded ``np.func()`` calls (use ``xp.func()`` instead).
- Backend-specific methods like ``torch.numpy()`` or ``cupy.get()`` (both return a NumPy array).
- In-place modification when supporting JAX arrays. They are immutable.
- Assuming ``float64`` is the default dtype (PyTorch defaults to
  ``float32``).


Writing tests for array API code
--------------------------------

Tests that exercise array-API code paths use the
``ArrayAPITestMixin`` and the ``@array_backends`` decorator. Write the
test once; the decorator runs it across every supported backend:

.. code-block:: python

   from unittest import TestCase
   from skbio.util._testing import ArrayAPITestMixin, array_backends

   class TestMyFunction(TestCase, ArrayAPITestMixin):

       @array_backends("numpy", "torch", "jax", "cupy")
       def test_my_function(self, xp, device):
           data = [[1.0, 2.0], [3.0, 4.0]]
           arr = self.make_array(xp, device, data)
           result = my_function(arr)
           expected = self.make_array(xp, device, [[4.0, 6.0]])
           self.assert_close(result, expected)
           self.assert_type_preserved(result, xp, device)

The decorator injects ``xp`` (the array namespace) and ``device`` (the
target device) into the test method. ``ArrayAPITestMixin`` provides:

- ``make_array(xp, device, data)`` — create a test array on the right
  backend and device
- ``assert_close(actual, expected)`` — numerical comparison that works
  across backends
- ``assert_type_preserved(result, xp, device)`` — verify the result
  came back on the expected backend and device

**When to use ``@array_backends``:** decorate tests that exercise the
array computation itself. Tests for input validation, error handling,
or pure-Python logic don't need it — those should remain standard
``TestCase`` methods.


Running tests locally
---------------------

Two environment variables control which backend and device tests use:

- ``SKBIO_ARRAY_BACKEND``: ``numpy`` (default), ``torch``, ``jax``,
  ``cupy``, or ``all``
- ``SKBIO_DEVICE``: ``cpu`` (default), ``cuda``, or ``gpu``

Common invocations:

.. code-block:: bash

   # Default: numpy on CPU
   pytest skbio/stats/composition/tests/test_base.py

   # Single backend on CPU (no GPU needed)
   SKBIO_ARRAY_BACKEND=torch pytest skbio/stats/composition/tests/test_base.py

   # All available backends on CPU and GPU (if available)
   SKBIO_ARRAY_BACKEND=all pytest skbio/stats/composition/tests/test_base.py

   # Run only array-API tests across the package
   pytest -m array_api skbio/

   # GPU run (requires CUDA-capable hardware)
   SKBIO_ARRAY_BACKEND=torch SKBIO_DEVICE=cuda pytest -m array_api skbio/

Tests for unavailable backends skip cleanly — you don't need every
backend installed. CuPy is GPU-only and will skip on machines without
CUDA.

The ``-m array_api`` marker filters to tests decorated with
``@array_backends``. Without it, pytest runs the full suite.


Continuous integration
----------------------

GPU testing happens in the ``Array API Compatibility`` workflow
(``.github/workflows/array-api.yml``), which runs on a T4 GPU runner
against the ``torch``, ``jax``, and ``cupy`` backends with CUDA.

**Triggers:**

- **Weekly cron**: every Monday at 06:00 UTC, catching regressions
  from upstream changes in array libraries
- **Manual dispatch**: via the Actions tab on GitHub
- **`gpu-ci` label on a PR**: applies the label to a PR to trigger
  the workflow; pushes to the labeled PR will re-trigger it

Reviewers should apply the ``gpu-ci`` label to any PR that modifies
array-API code paths. Most PRs do not need GPU validation and the
workflow has a non-trivial cost, so it is not run on every PR
automatically.

**What the workflow verifies:**

- Each backend installs successfully with CUDA support
- The GPU is visible to each backend at import time
- Decorated tests pass with arrays placed on CUDA, and results remain
  on CUDA after computation

**What it does not verify:**

- Code paths not covered by ``@array_backends``-decorated tests
- Subtle silent CPU fallback within individual operations
- Correctness on GPU architectures other than T4


Backend-specific gotchas
------------------------

A few cross-backend differences come up often enough to call out:

**PyTorch**

- Defaults to ``float32``, not ``float64``. Specify dtypes explicitly
  if precision matters.
- Tensors with ``requires_grad=True`` cannot be converted directly to
  NumPy; ``_to_numpy`` handles this by calling ``.detach()``.

**JAX**

- Arrays are immutable. In-place updates (``arr[i] = x``) raise an
  error; use ``arr.at[i].set(x)`` instead.
- Silently falls back to CPU if no GPU is detected at import time.
  Always check ``jax.devices("gpu")`` before assuming GPU execution.
- Uses ``"gpu"`` as its device string where Torch and CuPy use
  ``"cuda"``. These are treated as equivalent by
  ``_device_specs_match``.

**CuPy**

- GPU-only — there is no CPU mode. ``cupy.ndarray`` is always on a
  CUDA device.
- Requires a CUDA driver matching the installed CuPy build (e.g.
  ``cupy-cuda12x`` requires a CUDA 12-capable driver).

**General**

- Integer division and dtype promotion rules differ subtly across
  backends. When in doubt, cast explicitly.
- Some array API standard functions are unimplemented in certain
  backends. Check the library's documentation if you hit
  ``AttributeError`` on ``xp``.
