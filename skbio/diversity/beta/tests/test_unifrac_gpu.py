# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import importlib.util
import sys
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

from skbio.diversity.beta._unifrac_gpu import (
    _KERNEL_CACHE,
    _TILE_CANDIDATES,
    _TILE_CONFIG_CACHE,
    _build_pair_index,
    _get_tile_config,
    _make_unifrac_kernel,
    _probe_tile_config,
    _probe_tile_config_candidate,
    _detect_gpu_backend,
    _generalized_unifrac_gpu,
    _get_cuda_module,
    _unweighted_unifrac_gpu,
    _weighted_unifrac_gpu,
)
from skbio.diversity.beta.tests._fixtures import QiimeTinyTestMixin


def _cleanup_kernel_cache(fake_cuda):
    """Drop every `_KERNEL_CACHE` entry a fake cuda module may have created.

    Entries are keyed by ``(id(module), tile, node_chunk, variance_adjust)``,
    and a fake module's id can be reused by a later object once it is
    garbage collected, so its entries must not outlive the test that made
    them.
    """
    for tile, node_chunk in _TILE_CANDIDATES:
        for variance_adjust in (True, False):
            _KERNEL_CACHE.pop(
                (id(fake_cuda), tile, node_chunk, variance_adjust), None
            )


# Measured max abs deviation, GPU vs CPU-numba, across all 12
# method/normalized/variance_adjust combinations on the qiime-191-tt table:
# 3.331e-16 (ordinary float64 noise), confirmed on real NVIDIA and AMD
# hardware -- see
# docs/superpowers/specs/2026-10-02-ssu-unifrac-numba-phase1-design.md.
# 1e-10 gives a huge (~3e5x) safety margin.
GPU_CPU_TOLERANCE = 1e-10

# unweighted_unifrac's normalized=False and variance_adjust=True, and
# weighted_unifrac's variance_adjust=True, are GPU-only in this release (the
# corresponding CPU/numba kernel support was reverted; see CHANGELOG.md), so
# those combinations are checked against the qiime-191-tt SSU-fixture
# distance matrices directly instead of against the CPU numba kernel.
# Measured max abs deviation, GPU vs ssu-ascii-fixture, is far looser than
# the GPU-vs-CPU-numba figure above; see test_unifrac.py's
# SSU_FIXTURE_TOLERANCE for the same figure used on the CPU side.
GPU_FIXTURE_TOLERANCE = 1.5e-6


class GpuBackendTests(TestCase):

    def test_detect_gpu_backend_returns_known_value(self):
        self.assertIn(_detect_gpu_backend(), ("cuda", "hip", None))

    def test_get_cuda_module_raises_cleanly_when_no_gpu(self):
        backend = _detect_gpu_backend()
        if backend is not None:
            self.skipTest("a GPU backend is available in this environment")
        with self.assertRaises(ImportError):
            _get_cuda_module()


class NumbaOptionalImportTests(QiimeTinyTestMixin, TestCase):
    """`_unifrac_gpu` must stay importable, and `engine='gpu'` reachable via
    the array-API fallback, on a system with no numba installed at all.

    numba is an optional dependency project-wide (see
    `skbio.diversity.beta._unifrac.NUMBA_AVAILABLE`), so the module-level
    import of this module cannot assume numba is present -- it must defer
    any numba import to the point where a GPU kernel is actually compiled,
    which only happens after a real GPU backend has already been detected.

    Simulates "numba not installed" the same way
    `skbio.util.tests.test_plotting` simulates "matplotlib not installed":
    setting `sys.modules['numba'] = None` makes any subsequent `import
    numba` (including submodule imports like `from numba import cuda`)
    raise ``ImportError`` without numba actually being uninstalled.

    Executes a throwaway *copy* of ``_unifrac_gpu``'s source (via
    ``importlib.util``) rather than reloading the real, already-imported
    module in place: other test modules hold direct references to the real
    module's mutable globals (e.g. ``_unavailable_backends``,
    ``_KERNEL_CACHE``), and `importlib.reload` would rebind those names to
    new objects in the real module's namespace without updating anyone
    else's already-imported reference to the old one -- silently
    desynchronizing this test's view of that state from everyone else's.
    A separate module object sidesteps that entirely.
    """

    def _load_unifrac_gpu_without_numba(self):
        """Execute a fresh copy of `_unifrac_gpu` with `numba` unimportable.

        Registers restoration of the real `numba` via `addCleanup`, so this
        always runs even if the test body raises.
        """
        import skbio.diversity.beta._unifrac_gpu as gpu_mod

        numba_backup = {
            name: mod for name, mod in sys.modules.items()
            if name == "numba" or name.startswith("numba.")
        }
        for name in numba_backup:
            del sys.modules[name]
        sys.modules["numba"] = None

        def _restore():
            del sys.modules["numba"]
            sys.modules.update(numba_backup)

        self.addCleanup(_restore)

        spec = importlib.util.spec_from_file_location(
            "_unifrac_gpu_numba_absent_test_copy", gpu_mod.__file__
        )
        fresh_mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(fresh_mod)
        return fresh_mod

    def test_import_does_not_require_numba(self):
        # The real regression: a module-level `from numba import float64`
        # would raise ModuleNotFoundError right here, before engine='gpu'
        # ever gets a chance to fall back to the array-API path.
        fresh_mod = self._load_unifrac_gpu_without_numba()
        self.assertIsNone(fresh_mod._detect_gpu_backend())

    def test_gpu_or_xp_dispatch_works_without_numba(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        fresh_mod = self._load_unifrac_gpu_without_numba()
        # _detect_gpu_backend() is None (numba unimportable), so this must use
        # the array-API fallback, not raise.
        obs = fresh_mod._weighted_unifrac_gpu_or_xp(
            table, taxa, tree, normalized=True, validate=True
        )
        self.assertEqual(len(obs), len(table) * (len(table) - 1) // 2)


class BuildPairIndexTests(TestCase):
    """GPU-independent coverage for `_build_pair_index` (no GPU required)."""

    def test_matches_scipy_pdist_condensed_order(self):
        from scipy.spatial.distance import squareform

        n = 5
        pair_i, pair_j = _build_pair_index(n)
        n_pairs = n * (n - 1) // 2
        self.assertEqual(pair_i.shape, (n_pairs,))
        self.assertEqual(pair_j.shape, (n_pairs,))

        # Reconstruct which (i, j) scipy's condensed order expects at each
        # index via squareform, and check _build_pair_index agrees exactly.
        square = np.arange(n * n).reshape(n, n)
        condensed = squareform(square, checks=False)
        for idx in range(n_pairs):
            i, j = pair_i[idx], pair_j[idx]
            self.assertLess(i, j)
            self.assertEqual(condensed[idx], square[i, j])

    def test_dtype_and_trivial_cases(self):
        for n in (0, 1, 2, 3):
            pair_i, pair_j = _build_pair_index(n)
            expected_len = n * (n - 1) // 2
            self.assertEqual(len(pair_i), expected_len)
            self.assertEqual(len(pair_j), expected_len)
            self.assertEqual(pair_i.dtype, np.int32)
            self.assertEqual(pair_j.dtype, np.int32)


class _FakeDeviceArray:
    """Stands in for a `numba.cuda` device array, just enough to support
    `_probe_tile_config_candidate`'s `copy_to_host()` call at the end of a
    (faked) launch."""

    def __init__(self, arr):
        self._arr = arr

    def copy_to_host(self):
        return self._arr


class _ProbingFakeCuda:
    """Stands in for `numba.cuda`/`numba.hip`, good enough to drive
    `_probe_tile_config`'s compile-launch-catch control flow without real
    GPU hardware.

    Each call to `.jit(func)` corresponds to one candidate attempt (since
    `_make_unifrac_kernel` is keyed by `(id(cuda), tile, node_chunk)`, a
    fresh cache miss -- and hence a fresh `.jit()` call -- happens for
    every distinct candidate tried). The returned dispatcher's
    `[grid, block](...)` call raises the next scripted outcome in
    ``launch_outcomes`` instead of actually running the kernel body (which
    uses `cuda.shared.array`/`threadIdx`/`syncthreads`, none of which this
    fake implements) -- the probe only cares whether the launch raises, not
    what the kernel computes.
    """

    def __init__(self, launch_outcomes):
        # One entry per expected `.jit()` call: None means "succeeds",
        # an exception instance means "raise this instead".
        self._outcomes = list(launch_outcomes)
        self.jit_call_count = 0

    def to_device(self, arr):
        return np.asarray(arr)

    def device_array(self, shape, dtype):
        return _FakeDeviceArray(np.zeros(shape, dtype=dtype))

    def jit(self, func):
        outcome_index = self.jit_call_count
        self.jit_call_count += 1
        outcomes = self._outcomes

        class _Dispatcher:
            @staticmethod
            def __getitem__(grid_block):
                def _launch(*args, **kwargs):
                    outcome = outcomes[outcome_index]
                    if outcome is not None:
                        raise outcome

                return _launch

        return _Dispatcher()


class ProbeTileConfigCandidateTests(TestCase):
    """`_probe_tile_config_candidate` drives a real compile+launch of the
    actual kernel against tiny dummy data (no GPU required, via
    `_ProbingFakeCuda`)."""

    def test_succeeds_silently_when_launch_does_not_raise(self):
        fake = _ProbingFakeCuda(launch_outcomes=[None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        # Should not raise.
        _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertEqual(fake.jit_call_count, 1)

    def test_succeeds_silently_for_the_variance_adjust_false_variant_too(self):
        fake = _ProbingFakeCuda(launch_outcomes=[None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        _probe_tile_config_candidate(fake, 16, 32, False)
        self.assertEqual(fake.jit_call_count, 1)

    def test_propagates_the_launch_exception(self):
        fake = _ProbingFakeCuda(launch_outcomes=[RuntimeError("boom")])
        self.addCleanup(_cleanup_kernel_cache, fake)
        with self.assertRaisesRegex(RuntimeError, "boom"):
            _probe_tile_config_candidate(fake, 32, 32, True)


class LaunchUnifracKernelEvictsFailedCacheEntryTests(TestCase):
    """`_make_unifrac_kernel` caches a compiled kernel before it is ever
    launched, so a candidate that compiles fine but fails at launch time
    (the NVIDIA tile-probe failure mode) must not leave a dead entry behind
    in `_KERNEL_CACHE` -- a small per-process leak of a reference to the
    CUDA module and a broken Dispatcher. `_launch_unifrac_kernel` (driven
    here via `_probe_tile_config_candidate`, same as the probe) must evict
    that entry before re-raising.
    """

    def test_failed_launch_does_not_linger_in_kernel_cache(self):
        fake = _ProbingFakeCuda(launch_outcomes=[RuntimeError("boom")])
        self.addCleanup(_cleanup_kernel_cache, fake)
        cache_key = (id(fake), 16, 32, True)
        with self.assertRaisesRegex(RuntimeError, "boom"):
            _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertNotIn(cache_key, _KERNEL_CACHE)

    def test_failed_launch_does_not_linger_in_kernel_cache_variance_adjust_false(self):
        # Same as above, for the variance_adjust=False variant -- the key
        # shape behind the real A10 regression (see `_probe_tile_config`'s
        # docstring): this variant's cache entry must not linger either.
        fake = _ProbingFakeCuda(launch_outcomes=[RuntimeError("boom")])
        self.addCleanup(_cleanup_kernel_cache, fake)
        cache_key = (id(fake), 16, 32, False)
        with self.assertRaisesRegex(RuntimeError, "boom"):
            _probe_tile_config_candidate(fake, 16, 32, False)
        self.assertNotIn(cache_key, _KERNEL_CACHE)

    def test_successful_launch_still_caches_the_kernel(self):
        # Sanity check for the other side of the fix: a candidate that
        # actually works must still end up cached (eviction must only
        # trigger on failure).
        fake = _ProbingFakeCuda(launch_outcomes=[None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        cache_key = (id(fake), 16, 32, True)
        _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertIn(cache_key, _KERNEL_CACHE)

    def test_later_retry_after_a_failed_launch_recompiles_rather_than_reusing(self):
        # With the dead entry evicted, a later attempt at the same
        # (tile, node_chunk) must recompile (a fresh `.jit()` call) instead
        # of returning the broken cached Dispatcher.
        fake = _ProbingFakeCuda(launch_outcomes=[RuntimeError("boom"), None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        with self.assertRaisesRegex(RuntimeError, "boom"):
            _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertEqual(fake.jit_call_count, 1)
        _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertEqual(fake.jit_call_count, 2)


class ProbeTileConfigTests(TestCase):
    """`_probe_tile_config` walks `_TILE_CANDIDATES` in order and returns
    the first one that actually works, handling both known failure shapes
    (an NVIDIA-style launch-time error, an AMD-style compile-time
    RuntimeError) plus anything else via a broad catch, since these are
    controlled probes against synthetic data, not user input.
    """

    def test_first_candidate_succeeding_is_used_without_trying_others(self):
        # Both the variance_adjust=True and =False probes of the first
        # candidate must succeed -- 2 jit() calls, not 1 (see
        # `_probe_tile_config`'s docstring for why both are probed).
        fake = _ProbingFakeCuda(launch_outcomes=[None, None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "hip")
        self.assertEqual(result, _TILE_CANDIDATES[0])
        self.assertEqual(fake.jit_call_count, 2)

    def test_nvidia_style_launch_error_falls_through_to_next_candidate(self):
        # Simulates the real NVIDIA failure mode: the most aggressive
        # candidate(s) raise a launch-time resource error on their
        # variance_adjust=True probe, and the probe moves on to try the
        # next candidate instead of giving up.
        cuda_error = type("CudaAPIError", (Exception,), {})(
            "CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES"
        )
        fake = _ProbingFakeCuda(launch_outcomes=[cuda_error, cuda_error, None, None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "cuda")
        self.assertEqual(result, _TILE_CANDIDATES[2])
        self.assertEqual(fake.jit_call_count, 4)

    def test_amd_style_compile_error_falls_through_to_next_candidate(self):
        # Simulates the real AMD failure mode: a compile-time RuntimeError
        # from rocm.amd_comgr (numba.hip compiles lazily, on first
        # invocation, so this still surfaces at "launch" time from this
        # function's point of view).
        comgr_error = RuntimeError(
            "AMD_COMGR_ACTION_CODEGEN_BC_TO_RELOCATABLE failed: "
            "local memory (66560) exceeds limit (65536)"
        )
        fake = _ProbingFakeCuda(launch_outcomes=[comgr_error, None, None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "hip")
        self.assertEqual(result, _TILE_CANDIDATES[1])
        self.assertEqual(fake.jit_call_count, 3)

    def test_mixed_failure_types_are_all_caught(self):
        # A third, unanticipated exception type should also be caught by
        # the deliberately-broad probe, not just the two known shapes.
        outcomes = [
            type("CudaAPIError", (Exception,), {})("launch failed"),
            RuntimeError("AMD_COMGR_ACTION_CODEGEN_BC_TO_RELOCATABLE failed"),
            ValueError("some other, unanticipated failure"),
            None,
            None,
        ]
        fake = _ProbingFakeCuda(launch_outcomes=outcomes)
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "cuda")
        self.assertEqual(result, _TILE_CANDIDATES[3])
        self.assertEqual(fake.jit_call_count, 5)

    def test_raises_clear_error_when_every_candidate_fails(self):
        fail = RuntimeError("nope")
        fake = _ProbingFakeCuda(launch_outcomes=[fail] * len(_TILE_CANDIDATES))
        self.addCleanup(_cleanup_kernel_cache, fake)
        with self.assertRaisesRegex(RuntimeError, "No safe"):
            _probe_tile_config(fake, "cuda")
        self.assertEqual(fake.jit_call_count, len(_TILE_CANDIDATES))

    def test_variance_adjust_false_failing_rejects_a_candidate_whose_true_variant_passed(
        self,
    ):
        # Regression for the real A10 CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES
        # crash: variance_adjust is a compile-time constant (see A2), so its
        # True/False kernel variants are two separately-compiled artifacts
        # whose register footprints are not reliably ordered. A tile config
        # whose variance_adjust=True probe passes but whose
        # variance_adjust=False probe fails must still be rejected, not
        # cached and handed to real variance_adjust=False calls later.
        fail = RuntimeError("CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES")
        fake = _ProbingFakeCuda(launch_outcomes=[None, fail, None, None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "cuda")
        self.assertEqual(result, _TILE_CANDIDATES[1])
        self.assertEqual(fake.jit_call_count, 4)


class GetTileConfigTests(TestCase):
    """`_get_tile_config` probes once per backend and caches the result,
    rather than looking up a fixed per-backend dict (the hardware-specific
    safe config cannot be known statically, see `_unifrac_gpu.py`)."""

    def setUp(self):
        # Each test starts from a clean cache so probes/call-counts are
        # deterministic and tests cannot see each other's cached results.
        self._orig_cache = dict(_TILE_CONFIG_CACHE)
        _TILE_CONFIG_CACHE.clear()
        self.addCleanup(self._restore_cache)

    def _restore_cache(self):
        _TILE_CONFIG_CACHE.clear()
        _TILE_CONFIG_CACHE.update(self._orig_cache)

    def test_probes_and_caches_on_first_call(self):
        with patch(
            "skbio.diversity.beta._unifrac_gpu._probe_tile_config",
            return_value=(16, 32),
        ) as mock_probe:
            result = _get_tile_config("fake-cuda-module", "cuda")
        self.assertEqual(result, (16, 32))
        mock_probe.assert_called_once_with("fake-cuda-module", "cuda")
        self.assertEqual(_TILE_CONFIG_CACHE["cuda"], (16, 32))

    def test_second_call_for_the_same_backend_does_not_reprobe(self):
        with patch(
            "skbio.diversity.beta._unifrac_gpu._probe_tile_config",
            return_value=(32, 32),
        ) as mock_probe:
            first = _get_tile_config("fake-cuda-module", "hip")
            second = _get_tile_config("fake-cuda-module", "hip")
        self.assertEqual(first, (32, 32))
        self.assertEqual(second, (32, 32))
        mock_probe.assert_called_once()

    def test_distinct_backends_are_probed_and_cached_independently(self):
        def fake_probe(cuda, backend):
            return {"cuda": (16, 32), "hip": (32, 32)}[backend]

        with patch(
            "skbio.diversity.beta._unifrac_gpu._probe_tile_config",
            side_effect=fake_probe,
        ) as mock_probe:
            cuda_result = _get_tile_config("fake-cuda-module", "cuda")
            hip_result = _get_tile_config("fake-hip-module", "hip")
            # Calling again for either backend must not re-probe.
            _get_tile_config("fake-cuda-module", "cuda")
            _get_tile_config("fake-hip-module", "hip")
        self.assertEqual(cuda_result, (16, 32))
        self.assertEqual(hip_result, (32, 32))
        self.assertEqual(mock_probe.call_count, 2)

    def test_a_failed_probe_is_not_cached_and_raises(self):
        with patch(
            "skbio.diversity.beta._unifrac_gpu._probe_tile_config",
            side_effect=RuntimeError("No safe tile configuration"),
        ) as mock_probe:
            with self.assertRaises(RuntimeError):
                _get_tile_config("fake-cuda-module", "cuda")
        self.assertNotIn("cuda", _TILE_CONFIG_CACHE)
        mock_probe.assert_called_once()


class MakeUnifracKernelTests(TestCase):
    """Kernel compilation is memoized per backend (no GPU required)."""

    def _fake_cuda_module(self):
        class FakeCuda:
            """Stands in for numba.cuda; `jit` just returns the function."""

            @staticmethod
            def jit(func):
                return func

        return FakeCuda()

    def test_kernel_is_compiled_once_per_backend(self):
        fake = self._fake_cuda_module()
        cache_key = (id(fake), 32, 32, True)
        self.addCleanup(_KERNEL_CACHE.pop, cache_key, None)
        first = _make_unifrac_kernel(fake, 32, 32, True)
        second = _make_unifrac_kernel(fake, 32, 32, True)
        self.assertIs(first, second)

    def test_distinct_backends_get_distinct_kernels(self):
        fake1 = self._fake_cuda_module()
        fake2 = self._fake_cuda_module()
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake1), 32, 32, True), None)
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake2), 32, 32, True), None)
        self.assertIsNot(
            _make_unifrac_kernel(fake1, 32, 32, True),
            _make_unifrac_kernel(fake2, 32, 32, True),
        )

    def test_distinct_tile_configs_get_distinct_kernels(self):
        fake = self._fake_cuda_module()
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 32, 32, True), None)
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 16, 32, True), None)
        self.assertIsNot(
            _make_unifrac_kernel(fake, 32, 32, True),
            _make_unifrac_kernel(fake, 16, 32, True),
        )

    def test_distinct_variance_adjust_values_get_distinct_kernels(self):
        # A2: variance_adjust is now closed over as a compile-time constant
        # (2 shared tiles declared instead of 4 when False) rather than
        # passed as a runtime kernel argument, so it is its own axis of the
        # cache key, same as tile/node_chunk.
        fake = self._fake_cuda_module()
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 32, 32, True), None)
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 32, 32, False), None)
        self.assertIsNot(
            _make_unifrac_kernel(fake, 32, 32, True),
            _make_unifrac_kernel(fake, 32, 32, False),
        )

    def test_same_variance_adjust_value_is_memoized(self):
        fake = self._fake_cuda_module()
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 32, 32, False), None)
        first = _make_unifrac_kernel(fake, 32, 32, False)
        second = _make_unifrac_kernel(fake, 32, 32, False)
        self.assertIs(first, second)


class _GpuRequiredTestCase(QiimeTinyTestMixin, TestCase):
    """Base for the fused-kernel correctness tests: needs real GPU hardware.

    Everything else in this module is exercised without a GPU, via fake cuda
    modules; the classes below launch the real kernel, so they skip when no
    backend is detected.
    """

    def setUp(self):
        if _detect_gpu_backend() is None:
            self.skipTest("no GPU backend available")


class WeightedUnifracGpuTests(_GpuRequiredTestCase):

    def test_weighted_unifrac_gpu_matches_cpu_unnormalized(self):
        from skbio.diversity.beta._unifrac import _weighted_unifrac_pdist_numba
        table, taxa, tree, _ = self._load_qiime_191_tt()
        gpu = _weighted_unifrac_gpu(
            table, taxa, tree,
            normalized=False, variance_adjust=False, validate=True,
        )
        cpu = _weighted_unifrac_pdist_numba(
            table, taxa, tree, normalized=False, validate=True
        )
        np.testing.assert_allclose(gpu, cpu, rtol=0, atol=GPU_CPU_TOLERANCE)

    def test_weighted_unifrac_gpu_matches_cpu_normalized(self):
        from skbio.diversity.beta._unifrac import _weighted_unifrac_pdist_numba
        table, taxa, tree, _ = self._load_qiime_191_tt()
        gpu = _weighted_unifrac_gpu(
            table, taxa, tree,
            normalized=True, variance_adjust=False, validate=True,
        )
        cpu = _weighted_unifrac_pdist_numba(
            table, taxa, tree, normalized=True, validate=True
        )
        np.testing.assert_allclose(gpu, cpu, rtol=0, atol=GPU_CPU_TOLERANCE)

    def test_weighted_unifrac_gpu_matches_fixture_variance_adjusted_unnormalized(self):
        # weighted_unifrac's variance_adjust is GPU-only in this release (the
        # CPU/numba kernel's variance_adjust support was reverted; see
        # CHANGELOG.md), so this checks the GPU result against the
        # qiime-191-tt SSU fixture directly rather than against the CPU
        # kernel.
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        gpu = _weighted_unifrac_gpu(
            table, taxa, tree,
            normalized=False, variance_adjust=True, validate=True,
        )
        self._assert_matches_dm_fixture(
            gpu, sample_ids, 'weighted_unifrac_vaw_dm.txt',
            GPU_FIXTURE_TOLERANCE)

    def test_weighted_unifrac_gpu_matches_fixture_variance_adjusted_normalized(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        gpu = _weighted_unifrac_gpu(
            table, taxa, tree,
            normalized=True, variance_adjust=True, validate=True,
        )
        self._assert_matches_dm_fixture(
            gpu, sample_ids, 'weighted_normalized_unifrac_vaw_dm.txt',
            GPU_FIXTURE_TOLERANCE)


class UnweightedUnifracGpuTests(_GpuRequiredTestCase):
    # unweighted_unifrac's normalized and variance_adjust kwargs are GPU-only
    # in this release (the CPU/numba kernel was reverted to its pre-PR,
    # always-normalized, no-variance_adjust behavior; see CHANGELOG.md), so
    # all four combinations here are checked against the qiime-191-tt SSU
    # fixture distance matrices directly rather than against the CPU kernel.

    def test_unweighted_unifrac_gpu_matches_fixture_unnormalized(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        gpu = _unweighted_unifrac_gpu(
            table, taxa, tree,
            normalized=False, variance_adjust=False, validate=True,
        )
        self._assert_matches_dm_fixture(
            gpu, sample_ids, 'unweighted_unnormalized_unifrac_dm.txt',
            GPU_FIXTURE_TOLERANCE)

    def test_unweighted_unifrac_gpu_matches_fixture_normalized(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        gpu = _unweighted_unifrac_gpu(
            table, taxa, tree,
            normalized=True, variance_adjust=False, validate=True,
        )
        self._assert_matches_dm_fixture(
            gpu, sample_ids, 'unweighted_unifrac_dm.txt',
            GPU_FIXTURE_TOLERANCE)

    def test_unweighted_unifrac_gpu_matches_fixture_variance_adjusted_unnormalized(
        self,
    ):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        gpu = _unweighted_unifrac_gpu(
            table, taxa, tree,
            normalized=False, variance_adjust=True, validate=True,
        )
        self._assert_matches_dm_fixture(
            gpu, sample_ids, 'unweighted_unnormalized_unifrac_vaw_dm.txt',
            GPU_FIXTURE_TOLERANCE)

    def test_unweighted_unifrac_gpu_matches_fixture_variance_adjusted_normalized(
        self,
    ):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        gpu = _unweighted_unifrac_gpu(
            table, taxa, tree,
            normalized=True, variance_adjust=True, validate=True,
        )
        self._assert_matches_dm_fixture(
            gpu, sample_ids, 'unweighted_unifrac_vaw_dm.txt',
            GPU_FIXTURE_TOLERANCE)


class GeneralizedUnifracGpuTests(_GpuRequiredTestCase):

    def test_generalized_unifrac_gpu_matches_cpu_alpha_half(self):
        from skbio.diversity.beta._unifrac import _generalized_unifrac_pdist_numba
        table, taxa, tree, _ = self._load_qiime_191_tt()
        gpu = _generalized_unifrac_gpu(
            table, taxa, tree, alpha=0.5, variance_adjust=False, validate=True,
        )
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=0.5, variance_adjust=False, validate=True
        )
        np.testing.assert_allclose(gpu, cpu, rtol=0, atol=GPU_CPU_TOLERANCE)

    def test_generalized_unifrac_gpu_matches_cpu_alpha_one(self):
        from skbio.diversity.beta._unifrac import _generalized_unifrac_pdist_numba
        table, taxa, tree, _ = self._load_qiime_191_tt()
        gpu = _generalized_unifrac_gpu(
            table, taxa, tree, alpha=1.0, variance_adjust=False, validate=True,
        )
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=1.0, variance_adjust=False, validate=True
        )
        np.testing.assert_allclose(gpu, cpu, rtol=0, atol=GPU_CPU_TOLERANCE)

    def test_generalized_unifrac_gpu_matches_cpu_variance_adjusted(self):
        from skbio.diversity.beta._unifrac import _generalized_unifrac_pdist_numba
        table, taxa, tree, _ = self._load_qiime_191_tt()
        gpu = _generalized_unifrac_gpu(
            table, taxa, tree, alpha=0.5, variance_adjust=True, validate=True,
        )
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=0.5, variance_adjust=True, validate=True
        )
        np.testing.assert_allclose(gpu, cpu, rtol=0, atol=GPU_CPU_TOLERANCE)

    def test_generalized_unifrac_gpu_matches_cpu_variance_adjusted_alpha_one(self):
        from skbio.diversity.beta._unifrac import _generalized_unifrac_pdist_numba
        table, taxa, tree, _ = self._load_qiime_191_tt()
        gpu = _generalized_unifrac_gpu(
            table, taxa, tree, alpha=1.0, variance_adjust=True, validate=True,
        )
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=1.0, variance_adjust=True, validate=True
        )
        np.testing.assert_allclose(gpu, cpu, rtol=0, atol=GPU_CPU_TOLERANCE)


if __name__ == "__main__":
    main()
