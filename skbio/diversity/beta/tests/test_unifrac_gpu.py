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
from skbio.util import numba_code


def _cleanup_kernel_cache(fake_cuda):
    """Drop every `_KERNEL_CACHE` entry a fake cuda module may have created.

    Keys embed ``id(module)``, which a later object can reuse once this one is
    collected, so the entries must not outlive the test that made them.
    """
    for tile, node_chunk in _TILE_CANDIDATES:
        for variance_adjust in (True, False):
            _KERNEL_CACHE.pop(
                (id(fake_cuda), tile, node_chunk, variance_adjust), None
            )


# Max abs deviation, GPU vs CPU-numba, on the qiime-191-tt table is 3.3e-16
# (ordinary float64 noise) on both NVIDIA and AMD hardware.
GPU_CPU_TOLERANCE = 1e-10

# unweighted_unifrac's normalized=False and variance_adjust=True, and
# weighted_unifrac's variance_adjust=True, have no CPU/numba counterpart in
# this release (see CHANGELOG.md), so they are compared against the
# qiime-191-tt SSU fixture matrices, which agree only to ~1e-6. Matches
# test_unifrac.py's SSU_FIXTURE_TOLERANCE.
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
    the array-API fallback, with no numba installed (it is optional
    project-wide).

    ``sys.modules['numba'] = None`` makes any subsequent ``import numba``
    raise ``ImportError`` without uninstalling numba. A throwaway *copy* of
    the module is executed rather than `importlib.reload`, which would rebind
    the real module's mutable globals (``_unavailable_backends``,
    ``_KERNEL_CACHE``) that other test modules already hold references to.
    """

    def _load_unifrac_gpu_without_numba(self):
        """Execute a fresh copy of `_unifrac_gpu` with `numba` unimportable."""
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
        # A module-level `from numba import float64` would raise right here.
        fresh_mod = self._load_unifrac_gpu_without_numba()
        self.assertIsNone(fresh_mod._detect_gpu_backend())

    def test_gpu_or_xp_dispatch_works_without_numba(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        fresh_mod = self._load_unifrac_gpu_without_numba()
        # _detect_gpu_backend() is None here, so this takes the array-API path.
        obs = fresh_mod._weighted_unifrac_gpu_or_xp(
            table, taxa, tree, normalized=True, validate=True
        )
        self.assertEqual(len(obs), len(table) * (len(table) - 1) // 2)


class _FakeDeviceArray:
    """Stands in for a `numba.cuda` device array (just `copy_to_host()`)."""

    def __init__(self, arr):
        self._arr = arr

    def copy_to_host(self):
        return self._arr


class _ProbingFakeCuda:
    """Stands in for `numba.cuda`/`numba.hip`, enough to drive
    `_probe_tile_config`'s compile-launch-catch flow without GPU hardware.

    One `.jit(func)` call per candidate attempt (the kernel cache is keyed by
    ``id(cuda)`` plus the tile config). The dispatcher's `[grid, block](...)`
    call raises the next scripted outcome instead of running the kernel body,
    which needs `cuda.shared.array`/`threadIdx`/`syncthreads`; the probe only
    cares whether a launch raises.
    """

    def __init__(self, launch_outcomes):
        # One entry per expected `.jit()` call: None succeeds, an exception
        # instance is raised at launch.
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
    """`_probe_tile_config_candidate` compiles and launches the real kernel
    against tiny dummy data (no GPU required, via `_ProbingFakeCuda`)."""

    @numba_code
    def test_succeeds_silently_when_launch_does_not_raise(self):
        fake = _ProbingFakeCuda(launch_outcomes=[None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertEqual(fake.jit_call_count, 1)

    @numba_code
    def test_succeeds_silently_for_the_variance_adjust_false_variant_too(self):
        fake = _ProbingFakeCuda(launch_outcomes=[None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        _probe_tile_config_candidate(fake, 16, 32, False)
        self.assertEqual(fake.jit_call_count, 1)

    @numba_code
    def test_propagates_the_launch_exception(self):
        fake = _ProbingFakeCuda(launch_outcomes=[RuntimeError("boom")])
        self.addCleanup(_cleanup_kernel_cache, fake)
        with self.assertRaisesRegex(RuntimeError, "boom"):
            _probe_tile_config_candidate(fake, 32, 32, True)


class LaunchUnifracKernelEvictsFailedCacheEntryTests(TestCase):
    """`_make_unifrac_kernel` caches a kernel before it is ever launched, so a
    candidate that compiles but fails at launch (the NVIDIA tile-probe failure
    mode) would otherwise leave a broken Dispatcher in `_KERNEL_CACHE`.
    """

    @numba_code
    def test_failed_launch_does_not_linger_in_kernel_cache(self):
        fake = _ProbingFakeCuda(launch_outcomes=[RuntimeError("boom")])
        self.addCleanup(_cleanup_kernel_cache, fake)
        cache_key = (id(fake), 16, 32, True)
        with self.assertRaisesRegex(RuntimeError, "boom"):
            _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertNotIn(cache_key, _KERNEL_CACHE)

    @numba_code
    def test_failed_launch_does_not_linger_in_kernel_cache_variance_adjust_false(self):
        fake = _ProbingFakeCuda(launch_outcomes=[RuntimeError("boom")])
        self.addCleanup(_cleanup_kernel_cache, fake)
        cache_key = (id(fake), 16, 32, False)
        with self.assertRaisesRegex(RuntimeError, "boom"):
            _probe_tile_config_candidate(fake, 16, 32, False)
        self.assertNotIn(cache_key, _KERNEL_CACHE)

    @numba_code
    def test_successful_launch_still_caches_the_kernel(self):
        fake = _ProbingFakeCuda(launch_outcomes=[None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        cache_key = (id(fake), 16, 32, True)
        _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertIn(cache_key, _KERNEL_CACHE)

    @numba_code
    def test_later_retry_after_a_failed_launch_recompiles_rather_than_reusing(self):
        fake = _ProbingFakeCuda(launch_outcomes=[RuntimeError("boom"), None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        with self.assertRaisesRegex(RuntimeError, "boom"):
            _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertEqual(fake.jit_call_count, 1)
        _probe_tile_config_candidate(fake, 16, 32, True)
        self.assertEqual(fake.jit_call_count, 2)


class ProbeTileConfigTests(TestCase):
    """`_probe_tile_config` walks `_TILE_CANDIDATES` in order and returns the
    first that works, catching any exception shape: these are controlled
    probes against synthetic data, not user input.
    """

    @numba_code
    def test_first_candidate_succeeding_is_used_without_trying_others(self):
        # Each candidate is probed for both variance_adjust variants, so a
        # single accepted candidate means 2 jit() calls, not 1.
        fake = _ProbingFakeCuda(launch_outcomes=[None, None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "hip")
        self.assertEqual(result, _TILE_CANDIDATES[0])
        self.assertEqual(fake.jit_call_count, 2)

    @numba_code
    def test_nvidia_style_launch_error_falls_through_to_next_candidate(self):
        cuda_error = type("CudaAPIError", (Exception,), {})(
            "CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES"
        )
        fake = _ProbingFakeCuda(launch_outcomes=[cuda_error, cuda_error, None, None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "cuda")
        self.assertEqual(result, _TILE_CANDIDATES[2])
        self.assertEqual(fake.jit_call_count, 4)

    @numba_code
    def test_amd_style_compile_error_falls_through_to_next_candidate(self):
        # numba.hip compiles lazily, so a rocm.amd_comgr compile error still
        # surfaces at launch time from the probe's point of view.
        comgr_error = RuntimeError(
            "AMD_COMGR_ACTION_CODEGEN_BC_TO_RELOCATABLE failed: "
            "local memory (66560) exceeds limit (65536)"
        )
        fake = _ProbingFakeCuda(launch_outcomes=[comgr_error, None, None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "hip")
        self.assertEqual(result, _TILE_CANDIDATES[1])
        self.assertEqual(fake.jit_call_count, 3)

    @numba_code
    def test_mixed_failure_types_are_all_caught(self):
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

    @numba_code
    def test_raises_clear_error_when_every_candidate_fails(self):
        fail = RuntimeError("nope")
        fake = _ProbingFakeCuda(launch_outcomes=[fail] * len(_TILE_CANDIDATES))
        self.addCleanup(_cleanup_kernel_cache, fake)
        with self.assertRaisesRegex(RuntimeError, "No safe"):
            _probe_tile_config(fake, "cuda")
        self.assertEqual(fake.jit_call_count, len(_TILE_CANDIDATES))

    @numba_code
    def test_variance_adjust_false_failing_rejects_a_candidate_whose_true_variant_passed(
        self,
    ):
        # variance_adjust is a compile-time constant, so its True/False
        # variants are separately compiled artifacts whose register footprints
        # are not reliably ordered: either one failing rejects the candidate.
        fail = RuntimeError("CUDA_ERROR_LAUNCH_OUT_OF_RESOURCES")
        fake = _ProbingFakeCuda(launch_outcomes=[None, fail, None, None])
        self.addCleanup(_cleanup_kernel_cache, fake)
        result = _probe_tile_config(fake, "cuda")
        self.assertEqual(result, _TILE_CANDIDATES[1])
        self.assertEqual(fake.jit_call_count, 4)


class GetTileConfigTests(TestCase):
    """`_get_tile_config` probes once per backend and caches the result; the
    safe config is hardware-specific and cannot be known statically."""

    def setUp(self):
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

    @numba_code
    def test_kernel_is_compiled_once_per_backend(self):
        fake = self._fake_cuda_module()
        cache_key = (id(fake), 32, 32, True)
        self.addCleanup(_KERNEL_CACHE.pop, cache_key, None)
        first = _make_unifrac_kernel(fake, 32, 32, True)
        second = _make_unifrac_kernel(fake, 32, 32, True)
        self.assertIs(first, second)

    @numba_code
    def test_distinct_backends_get_distinct_kernels(self):
        fake1 = self._fake_cuda_module()
        fake2 = self._fake_cuda_module()
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake1), 32, 32, True), None)
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake2), 32, 32, True), None)
        self.assertIsNot(
            _make_unifrac_kernel(fake1, 32, 32, True),
            _make_unifrac_kernel(fake2, 32, 32, True),
        )

    @numba_code
    def test_distinct_tile_configs_get_distinct_kernels(self):
        fake = self._fake_cuda_module()
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 32, 32, True), None)
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 16, 32, True), None)
        self.assertIsNot(
            _make_unifrac_kernel(fake, 32, 32, True),
            _make_unifrac_kernel(fake, 16, 32, True),
        )

    @numba_code
    def test_distinct_variance_adjust_values_get_distinct_kernels(self):
        # variance_adjust is closed over as a compile-time constant, not a
        # runtime argument, so it is its own axis of the cache key.
        fake = self._fake_cuda_module()
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 32, 32, True), None)
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 32, 32, False), None)
        self.assertIsNot(
            _make_unifrac_kernel(fake, 32, 32, True),
            _make_unifrac_kernel(fake, 32, 32, False),
        )

    @numba_code
    def test_same_variance_adjust_value_is_memoized(self):
        fake = self._fake_cuda_module()
        self.addCleanup(_KERNEL_CACHE.pop, (id(fake), 32, 32, False), None)
        first = _make_unifrac_kernel(fake, 32, 32, False)
        second = _make_unifrac_kernel(fake, 32, 32, False)
        self.assertIs(first, second)


class _GpuRequiredTestCase(QiimeTinyTestMixin, TestCase):
    """Base for the fused-kernel correctness tests, which launch the real
    kernel and so skip when no GPU backend is detected."""

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
    # The CPU/numba kernel is always-normalized with no variance_adjust, so
    # all four combinations here are checked against the SSU fixtures.

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
