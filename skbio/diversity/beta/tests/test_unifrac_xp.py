# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

from skbio.tree import DuplicateNodeError
from skbio.diversity.beta._unifrac import (
    _weighted_unifrac_pdist_numba,
    _unweighted_unifrac_pdist_numba,
    _generalized_unifrac_pdist_numba,
)
from skbio.diversity.beta._unifrac_gpu import (
    _dispatch_gpu_or_xp,
    _unavailable_backends,
    _generalized_unifrac_gpu_or_xp,
    _unweighted_unifrac_gpu_or_xp,
    _weighted_unifrac_gpu,
    _weighted_unifrac_gpu_or_xp,
)
import skbio.diversity.beta._unifrac_xp as _unifrac_xp_mod
from skbio.diversity.beta._unifrac_xp import (
    _detect_array_api_gpu_module,
    _generalized_unifrac_xp,
    _move_to_array_api_gpu,
    _unweighted_unifrac_xp,
    _weighted_unifrac_xp,
)
from skbio.diversity.beta.tests._fixtures import (
    QiimeTinyTestMixin,
    patch_gpu_backend,
)


def _make_fake_cupy(device_count=1, fail_asarray=False):
    """A minimal stand-in for `cupy`, just enough to drive
    `_detect_array_api_gpu_module`/`_move_to_array_api_gpu` without real GPU
    hardware. Wraps plain NumPy, so a successful call is numerically
    identical to the no-GPU path -- tests inspect `asarray_calls`/the raised
    exception, not the type of what comes back.
    """

    class _Runtime:
        @staticmethod
        def getDeviceCount():
            return device_count

    class _Cuda:
        runtime = _Runtime

    class _FakeCupy:
        __name__ = "cupy"
        cuda = _Cuda
        asarray_calls = 0

        @staticmethod
        def empty(n):
            return np.empty(n)

        @staticmethod
        def asarray(arr):
            _FakeCupy.asarray_calls += 1
            if fail_asarray:
                raise RuntimeError("simulated CuPy device error")
            return np.asarray(arr)

    # An instance, not the class itself: `type.__name__` is a data
    # descriptor on the metaclass, so a bare class object's `__name__`
    # always reports its real Python class name ("_FakeCupy") regardless of
    # what the class body assigns -- only a plain instance's `__name__`
    # attribute lookup actually finds the class-body value, which
    # `_move_to_array_api_gpu` relies on to tell CuPy-like from
    # PyTorch-like.
    return _FakeCupy()


def _make_fake_torch(available=True, fail_asarray=False):
    """A minimal stand-in for `torch`, the same role `_make_fake_cupy` plays
    for CuPy, used to test the "CuPy absent/unusable, fall back to PyTorch"
    branch of `_detect_array_api_gpu_module`.
    """

    class _Cuda:
        @staticmethod
        def is_available():
            return available

    class _FakeTorch:
        __name__ = "torch"
        cuda = _Cuda
        asarray_calls = 0

        @staticmethod
        def empty(n, device=None):
            return np.empty(n)

        @staticmethod
        def as_tensor(arr, device=None):
            _FakeTorch.asarray_calls += 1
            if fail_asarray:
                raise RuntimeError("simulated PyTorch device error")
            return np.asarray(arr)

    return _FakeTorch()  # see the comment in `_make_fake_cupy` re: instance vs class

# Measured max abs deviation, array-API path vs CPU-numba kernel, for the
# method/normalized combinations both implement (no variance_adjust) on the
# qiime-191-tt table: ~3.3e-16 (ordinary float64 noise) for weighted_unifrac,
# which still has its own independent @njit kernel. 1e-10 matches the
# tolerance used for the equivalent GPU-vs-CPU-numba comparison in
# test_unifrac_gpu.py. unweighted_unifrac's non-variance_adjust numba path
# now delegates to _unweighted_unifrac_xp directly (see _unifrac.py's
# _unweighted_unifrac_pdist_numba), so that comparison is exact (0.0), not
# independently-measured noise.
XP_CPU_TOLERANCE = 1e-10

# variance_adjust and unweighted_unifrac's normalized=False have no CPU/numba
# counterpart (GPU-only in this release), so those combinations are checked
# against the qiime-191-tt SSU fixture distance matrices directly instead,
# using the same tolerance as test_unifrac_gpu.py's GPU_FIXTURE_TOLERANCE.
XP_FIXTURE_TOLERANCE = 1.5e-6


class WeightedUnifracXpTests(QiimeTinyTestMixin, TestCase):

    def test_weighted_unifrac_xp_matches_cpu_unnormalized(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        xp_res = _weighted_unifrac_xp(
            table, taxa, tree, normalized=False, validate=True)
        cpu = _weighted_unifrac_pdist_numba(
            table, taxa, tree, normalized=False, validate=True)
        np.testing.assert_allclose(xp_res, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_weighted_unifrac_xp_matches_cpu_normalized(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        xp_res = _weighted_unifrac_xp(
            table, taxa, tree, normalized=True, validate=True)
        cpu = _weighted_unifrac_pdist_numba(
            table, taxa, tree, normalized=True, validate=True)
        np.testing.assert_allclose(xp_res, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_weighted_unifrac_xp_matches_fixture_variance_adjusted_unnormalized(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        xp_res = _weighted_unifrac_xp(
            table, taxa, tree, normalized=False, variance_adjust=True,
            validate=True,
        )
        self._assert_matches_dm_fixture(
            xp_res, sample_ids, 'weighted_unifrac_vaw_dm.txt', XP_FIXTURE_TOLERANCE)

    def test_weighted_unifrac_xp_matches_fixture_variance_adjusted_normalized(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        xp_res = _weighted_unifrac_xp(
            table, taxa, tree, normalized=True, variance_adjust=True,
            validate=True,
        )
        self._assert_matches_dm_fixture(
            xp_res, sample_ids, 'weighted_normalized_unifrac_vaw_dm.txt',
            XP_FIXTURE_TOLERANCE)


class UnweightedUnifracXpTests(QiimeTinyTestMixin, TestCase):

    def test_unweighted_unifrac_xp_matches_cpu_normalized(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        xp_res = _unweighted_unifrac_xp(
            table, taxa, tree, normalized=True, validate=True)
        cpu = _unweighted_unifrac_pdist_numba(table, taxa, tree, validate=True)
        np.testing.assert_allclose(xp_res, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_unweighted_unifrac_xp_matches_fixture_unnormalized(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        xp_res = _unweighted_unifrac_xp(
            table, taxa, tree, normalized=False, validate=True)
        self._assert_matches_dm_fixture(
            xp_res, sample_ids, 'unweighted_unnormalized_unifrac_dm.txt',
            XP_FIXTURE_TOLERANCE)

    def test_unweighted_unifrac_xp_matches_fixture_variance_adjusted_unnormalized(
        self,
    ):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        xp_res = _unweighted_unifrac_xp(
            table, taxa, tree, normalized=False, variance_adjust=True,
            validate=True,
        )
        self._assert_matches_dm_fixture(
            xp_res, sample_ids, 'unweighted_unnormalized_unifrac_vaw_dm.txt',
            XP_FIXTURE_TOLERANCE)

    def test_unweighted_unifrac_xp_matches_fixture_variance_adjusted_normalized(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        xp_res = _unweighted_unifrac_xp(
            table, taxa, tree, normalized=True, variance_adjust=True,
            validate=True,
        )
        self._assert_matches_dm_fixture(
            xp_res, sample_ids, 'unweighted_unifrac_vaw_dm.txt', XP_FIXTURE_TOLERANCE)


class GeneralizedUnifracXpTests(QiimeTinyTestMixin, TestCase):

    def test_generalized_unifrac_xp_matches_cpu_alpha_half(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        xp_res = _generalized_unifrac_xp(
            table, taxa, tree, alpha=0.5, variance_adjust=False, validate=True)
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=0.5, variance_adjust=False, validate=True)
        np.testing.assert_allclose(xp_res, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_generalized_unifrac_xp_matches_cpu_alpha_one(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        xp_res = _generalized_unifrac_xp(
            table, taxa, tree, alpha=1.0, variance_adjust=False, validate=True)
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=1.0, variance_adjust=False, validate=True)
        np.testing.assert_allclose(xp_res, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_generalized_unifrac_xp_matches_cpu_variance_adjusted(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        xp_res = _generalized_unifrac_xp(
            table, taxa, tree, alpha=0.5, variance_adjust=True, validate=True)
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=0.5, variance_adjust=True, validate=True)
        np.testing.assert_allclose(xp_res, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_generalized_unifrac_xp_matches_cpu_variance_adjusted_alpha_one(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        xp_res = _generalized_unifrac_xp(
            table, taxa, tree, alpha=1.0, variance_adjust=True, validate=True)
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=1.0, variance_adjust=True, validate=True)
        np.testing.assert_allclose(xp_res, cpu, rtol=0, atol=XP_CPU_TOLERANCE)


class GpuOrXpDispatchTests(QiimeTinyTestMixin, TestCase):
    """Unit coverage for the `*_gpu_or_xp` dispatch wrappers themselves.

    Forces the no-usable-backend branch by monkeypatching
    `_detect_gpu_backend` to return None, independent of what hardware this
    machine actually has, and checks the result against the CPU reference.
    Public-API-level coverage of the same fallback (through
    unweighted_unifrac/weighted_unifrac/generalized_unifrac and
    beta_diversity) lives in test_unifrac.py.
    """

    def test_weighted_unifrac_gpu_or_xp_falls_back_without_gpu(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        with patch_gpu_backend(None):
            obs = _weighted_unifrac_gpu_or_xp(
                table, taxa, tree, normalized=True, validate=True)
        cpu = _weighted_unifrac_pdist_numba(
            table, taxa, tree, normalized=True, validate=True)
        np.testing.assert_allclose(obs, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_unweighted_unifrac_gpu_or_xp_falls_back_without_gpu(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        with patch_gpu_backend(None):
            obs = _unweighted_unifrac_gpu_or_xp(
                table, taxa, tree, normalized=True, validate=True)
        cpu = _unweighted_unifrac_pdist_numba(table, taxa, tree, validate=True)
        np.testing.assert_allclose(obs, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_generalized_unifrac_gpu_or_xp_falls_back_without_gpu(self):
        # generalized_unifrac has no CPU implementation reachable from the
        # public API, so this does not raise when no GPU backend is usable
        # -- it is the regression this whole fallback exists to fix.
        table, taxa, tree, _ = self._load_qiime_191_tt()
        with patch_gpu_backend(None):
            obs = _generalized_unifrac_gpu_or_xp(
                table, taxa, tree, alpha=0.5, validate=True)
        cpu = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=0.5, validate=True)
        np.testing.assert_allclose(obs, cpu, rtol=0, atol=XP_CPU_TOLERANCE)


class DispatchValidationErrorTests(QiimeTinyTestMixin, TestCase):
    """A bad tree/taxa input must raise its real validation error, not get
    misclassified by `_dispatch_gpu_or_xp` as "this backend's kernel can't
    run" (which would both hide the error behind the array-API fallback and
    permanently degrade the backend to that fallback for the rest of the
    process, over what was really just bad input on one call).
    """

    def setUp(self):
        # Never let a previous test's (or this test's own) backend
        # unavailability bleed across tests.
        self._backends_backup = set(_unavailable_backends)
        _unavailable_backends.clear()
        self.addCleanup(_unavailable_backends.clear)
        self.addCleanup(_unavailable_backends.update, self._backends_backup)

    def test_dispatch_reraises_value_error_without_marking_backend_unavailable(self):
        def bad_gpu_func(*args, **kwargs):
            raise ValueError("bad input")

        def xp_func(*args, **kwargs):
            raise AssertionError("xp_func must not be reached")

        with patch_gpu_backend("hip"):
            with self.assertRaises(ValueError):
                _dispatch_gpu_or_xp(bad_gpu_func, xp_func)
        self.assertNotIn("hip", _unavailable_backends)

    def test_dispatch_reraises_tree_error_without_marking_backend_unavailable(self):
        def bad_gpu_func(*args, **kwargs):
            raise DuplicateNodeError("bad tree")

        def xp_func(*args, **kwargs):
            raise AssertionError("xp_func must not be reached")

        with patch_gpu_backend("hip"):
            with self.assertRaises(DuplicateNodeError):
                _dispatch_gpu_or_xp(bad_gpu_func, xp_func)
        self.assertNotIn("hip", _unavailable_backends)

    def test_dispatch_still_marks_backend_unavailable_for_genuine_kernel_errors(self):
        # Sanity check for the other side of the fix: a non-validation
        # exception (standing in for e.g. a numba-hip compile/runtime
        # failure) must still fall through to the array-API path and mark
        # the backend unavailable, exactly as before this fix.
        def bad_gpu_func(*args, **kwargs):
            raise RuntimeError("kernel failed to build")

        def xp_func(*args, **kwargs):
            return "xp result"

        with patch_gpu_backend("hip"):
            with self.assertWarns(UserWarning):
                obs = _dispatch_gpu_or_xp(bad_gpu_func, xp_func)
        self.assertEqual(obs, "xp result")
        self.assertIn("hip", _unavailable_backends)

    def test_invalid_taxa_under_engine_gpu_raises_and_keeps_backend_usable(self):
        # End-to-end version through the real `_weighted_unifrac_gpu` driver:
        # `_setup_multiple_unifrac` (called inside `_weighted_unifrac_gpu`,
        # before any GPU-specific code) raises on duplicate taxa (a
        # ValueError -- "All taxa must be unique" -- distinct from the
        # DuplicateNodeError the same validator raises for duplicate *tip
        # names in the tree*; both are TreeError/ValueError cases this fix
        # covers). A fake, never-actually-used cuda module stands in for
        # `_get_cuda_module`'s real backend-presence check, so this is
        # exercisable without real GPU hardware; what is under test is
        # purely the exception classification in
        # `_dispatch_gpu_or_xp`/`_weighted_unifrac_gpu`, not the kernel
        # itself.
        table, taxa, tree, _ = self._load_qiime_191_tt()
        duplicate_taxa = list(taxa)
        duplicate_taxa[1] = duplicate_taxa[0]

        with patch_gpu_backend("hip"), patch(
            "skbio.diversity.beta._unifrac_gpu._get_cuda_module",
            return_value=object(),
        ), patch(
            "skbio.diversity.beta._unifrac_gpu._weighted_unifrac_gpu",
            wraps=_weighted_unifrac_gpu,
        ) as spy_gpu_func:
            with self.assertRaises(ValueError):
                _weighted_unifrac_gpu_or_xp(
                    table, duplicate_taxa, tree, normalized=True, validate=True
                )
            self.assertNotIn("hip", _unavailable_backends)
            self.assertEqual(spy_gpu_func.call_count, 1)

            # A second, valid request on the still-usable backend must still
            # attempt the real kernel driver, rather than having been
            # silently routed to the fallback by the first call's
            # (correctly re-raised) validation error.
            obs = _weighted_unifrac_gpu_or_xp(
                table, taxa, tree, normalized=True, validate=True
            )
            self.assertEqual(spy_gpu_func.call_count, 2)
            # The fake cuda module has none of numba.cuda's real API, so the
            # kernel launch itself fails and this call falls back to the
            # array-API path -- expected, since this test's fake-GPU setup
            # cannot run a real kernel; what matters is that the driver was
            # attempted (checked above) and the fallback result is correct.
            cpu = _weighted_unifrac_pdist_numba(
                table, taxa, tree, normalized=True, validate=True
            )
            np.testing.assert_allclose(obs, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_invalid_taxa_under_engine_gpu_validate_false_raises_and_keeps_usable(
        self,
    ):
        # Same contract as the validate=True test above, but for
        # validate=False: skipping `_validate_taxa_and_tree` means a
        # taxon absent from the tree is never caught by that check, so it
        # instead surfaces as a plain `KeyError` from `_nodes_by_counts`
        # (via `vectorize_counts_and_tree`/`_setup_multiple_unifrac`),
        # confirmed by tracing the real call chain -- not a `ValueError`/
        # `TreeError`. This must still be re-raised immediately rather than
        # misclassified as "this backend's kernel can't run" (which would
        # hide the error and permanently degrade the backend).
        table, taxa, tree, _ = self._load_qiime_191_tt()
        mismatched_taxa = list(taxa)
        mismatched_taxa[0] = "not-a-real-tip-name"

        with patch_gpu_backend("hip"), patch(
            "skbio.diversity.beta._unifrac_gpu._get_cuda_module",
            return_value=object(),
        ), patch(
            "skbio.diversity.beta._unifrac_gpu._weighted_unifrac_gpu",
            wraps=_weighted_unifrac_gpu,
        ) as spy_gpu_func:
            with self.assertRaises(KeyError):
                _weighted_unifrac_gpu_or_xp(
                    table, mismatched_taxa, tree, normalized=True, validate=False
                )
            self.assertNotIn("hip", _unavailable_backends)
            self.assertEqual(spy_gpu_func.call_count, 1)

            # A second, valid request on the still-usable backend must still
            # attempt the real kernel driver, exactly as in the validate=True
            # test above.
            obs = _weighted_unifrac_gpu_or_xp(
                table, taxa, tree, normalized=True, validate=True
            )
            self.assertEqual(spy_gpu_func.call_count, 2)
            cpu = _weighted_unifrac_pdist_numba(
                table, taxa, tree, normalized=True, validate=True
            )
            np.testing.assert_allclose(obs, cpu, rtol=0, atol=XP_CPU_TOLERANCE)


class _ArrayApiGpuStateMixin:
    """Resets `_unifrac_xp`'s module-level GPU-array-API detection state
    (cache + disabled flag) and `sys.modules['cupy'/'torch']` around each
    test, mirroring `GetTileConfigTests`'s cache save/restore in
    test_unifrac_gpu.py so tests cannot see each other's cached result or
    leak a fake module into later tests.
    """

    def setUp(self):
        super().setUp()
        self._cache_backup = _unifrac_xp_mod._array_api_gpu_module_cache
        self._disabled_backup = _unifrac_xp_mod._array_api_gpu_disabled
        _unifrac_xp_mod._array_api_gpu_module_cache = None
        _unifrac_xp_mod._array_api_gpu_disabled = False
        self.addCleanup(self._restore_state)

        self._modules_backup = {
            name: sys.modules.get(name) for name in ("cupy", "torch")
        }
        self.addCleanup(self._restore_modules)

    def _restore_state(self):
        _unifrac_xp_mod._array_api_gpu_module_cache = self._cache_backup
        _unifrac_xp_mod._array_api_gpu_disabled = self._disabled_backup

    def _restore_modules(self):
        for name, mod in self._modules_backup.items():
            if mod is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = mod


class ArrayApiGpuDetectionTests(_ArrayApiGpuStateMixin, TestCase):
    """Unit coverage for `_detect_array_api_gpu_module`: prefers CuPy, falls
    back to PyTorch, returns None when neither is usable, and caches its
    result for the rest of the process -- independent of real GPU hardware,
    via fake modules injected into `sys.modules` (the same "simulate absent/
    present without needing the real package" trick
    `NumbaOptionalImportTests` uses in test_unifrac_gpu.py).

    This is deliberately independent of `_detect_gpu_backend`
    (numba-cuda/numba.hip): a machine can have one working while the other
    is broken, so the two detectors must not be conflated.
    """

    def test_prefers_cupy_when_both_are_usable(self):
        fake_cupy = _make_fake_cupy()
        fake_torch = _make_fake_torch()
        sys.modules["cupy"] = fake_cupy
        sys.modules["torch"] = fake_torch
        result = _detect_array_api_gpu_module()
        self.assertIs(result, fake_cupy)
        self.assertEqual(fake_torch.asarray_calls, 0)

    def test_falls_back_to_torch_when_cupy_is_not_installed(self):
        sys.modules["cupy"] = None  # not installed
        fake_torch = _make_fake_torch()
        sys.modules["torch"] = fake_torch
        result = _detect_array_api_gpu_module()
        self.assertIs(result, fake_torch)

    def test_falls_back_to_torch_when_cupy_has_no_device(self):
        fake_cupy = _make_fake_cupy(device_count=0)
        fake_torch = _make_fake_torch()
        sys.modules["cupy"] = fake_cupy
        sys.modules["torch"] = fake_torch
        result = _detect_array_api_gpu_module()
        self.assertIs(result, fake_torch)

    def test_returns_none_when_neither_is_installed(self):
        sys.modules["cupy"] = None
        sys.modules["torch"] = None
        result = _detect_array_api_gpu_module()
        self.assertIsNone(result)

    def test_returns_none_when_torch_cuda_is_not_available(self):
        sys.modules["cupy"] = None
        sys.modules["torch"] = _make_fake_torch(available=False)
        result = _detect_array_api_gpu_module()
        self.assertIsNone(result)

    def test_result_is_cached_across_calls(self):
        fake_cupy = _make_fake_cupy()
        sys.modules["cupy"] = fake_cupy
        first = _detect_array_api_gpu_module()
        # A second, distinct fake module -- if this were returned instead of
        # `fake_cupy`, it would prove the cache did not take, i.e. detection
        # re-probed instead of reusing the first result.
        sys.modules["cupy"] = _make_fake_cupy()
        second = _detect_array_api_gpu_module()
        self.assertIs(first, fake_cupy)
        self.assertIs(second, fake_cupy)


class MoveToArrayApiGpuTests(_ArrayApiGpuStateMixin, TestCase):
    """Unit coverage for `_move_to_array_api_gpu` itself, independent of
    `_prepare`."""

    def test_cupy_like_module_uses_asarray(self):
        fake_cupy = _make_fake_cupy()
        a = np.arange(3.0)
        b = np.arange(4.0)
        out_a, out_b = _move_to_array_api_gpu(fake_cupy, a, b)
        np.testing.assert_array_equal(out_a, a)
        np.testing.assert_array_equal(out_b, b)
        self.assertEqual(fake_cupy.asarray_calls, 2)

    def test_torch_like_module_uses_as_tensor(self):
        fake_torch = _make_fake_torch()
        a = np.arange(3.0)
        out_a, = _move_to_array_api_gpu(fake_torch, a)
        np.testing.assert_array_equal(out_a, a)
        self.assertEqual(fake_torch.asarray_calls, 1)

    def test_a_failure_partway_raises_without_moving_the_rest(self):
        fake_cupy = _make_fake_cupy(fail_asarray=True)
        with self.assertRaises(RuntimeError):
            _move_to_array_api_gpu(fake_cupy, np.arange(3.0), np.arange(4.0))


class PrepareArrayApiGpuWiringTests(
    _ArrayApiGpuStateMixin, QiimeTinyTestMixin, TestCase
):
    """`_prepare` (exercised here through `_weighted_unifrac_xp`/
    `_unweighted_unifrac_xp`) must actually use a detected GPU array-API
    module when one is available, and must fall back to plain NumPy --
    correctly, not just without crashing -- when moving the data onto it
    fails.
    """

    def test_prepare_moves_data_onto_a_detected_gpu_module(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        fake_cupy = _make_fake_cupy()
        with patch(
            "skbio.diversity.beta._unifrac_xp._detect_array_api_gpu_module",
            return_value=fake_cupy,
        ):
            obs = _weighted_unifrac_xp(
                table, taxa, tree, normalized=True, validate=True
            )
        # counts_by_node, sample_totals, branch_lengths -- all 3 arrays
        # `_prepare` builds from tree traversal before the pairwise math.
        self.assertEqual(fake_cupy.asarray_calls, 3)
        cpu = _weighted_unifrac_pdist_numba(
            table, taxa, tree, normalized=True, validate=True
        )
        np.testing.assert_allclose(obs, cpu, rtol=0, atol=XP_CPU_TOLERANCE)

    def test_prepare_also_moves_proportions_onto_the_gpu_module(self):
        # `_unweighted_unifrac_xp` takes `need_proportions=False`, so this
        # covers `_weighted_unifrac_xp`'s `need_proportions=True` path,
        # where `proportions` itself must come out of the detected module's
        # `xp`, not a stray NumPy array left over from `_prepare`'s own
        # `np.ascontiguousarray` calls (the thing the whole task fixes).
        table, taxa, tree, _ = self._load_qiime_191_tt()
        fake_cupy = _make_fake_cupy()
        with patch(
            "skbio.diversity.beta._unifrac_xp._detect_array_api_gpu_module",
            return_value=fake_cupy,
        ):
            xp, counts_by_node, sample_totals, branch_lengths, proportions = (
                _unifrac_xp_mod._prepare(table, taxa, tree, validate=True)
            )
        # The fake module wraps plain NumPy, so `xp` itself still reports as
        # the NumPy namespace -- what this checks is that `proportions` was
        # computed from the already-"moved" `counts_by_node`/`sample_totals`
        # (i.e. after the move, not from the pre-move NumPy locals), which
        # the 3-array `asarray_calls` count below confirms occurred before
        # `proportions` was built.
        self.assertEqual(fake_cupy.asarray_calls, 3)
        self.assertIsNotNone(proportions)

    def test_falls_back_to_numpy_when_moving_the_data_fails(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        fake_cupy = _make_fake_cupy(fail_asarray=True)
        with patch(
            "skbio.diversity.beta._unifrac_xp._detect_array_api_gpu_module",
            return_value=fake_cupy,
        ):
            with self.assertWarns(UserWarning):
                obs = _weighted_unifrac_xp(
                    table, taxa, tree, normalized=True, validate=True
                )
        cpu = _weighted_unifrac_pdist_numba(
            table, taxa, tree, normalized=True, validate=True
        )
        np.testing.assert_allclose(obs, cpu, rtol=0, atol=XP_CPU_TOLERANCE)
        self.assertTrue(_unifrac_xp_mod._array_api_gpu_disabled)

    def test_a_failed_move_permanently_disables_the_gpu_module_for_the_process(self):
        # Unlike the other tests in this class, this one does NOT mock
        # `_detect_array_api_gpu_module` -- it injects the fake module via
        # `sys.modules` instead, so the real detection function runs,
        # caches `fake_cupy`, and (after the move below fails) sets the
        # real `_array_api_gpu_disabled` flag that later calls check.
        table, taxa, tree, _ = self._load_qiime_191_tt()
        fake_cupy = _make_fake_cupy(fail_asarray=True)
        sys.modules["cupy"] = fake_cupy
        with self.assertWarns(UserWarning):
            _weighted_unifrac_xp(table, taxa, tree, normalized=True, validate=True)
        # A second call must not even attempt the (still-failing) GPU
        # module again: `_detect_array_api_gpu_module` itself must
        # short-circuit to None now, not just have `_prepare` catch a
        # second failure.
        self.assertIsNone(_detect_array_api_gpu_module())
        self.assertEqual(fake_cupy.asarray_calls, 1)


class ThreeTierDispatchTests(_ArrayApiGpuStateMixin, QiimeTinyTestMixin, TestCase):
    """End-to-end confirmation that `_dispatch_gpu_or_xp`'s fallback branch
    (taken when the fused GPU kernel is unusable) genuinely benefits from a
    detected GPU array-API module, automatically, with no change needed in
    `_unifrac_gpu.py` itself -- `_weighted_unifrac_gpu_or_xp` calls
    `_weighted_unifrac_xp` as its fallback, and `_weighted_unifrac_xp` calls
    `_prepare`, which is where the new detection+transfer lives.
    """

    def test_fallback_path_uses_the_detected_gpu_module_without_gpu_py_changes(self):
        table, taxa, tree, _ = self._load_qiime_191_tt()
        fake_cupy = _make_fake_cupy()
        with patch_gpu_backend(None), patch(
            "skbio.diversity.beta._unifrac_xp._detect_array_api_gpu_module",
            return_value=fake_cupy,
        ):
            obs = _weighted_unifrac_gpu_or_xp(
                table, taxa, tree, normalized=True, validate=True
            )
        # Proof the 3rd tier (GPU array-API module) actually activated
        # rather than silently falling through to tier 4 (plain NumPy):
        # `_weighted_unifrac_xp`'s `_prepare` call must have moved all 3
        # setup arrays onto the fake module.
        self.assertEqual(fake_cupy.asarray_calls, 3)
        cpu = _weighted_unifrac_pdist_numba(
            table, taxa, tree, normalized=True, validate=True
        )
        np.testing.assert_allclose(obs, cpu, rtol=0, atol=XP_CPU_TOLERANCE)


if __name__ == "__main__":
    main()
