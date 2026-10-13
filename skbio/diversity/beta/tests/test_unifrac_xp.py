# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

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
from skbio.diversity.beta._unifrac_xp import (
    _generalized_unifrac_xp,
    _unweighted_unifrac_xp,
    _weighted_unifrac_xp,
)
from skbio.diversity.beta.tests._fixtures import (
    QiimeTinyTestMixin,
    patch_gpu_backend,
)

# Measured max abs deviation, array-API path vs CPU-numba kernel, for the
# method/normalized combinations both implement (no variance_adjust) on the
# qiime-191-tt table: ~3.3e-16 (ordinary float64 noise). 1e-10 matches the
# tolerance used for the equivalent GPU-vs-CPU-numba comparison in
# test_unifrac_gpu.py.
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


if __name__ == "__main__":
    main()
