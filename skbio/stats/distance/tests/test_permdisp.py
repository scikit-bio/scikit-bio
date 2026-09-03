# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial
from unittest import TestCase, main, skipIf
import platform

import numpy as np
import numpy.testing as npt
import pandas as pd
from pandas.testing import assert_series_equal
from scipy.stats import f_oneway, ConstantInputWarning

from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa, OrdinationResults
from skbio.stats.distance import permdisp
from skbio.stats.distance._permdisp import _compute_groups
from skbio.stats.distance._cutils import geomedian_axis_one
from skbio.util import get_data_path, numba_code

IS_INTEL_MAC = platform.system() == "Darwin" and platform.machine() == "x86_64"

class PERMDISPTests(TestCase):

    def setUp(self):
        # test with 2 groups of equal size
        # when assigned different labels, results should be the same
        self.grouping_eq = ['foo', 'foo', 'foo', 'bar', 'bar', 'bar']
        self.grouping_eq_relab = ['pyt', 'pyt', 'pyt', 'hon', 'hon', 'hon']
        self.exp_index = ['method name', 'test statistic name', 'sample size',
                          'number of groups', 'test statistic', 'p-value',
                          'number of permutations']
        # test with 3 groups of different sizes
        # when assigned different labels results should be the same
        self.grouping_uneq = ['foo', 'foo', 'bar', 'bar', 'bar',
                              'qw', 'qw', 'qw', 'qw']

        self.grouping_uneq_relab = [12, 12, 7, 7, 7, 23, 23, 23, 23]

        self.grouping_un_mixed = ['a', 'a', 7, 7, 7, 'b', 'b', 'b', 'b']

        eq_ids = ['s1', 's2', 's3', 's4', 's5', 's6']
        uneq_ids = ['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9']
        # matrix for equal grouping
        self.eq_mat = DistanceMatrix([[0, 4, 0, 0, 4, 2],
                                      [4, 0, 2, 0, 3, 1],
                                      [0, 2, 0, 5, 2, 5],
                                      [0, 0, 5, 0, 0, 2],
                                      [4, 3, 2, 0, 0, 2],
                                      [2, 1, 5, 2, 2, 0]], eq_ids)

        # matrix for unequal grouping
        self.uneq_mat = DistanceMatrix([[0, 0, 4, 0, 0, 3, 5, 3, 0],
                                        [0, 0, 0, 3, 4, 5, 3, 0, 3],
                                        [4, 0, 0, 4, 3, 1, 0, 5, 2],
                                        [0, 3, 4, 0, 0, 2, 1, 3, 5],
                                        [0, 4, 3, 0, 0, 1, 1, 5, 0],
                                        [3, 5, 1, 2, 1, 0, 2, 0, 5],
                                        [5, 3, 0, 1, 1, 2, 0, 4, 3],
                                        [3, 0, 5, 3, 5, 0, 4, 0, 4],
                                        [0, 3, 2, 5, 0, 5, 3, 4, 0]], uneq_ids)

        # null matrix for equal grouping
        self.null_mat = DistanceMatrix([[0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0]], eq_ids)

        unif_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
                    'PC.634', 'PC.635', 'PC.636']

        self.unifrac_dm = DistanceMatrix(
            [[0.0, 0.595483768391, 0.618074717633, 0.582763100909,
              0.566949022108, 0.714717232268, 0.772001731764, 0.690237118413,
              0.740681707488],
             [0.595483768391, 0.0, 0.581427669668, 0.613726772383,
              0.65945132763, 0.745176523638, 0.733836123821, 0.720305073505,
              0.680785600439],
             [0.618074717633, 0.581427669668, 0.0, 0.672149021573,
              0.699416863323, 0.71405573754, 0.759178215168, 0.689701276341,
              0.725100672826],
             [0.582763100909, 0.613726772383, 0.672149021573, 0.0,
              0.64756120797, 0.666018240373, 0.66532968784, 0.650464714994,
              0.632524644216],
             [0.566949022108, 0.65945132763, 0.699416863323, 0.64756120797,
              0.0, 0.703720200713, 0.748240937349, 0.73416971958,
              0.727154987937],
             [0.714717232268, 0.745176523638, 0.71405573754, 0.666018240373,
              0.703720200713, 0.0, 0.707316869557, 0.636288883818,
              0.699880573956],
             [0.772001731764, 0.733836123821, 0.759178215168, 0.66532968784,
              0.748240937349, 0.707316869557, 0.0, 0.565875193399,
              0.560605525642],
             [0.690237118413, 0.720305073505, 0.689701276341, 0.650464714994,
              0.73416971958, 0.636288883818, 0.565875193399, 0.0,
              0.575788039321],
             [0.740681707488, 0.680785600439, 0.725100672826, 0.632524644216,
              0.727154987937, 0.699880573956, 0.560605525642, 0.575788039321,
              0.0]], unif_ids)
        
        # condensed form
        self.unifrac_dm_condensed = DistanceMatrix(self.unifrac_dm)

        self.unif_grouping = ['Control', 'Control', 'Control', 'Control',
                              'Control', 'Fast', 'Fast', 'Fast', 'Fast']

        self.assert_series_equal = partial(assert_series_equal,
                                           check_index_type=True,
                                           check_series_type=True)

    def test_centroids_eq_groups(self):
        exp = [[1.2886811963240687, 1.890538910062923, 1.490527658097728],
               [2.17349240061718, 2.3192679626679946, 2.028338553903792]]
        exp_stat, _ = f_oneway(*exp)

        dm = pcoa(self.eq_mat, warn_neg_eigval=False)
        dm = dm.samples

        obs = _compute_groups(dm, 'centroid', self.grouping_eq)
        self.assertAlmostEqual(obs, exp_stat, places=6)

        obs_relab = _compute_groups(dm, 'centroid', self.grouping_eq_relab)
        self.assertAlmostEqual(obs_relab, obs, places=6)

    def test_centroids_uneq_groups(self):
        """
        the expected result here was calculated by hand
        """
        exp = [[2.5847022428144935, 2.285624595858895,
                1.7022431146340287],
               [1.724817266046108, 1.724817266046108],
               [2.4333280644972795, 2.389000390879655,
                2.8547180589306036, 3.218568759338847]]
        exp_stat, _ = f_oneway(*exp)

        dm = pcoa(self.uneq_mat, warn_neg_eigval=False)
        dm = dm.samples

        obs = _compute_groups(dm, 'centroid', self.grouping_uneq)
        self.assertAlmostEqual(obs, exp_stat, places=6)

        obs_relab = _compute_groups(dm, 'centroid', self.grouping_uneq_relab)
        self.assertAlmostEqual(obs, obs_relab, places=6)

    def test_centroids_mixedgroups(self):
        exp = [[2.5847022428144935, 2.285624595858895,
                1.7022431146340287],
               [1.724817266046108, 1.724817266046108],
               [2.4333280644972795, 2.389000390879655,
                2.8547180589306036, 3.218568759338847]]
        dm = pcoa(self.uneq_mat, warn_neg_eigval=False)
        dm = dm.samples

        exp_stat, _ = f_oneway(*exp)

        obs_mixed = _compute_groups(dm, 'centroid', self.grouping_un_mixed)
        self.assertAlmostEqual(exp_stat, obs_mixed, places=6)

    def test_centroids_null(self):
        dm = pcoa(self.null_mat)
        dm = dm.samples

        # TODO: SciPy 1.17+ changed the behavior of f_oneway so that a
        # ConstantInputWarning may no longer be emitted for this case.
        # Once the expected behavior is clarified, consider reintroducing
        # an explicit warning assertion here.
        obs_null = _compute_groups(dm, 'centroid', self.grouping_eq)
        np.isnan(obs_null)

    def test_centroid_normal(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMDISP', 'F-value', 9, 2, 0.244501519876,
                              0.53, 99],
                        name='PERMDISP results')

        grouping = ['Control', 'Control', 'Control', 'Control', 'Control',
                    'Fast', 'Fast', 'Fast', 'Fast']

        obs = permdisp(self.unifrac_dm, grouping, test='centroid',
                       permutations=99, seed=42,
                       dimensions=self.unifrac_dm.shape[0])

        self.assert_series_equal(obs, exp)

    def test_centroid_normal_condensed(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMDISP', 'F-value', 9, 2, 0.244501519876,
                              0.53, 99],
                        name='PERMDISP results')

        grouping = ['Control', 'Control', 'Control', 'Control', 'Control',
                    'Fast', 'Fast', 'Fast', 'Fast']

        obs = permdisp(self.unifrac_dm_condensed, grouping, test='centroid',
                       permutations=99, seed=42,
                       dimensions=self.unifrac_dm_condensed.shape[0])

        self.assert_series_equal(obs, exp)

    @skipIf(IS_INTEL_MAC, "See issue #2382.")
    def test_median_normal(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMDISP', 'F-value', 9, 2, 0.139475441876,
                              0.51, 99],
                        name='PERMDISP results')

        obs = permdisp(self.unifrac_dm, self.unif_grouping, test='median',
                       permutations=99, seed=42,
                       dimensions=self.unifrac_dm.shape[0])

        self.assert_series_equal(obs, exp)

        po = pcoa(self.unifrac_dm)

        obs2 = permdisp(po, self.unif_grouping, test='median',
                        permutations=99, seed=42)

        self.assert_series_equal(obs2, exp)

    @skipIf(IS_INTEL_MAC, "See issue #2382.")
    def test_median_normal_condensed(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMDISP', 'F-value', 9, 2, 0.139475441876,
                              0.51, 99],
                        name='PERMDISP results')

        obs = permdisp(self.unifrac_dm_condensed, self.unif_grouping, test='median',
                       permutations=99, seed=42,
                       dimensions=self.unifrac_dm_condensed.shape[0])

        self.assert_series_equal(obs, exp)

        po = pcoa(self.unifrac_dm_condensed)

        obs2 = permdisp(po, self.unif_grouping, test='median',
                        permutations=99, seed=42)

        self.assert_series_equal(obs2, exp)

    @skipIf(IS_INTEL_MAC, "See issue #2382.")
    def test_median_fsvd(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMDISP', 'F-value', 9, 2, 0.04078077215673714,
                              0.79, 99],
                        name='PERMDISP results')

        obs = permdisp(self.unifrac_dm, self.unif_grouping, test='median',
                       permutations=99,
                       method='fsvd', dimensions=3, seed=42)

        self.assert_series_equal(obs, exp)

        po = pcoa(self.unifrac_dm, method='fsvd', dimensions=3)
        obs = permdisp(po, self.unif_grouping, test='median',
                       permutations=99, seed=42)

        self.assert_series_equal(obs, exp)

    @skipIf(IS_INTEL_MAC, "See issue #2382.")
    def test_median_fsvd_condensed(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMDISP', 'F-value', 9, 2, 0.04078077215673714,
                              0.79, 99],
                        name='PERMDISP results')

        obs = permdisp(self.unifrac_dm_condensed, self.unif_grouping, test='median',
                       permutations=99,
                       method='fsvd', dimensions=3, seed=42)

        self.assert_series_equal(obs, exp)

        po = pcoa(self.unifrac_dm_condensed, method='fsvd', dimensions=3)
        obs = permdisp(po, self.unif_grouping, test='median',
                       permutations=99, seed=42)

        self.assert_series_equal(obs, exp)

    def test_fsvd_seed_is_reproducible(self):
        # fsvd draws a random projection, so permdisp has to pass its seed down
        # to pcoa for the result to be reproducible. The matrix has to be large
        # enough that different projections actually disagree; on a handful of
        # samples they all converge and the difference is invisible. No
        # permutations are needed, since the statistic alone depends on the
        # ordination.
        rng = np.random.default_rng(0)
        mat = rng.random((60, 60))
        mat = (mat + mat.T) / 2
        np.fill_diagonal(mat, 0.0)
        dm = DistanceMatrix(mat)
        grouping = ['a'] * 30 + ['b'] * 30

        for test in ('centroid', 'median'):
            first = permdisp(dm, grouping, test=test, permutations=0,
                             method='fsvd', seed=42)
            repeat = permdisp(dm, grouping, test=test, permutations=0,
                              method='fsvd', seed=42)
            npt.assert_allclose(repeat['test statistic'],
                                first['test statistic'])

    def test_not_distance_matrix(self):
        dm = []
        grouping = ['Control', 'Control', 'Control', 'Control', 'Control',
                    'Fast', 'Fast', 'Fast', 'Fast']

        npt.assert_raises(TypeError, permdisp, dm, grouping, permutations=0)

    def test_mismatched_group(self):

        gr = ['foo', 'bar']
        npt.assert_raises(ValueError, permdisp, self.unifrac_dm, gr)

    def test_single_group(self):

        gr = ['f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']
        npt.assert_raises(ValueError, permdisp, self.unifrac_dm, gr)

    def test_no_permuations(self):
        obs = permdisp(self.eq_mat, self.grouping_eq, permutations=0,
                       dimensions=self.eq_mat.shape[0],
                       warn_neg_eigval=False)

        pval = obs['p-value']
        np.isnan(pval)

    def test_geomedian(self):
        exp = np.array([2.01956244, 1.53164546, 2.60571752, 0.91424179,
                        1.76214416, 1.69943057])
        obs = np.array(geomedian_axis_one(self.eq_mat.data))
        npt.assert_almost_equal(obs, exp, decimal=6)

    def test_confirm_betadispr_results(self):
        mp_dm = DistanceMatrix.read(get_data_path('moving_pictures_dm.tsv'))
        mp_mf = pd.read_csv(get_data_path('moving_pictures_mf.tsv'), sep='\t')
        mp_mf.set_index('#SampleID', inplace=True)

        obs_med_mp = permdisp(mp_dm, mp_mf, column='BodySite', seed=42,
                              dimensions=mp_dm.shape[0])
        obs_cen_mp = permdisp(mp_dm, mp_mf, column='BodySite', test='centroid',
                              seed=42, dimensions=mp_dm.shape[0])

        exp_data_m = ['PERMDISP', 'F-value', 33, 4, 10.1956, 0.001, 999]
        exp_data_c = ['PERMDISP', 'F-value', 33, 4, 17.4242, 0.001, 999]
        exp_ind = ['method name', 'test statistic name', 'sample size',
                   'number of groups', 'test statistic', 'p-value',
                   'number of permutations']

        exp_med_mp = pd.Series(data=exp_data_m, index=exp_ind, dtype='object',
                               name='PERMDISP results')

        exp_cen_mp = pd.Series(data=exp_data_c, index=exp_ind, dtype='object',
                               name='PERMDISP results')

        self.assert_series_equal(exp_med_mp, obs_med_mp)

        self.assert_series_equal(exp_cen_mp, obs_cen_mp)

    def test_call_via_series(self):
        # test https://github.com/scikit-bio/scikit-bio/issues/1877
        # actual issue is with _base._preprocess_input_sng but permdisp is
        # indirectly affected
        dm = DistanceMatrix.read(get_data_path('frameSeries_dm.tsv'))
        grouping = pd.read_csv(get_data_path("frameSeries_grouping.tsv"),
                               sep="\t", index_col=0)

        obs_frame = permdisp(dm, grouping, column='tumor', seed=42,
                             warn_neg_eigval=False)

        obs_series = permdisp(dm, grouping['tumor'], seed=42, warn_neg_eigval=False)

        # in principle, both tests - if seed is the same - should return the
        # exact same results. However, they don't for the current example ...
        self.assert_series_equal(obs_frame, obs_series)

    def test_invalid_test(self):
        with self.assertRaises(ValueError) as cm:
            permdisp(self.unifrac_dm, [], test='invalid')
        self.assertEqual(str(cm.exception), "Test must be centroid or median.")

    def test_invalid_method(self):
        with self.assertRaises(ValueError) as cm:
            permdisp(self.unifrac_dm, [], method='invalid')
        self.assertEqual(str(cm.exception), "Method must be eigh or fsvd.")

    def test_no_mutation_of_ordination_results(self):
        po = pcoa(self.unifrac_dm, warn_neg_eigval=False)
        original_columns = po.samples.columns.tolist()
        
        permdisp(po, self.unif_grouping, permutations=0, warn_neg_eigval=False)
        
        self.assertEqual(po.samples.columns.tolist(), original_columns)
        self.assertNotIn("grouping", po.samples.columns.tolist())

    def test_permdisp_passes_dimensions_through(self):
        # Regression test for issue #2430. Previously `permdisp` silently
        # overrode `dimensions = 0` before calling pcoa, which tripped
        # pcoa's "EIGH: since no value for dimensions..." warning on every
        # call for distance matrices larger than 10 samples. `permdisp`
        # now passes `dimensions` through unchanged: normal usage (an
        # explicit positive value, or the default of 10) no longer leaks
        # the warning, while users who explicitly opt in with
        # `dimensions=0` still receive it (per #2456 review).
        import warnings
        from scipy.spatial.distance import pdist, squareform

        rng = np.random.default_rng(0)
        n = 12  # exceed pcoa's 10-sample warning threshold
        coords = rng.random((n, 5))
        big_dm = DistanceMatrix(
            squareform(pdist(coords)),
            ids=[f"s{i}" for i in range(n)],
        )
        grouping = ["A"] * 6 + ["B"] * 6

        # Case 1: default permdisp call (dimensions=10) — no warning.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            permdisp(big_dm, grouping, permutations=9, warn_neg_eigval=False)
        self.assertFalse(
            any("EIGH" in str(w.message) for w in caught),
            "permdisp with a positive dimensions value must not emit "
            "pcoa's 'EIGH: since no value for dimensions...' warning",
        )

        # Case 2: user explicitly passes dimensions=0 — warning preserved.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            permdisp(big_dm, grouping, permutations=9,
                     dimensions=0, warn_neg_eigval=False)
        self.assertTrue(
            any(issubclass(w.category, RuntimeWarning)
                and "EIGH" in str(w.message) for w in caught),
            "permdisp must surface pcoa's EIGH warning when the caller "
            "explicitly passes dimensions=0 on a large matrix",
        )


class PERMDISPEngineTests(TestCase):
    """The numba engine must agree with cython, statistic and p-value alike."""

    def setUp(self):
        # A few well separated groups over a random matrix, large enough that
        # the pcoa retains the default ten dimensions.
        rng = np.random.default_rng(0)
        mat = rng.random((30, 30))
        mat = (mat + mat.T) / 2
        np.fill_diagonal(mat, 0.0)
        self.dm = DistanceMatrix(mat)
        self.grouping = ['a'] * 10 + ['b'] * 10 + ['c'] * 10

    def _assert_engines_agree(self, test, **kwargs):
        obs = permdisp(self.dm, self.grouping, test=test, seed=42,
                       engine="numba", **kwargs)
        exp = permdisp(self.dm, self.grouping, test=test, seed=42,
                       engine="cython", **kwargs)
        npt.assert_allclose(obs['test statistic'], exp['test statistic'])
        # The p-value is the real check: it only matches if the numba driver
        # draws permutations in the same order as the cython path.
        self.assertEqual(obs['p-value'], exp['p-value'])

    @numba_code
    def test_centroid_matches_cython(self):
        self._assert_engines_agree("centroid", permutations=99)

    @numba_code
    def test_median_matches_cython(self):
        self._assert_engines_agree("median", permutations=99)

    @numba_code
    def test_p_value_matches_cython_across_seeds(self):
        # The p-value only matches if the numba driver consumes the random
        # stream exactly as the cython path does. A desynchronized stream
        # still agrees on most seeds, so one seed is not enough to catch it.
        for seed in range(10):
            obs = permdisp(self.dm, self.grouping, test="median",
                           permutations=99, seed=seed, engine="numba")
            exp = permdisp(self.dm, self.grouping, test="median",
                           permutations=99, seed=seed, engine="cython")
            self.assertEqual(obs['p-value'], exp['p-value'])

    @numba_code
    def test_uneven_groups_match_cython(self):
        self.grouping = ['a'] * 4 + ['b'] * 6 + ['c'] * 20
        self._assert_engines_agree("median", permutations=99)

    @numba_code
    def test_dimensions_passed_through(self):
        self._assert_engines_agree("median", permutations=99, dimensions=3)

    @numba_code
    def test_ordination_input_matches_cython(self):
        # permdisp also accepts an OrdinationResults directly, which skips the
        # internal pcoa call. The engine dispatch happens after that branch, so
        # both input forms have to agree.
        ordination = pcoa(self.dm, method="eigh", dimensions=10)
        for test in ("centroid", "median"):
            obs = permdisp(ordination, self.grouping, test=test,
                           permutations=99, seed=42, engine="numba")
            exp = permdisp(ordination, self.grouping, test=test,
                           permutations=99, seed=42, engine="cython")
            npt.assert_allclose(obs['test statistic'], exp['test statistic'])
            self.assertEqual(obs['p-value'], exp['p-value'])

    @numba_code
    def test_no_permutations(self):
        # permutations=0 skips the batch entirely and returns a NaN p-value.
        obs = permdisp(self.dm, self.grouping, permutations=0, engine="numba")
        exp = permdisp(self.dm, self.grouping, permutations=0, engine="cython")
        npt.assert_allclose(obs['test statistic'], exp['test statistic'])
        self.assertTrue(np.isnan(obs['p-value']))

    @numba_code
    def test_geomedian_matches_cython_kernel(self):
        # The numba geomedian is a port of geomedian_axis_one and has to
        # return the same center, including for the degenerate single-column
        # case that returns the mean without iterating.
        from skbio.stats.distance._permdisp import _geomedian_nb

        rng = np.random.default_rng(1)
        for shape in [(10, 30), (3, 50), (10, 1), (10, 2)]:
            X = np.ascontiguousarray(rng.random(shape))
            npt.assert_allclose(_geomedian_nb(X.copy()),
                                np.asarray(geomedian_axis_one(X.copy())))

        # Layouts where the iteration lands exactly on one of the samples.
        # That is the one path needing the zero-distance correction, and these
        # particular ones move to a different center without it, so they pin
        # the correction down rather than merely reaching it.
        degenerate = [np.array([[2.0, 1.0, -2.0, 2.0, -3.0, 0.0]]),
                      np.array([[0.0, 3.0, -2.0, -1.0, 2.0, -2.0]]),
                      np.array([[-2.0, -1.0, -1.0, 0.0, -1.0, 1.0]]),
                      np.array([[-1.0, 1.0, 0.0, 0.0, 0.0],
                                [0.0, 0.0, -1.0, 1.0, 0.0]])]
        for X in degenerate:
            X = np.ascontiguousarray(X)
            npt.assert_allclose(_geomedian_nb(X.copy()),
                                np.asarray(geomedian_axis_one(X.copy())))

    @numba_code
    def test_geomedian_all_points_coincide(self):
        # Every sample sitting on the same spot leaves the iteration with no
        # direction to move in, and the median is that spot. Not compared
        # against geomedian_axis_one here because the Cython kernel divides by
        # a zero weight total before reaching that test and raises instead.
        from skbio.stats.distance._permdisp import _geomedian_nb

        X = np.ascontiguousarray(np.tile(np.array([[1.5], [-2.0]]), (1, 6)))
        npt.assert_allclose(_geomedian_nb(X), np.array([1.5, -2.0]))

    @numba_code
    def test_degenerate_anova_matches_f_oneway(self):
        # The numba kernels inline the one-way ANOVA instead of calling
        # scipy.stats.f_oneway, so the two have to agree where the ratio stops
        # being finite. Putting each group on a shell of constant radius about
        # its own center makes every within-group distance identical, which
        # zeroes the within-group term: differing radii then leave the groups
        # apart and give inf, equal radii give nan.
        from skbio.stats.distance._permdisp import _permdisp_f_stat_centroid_nb

        def shell(cx, cy, r, k):
            th = np.linspace(0, 2 * np.pi, k, endpoint=False)
            return np.stack([cx + r * np.cos(th), cy + r * np.sin(th)], axis=1)

        codes = np.array([0] * 4 + [1] * 4)
        for radius in (2.0, 1.0):
            samples = np.vstack([shell(0.0, 0.0, 1.0, 4),
                                 shell(10.0, 10.0, radius, 4)])
            exp = _compute_groups(samples, "centroid", codes)
            obs = _permdisp_f_stat_centroid_nb(
                np.ascontiguousarray(samples),
                np.ascontiguousarray(codes, dtype=np.int32), 2)
            npt.assert_array_equal(obs, exp)

    @numba_code
    def test_tied_statistics_match_cython(self):
        # With few samples per group many permutations reproduce the observed
        # partition under a different labeling, which gives a statistic that
        # is exactly equal to the observed one rather than merely close. Those
        # ties count towards the p-value, so this pins the comparison down to
        # the same ">=" the cython path uses. Only the centroid test is checked
        # this way. _compute_groups builds its group list in np.unique label
        # order and hands it to scipy.stats.f_oneway, which is not invariant to
        # the order of its arguments; on the median distances that costs a last
        # bit, so a relabeled partition does not always reproduce the statistic
        # and the tie counts legitimately differ between the engines.
        rng = np.random.default_rng(0)
        mat = rng.random((6, 6))
        mat = (mat + mat.T) / 2
        np.fill_diagonal(mat, 0.0)
        dm = DistanceMatrix(mat)
        grouping = ['a'] * 3 + ['b'] * 3
        obs = permdisp(dm, grouping, test="centroid", permutations=999, seed=42,
                       engine="numba", dimensions=5, warn_neg_eigval=False)
        exp = permdisp(dm, grouping, test="centroid", permutations=999, seed=42,
                       engine="cython", dimensions=5, warn_neg_eigval=False)
        npt.assert_allclose(obs['test statistic'], exp['test statistic'])
        self.assertEqual(obs['p-value'], exp['p-value'])

    @numba_code
    def test_negative_permutations_raises(self):
        with self.assertRaisesRegex(ValueError, "greater than or equal to zero"):
            permdisp(self.dm, self.grouping, permutations=-1, engine="numba")

    @numba_code
    def test_crosses_chunk_boundary(self):
        # The driver batches the permutations, so more than one chunk's worth
        # exercises the multi-batch path. The RNG has to stay in step across
        # chunks for the p-value to still match the cython engine.
        from numba import get_num_threads
        from skbio.stats.distance._permdisp import _PERM_CHUNK_PER_THREAD

        permutations = _PERM_CHUNK_PER_THREAD * get_num_threads() + 5
        obs = permdisp(self.dm, self.grouping, permutations=permutations,
                       seed=42, engine="numba")
        exp = permdisp(self.dm, self.grouping, permutations=permutations,
                       seed=42, engine="cython")
        npt.assert_allclose(obs['test statistic'], exp['test statistic'])
        self.assertEqual(obs['p-value'], exp['p-value'])

    @numba_code
    def test_non_float_ordination_raises_like_cython(self):
        # geomedian_axis_one only has single and double precision signatures,
        # so the cython path raises on anything else. The numba kernels would
        # happily widen an integer array and return a number instead, which is
        # worse than the error, so the engine refuses the same inputs.
        ordination = pcoa(self.dm, method="eigh", dimensions=10)

        def recast(dtype):
            return OrdinationResults(
                short_method_name=ordination.short_method_name,
                long_method_name=ordination.long_method_name,
                eigvals=ordination.eigvals,
                samples=ordination.samples.astype(dtype),
                proportion_explained=ordination.proportion_explained)

        for dtype in ("float16", "int64"):
            for engine in ("cython", "numba"):
                with self.assertRaises(TypeError):
                    permdisp(recast(dtype), self.grouping, test="median",
                             permutations=0, engine=engine)
        # single and double precision stay accepted by both
        for dtype in ("float32", "float64"):
            for engine in ("cython", "numba"):
                permdisp(recast(dtype), self.grouping, test="median",
                         permutations=0, engine=engine)

    @numba_code
    def test_result_types_match_cython(self):
        # permdisp returns an object-dtype Series, so the type of each entry is
        # visible to the caller and must not depend on the engine.
        for permutations in (0, 99):
            obs = permdisp(self.dm, self.grouping, permutations=permutations,
                           seed=42, engine="numba")
            exp = permdisp(self.dm, self.grouping, permutations=permutations,
                           seed=42, engine="cython")
            for key in ("test statistic", "p-value"):
                self.assertIs(type(obs[key]), type(exp[key]))

    def test_bad_engine_raises(self):
        with self.assertRaisesRegex(ValueError, "engine='julia' is not supported"):
            permdisp(self.dm, self.grouping, permutations=0, engine="julia")


if __name__ == '__main__':
    main()
