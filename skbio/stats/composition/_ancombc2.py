# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# This implementation of ANCOM-BC is based on an analysis of the source code
# from the R package ANCOMBC:
# - https://github.com/FrederickHuangLin/ANCOMBC
#
# Which is licensed under Artistic-2.0:
# - https://www.bioconductor.org/packages/release/bioc/html/ANCOMBC.html
#
# We thank Dr. Huang Lin (@FrederickHuangLin) for his helpful advice.
# ----------------------------------------------------------------------------

from __future__ import annotations

from typing import Optional
from typing import TYPE_CHECKING
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import norm, chi2, f, t
from scipy.optimize import minimize
from patsy import dmatrix

from skbio.util import get_rng
from skbio.table._tabular import _ingest_table, _aggregate_features
from ._base import _check_composition
from ._utils import _check_metadata, _check_p_adjust, _type_cast_to_float

if TYPE_CHECKING:  # pragma: no cover
    from skbio.util._typing import SeedLike


def ancombc(
    table,
    metadata,
    formula,
    grouping=None,
    max_iter=100,
    tol=1e-5,
    alpha=0.05,
    p_adjust="holm",
):
    r"""Perform differential abundance test using ANCOM-BC.

    Analysis of compositions of microbiomes with bias correction (ANCOM-BC) [1]_ is a
    differential abundance testing method featuring the estimation and correction for
    the bias of differential sampling fractions.

    .. versionadded:: 0.7.1

    .. versionchanged:: 0.7.4
        Fixed a bug in the global test, which would produce inaccurate results. Please
        update the program. The main results are not impacted.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing strictly positive count or proportional abundance data of
        the samples. See :ref:`supported formats <table_like>`.
    metadata : pd.DataFrame or 2-D array_like
        Metadata of the samples. Rows correspond to samples and columns correspond
        to covariates (attributes). Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    formula : str or generic Formula object
        Formula defining the model using factors included in the metadata columns.
        Refer to `Patsy's documentation
        <https://patsy.readthedocs.io/en/latest/formulas.html>`_ on how to specify
        a formula.
    grouping : str, optional
        Metadata column defining sample groups for post-hoc analyses. Must be a term
        in ``formula`` and contain at least three groups. Other terms are treated as
        adjustment covariates. It has no impact on the main esult. If None (default),
        post-hoc analyses will be unavailable.
    max_iter : int, optional
        Maximum number of iterations for the bias estimation process. Default is 100.
    tol : float, optional
        Absolute convergence tolerance for the bias estimation process. Default is
        1e-5.
    alpha : float, optional
        Significance level for the statistical tests. Must be in the range of (0, 1).
        Default is 0.05.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are
        Holm-Bonferroni ("holm" or "holm-bonferroni") (default), Benjamini-Hochberg
        ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported by
        statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.

    Returns
    -------
    :class:`ANCOMBCResult`
        Result object with primary results and post-hoc analysis methods.

    See Also
    --------
    ancombc2 : ANCOM-BC2 with multi-group tests and sensitivity analysis.
    struc_zero : Standalone structural zero detection.

    References
    ----------
    .. [1] Lin, H. and Peddada, S.D., 2020. Analysis of compositions of microbiomes
       with bias correction. Nature Communications, 11(1), p.3514.

    Examples
    --------
    >>> from skbio.stats.composition import ancombc
    >>> import pandas as pd

    **A basic example**

    Let's create a data table with six samples and seven features (e.g., these may be
    microbial taxa):

    >>> samples = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    >>> features = ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7']
    >>> table = pd.DataFrame(
    ...     [[ 2,  0,  4,  7,  0,  0,  1],
    ...      [ 1,  1,  0,  6,  5,  1, 10],
    ...      [ 3,  0,  2,  9,  6,  0,  1],
    ...      [ 0, 12,  1,  2,  0,  3,  2],
    ...      [ 2,  2, 37,  0,  0,  7,  3],
    ...      [10,  1,  0,  0,  4,  4,  3]],
    ...     index=samples, columns=features)

    Then create a sampling grouping vector. In this example, there is a "healthy" group
    and a "sick" group.

    >>> grouping = ['healthy', 'healthy', 'healthy', 'sick', 'sick', 'sick']
    >>> metadata = pd.Series(grouping, index=samples, name='status').to_frame()

    Now run ``ancombc`` to determine if there are any features that are significantly
    different in abundance between the healthy and the sick groups. Note that a
    pseudocount of 1 is manually added to the table to remove zero values.

    >>> result = ancombc(table + 1, metadata, formula='status').res
    >>> result.round(3)
                              Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F1        Intercept         -0.045  0.218 -0.207   0.836   1.000   False
              status[T.sick]     0.126  0.589  0.214   0.831   1.000   False
    F2        Intercept         -0.873  0.137 -6.381   0.000   0.000    True
              status[T.sick]     1.241  0.538  2.307   0.021   0.105   False
    F3        Intercept         -0.202  0.475 -0.425   0.671   1.000   False
              status[T.sick]     0.561  0.964  0.582   0.561   1.000   False
    F4        Intercept          1.005  0.136  7.399   0.000   0.000    True
              status[T.sick]    -1.723  0.392 -4.399   0.000   0.000    True
    F5        Intercept          0.141  0.422  0.335   0.738   1.000   False
              status[T.sick]    -0.690  0.625 -1.104   0.270   1.000   False
    F6        Intercept         -0.873  0.137 -6.381   0.000   0.000    True
              status[T.sick]     1.480  0.160  9.255   0.000   0.000    True
    F7        Intercept          0.157  0.401  0.391   0.695   1.000   False
              status[T.sick]     0.049  0.405  0.121   0.904   1.000   False

    The covariate "status[T.sick]" stands for the effect of the "sick" group relative
    to the reference group, "healthy" (the first group in alphabetical order is
    automatically selected as the reference group). "Log2(FC)" represents the
    :func:`clr`-transformed fold change of abundance (positive/negative: more/less
    abundant in "sick" than in "healthy", respectively). A "True" in the "Signif"
    column indicates a significantly differentially abundant feature-covariate pair.
    This example shows that two features differ by healthy/sick status.

    >>> result.query('Covariate != "Intercept" & Signif == True').round(3)
                              Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F4        status[T.sick]    -1.723  0.392 -4.399     0.0     0.0    True
    F6        status[T.sick]     1.480  0.160  9.255     0.0     0.0    True

    **An advanced example**

    Now we will create a complex dataset with 15 samples grouped into three disease
    status: "mild", "moderate" and "severe", plus age as a confounder. Our goal is
    to identify features that are differentially abundant between disease status.

    >>> samples = [f'S{i}' for i in range(1, 16)]
    >>> features = [f'F{i}' for i in range(1, 9)]
    >>> data = [[ 2,  0,  7,  0,  0,  2,  3,  2],
    ...         [ 0,  2,  0,  1,  1,  2,  0,  0],
    ...         [ 3,  1,  0,  9,  0,  1,  1,  0],
    ...         [ 2,  0,  1,  8,  1,  0, 11, 46],
    ...         [ 1,  0,  1,  1,  0,  0,  2,  2],
    ...         [ 0,  3, 22,  1,  1,  1,  3,  0],
    ...         [ 1,  7, 16,  1,  0,  0,  2,  2],
    ...         [ 0,  5,  6,  1,  2,  1,  0,  1],
    ...         [ 1,  7,  0,  2,  1,  0,  3,  2],
    ...         [ 0,  6,  4,  2,  0,  2,  2,  1],
    ...         [ 3, 13,  7,  0,  0,  0,  3,  9],
    ...         [ 1,  8,  5,  1,  0,  0,  0,  0],
    ...         [ 0,  5, 14,  1,  0,  1,  0,  1],
    ...         [ 5, 26,  3,  2,  0,  3,  1,  3],
    ...         [ 0, 18,  7,  0,  0,  3,  1,  0]]
    >>> table = pd.DataFrame(data, index=samples, columns=features)
    >>> status = ['mild'] * 5 + ['moderate'] * 5 + ['severe'] * 5
    >>> age = [39, 19, 20, 31, 15, 37, 27, 47, 26, 23, 39, 48, 46, 33, 36]
    >>> metadata = pd.DataFrame({'status': status, 'age': age}, index=samples)

    Run ``ancombc``. This time, we specify two factors: "status" and "age" in the
    formula, such that the function will test the individual effects of each factor
    while controlling for the other.

    >>> res = ancombc(
    ...     table + 1, metadata, formula='status + age', grouping='status')
    >>> res_main = res.res
    >>> res_main.round(3)
                                  Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F1        Intercept             -0.006  0.346 -0.016   0.987   1.000   False
              status[T.moderate]    -0.337  0.245 -1.376   0.169   1.000   False
              status[T.severe]       0.278  0.350  0.792   0.428   1.000   False
              age                    0.001  0.011  0.117   0.907   1.000   False
    F2        Intercept              0.072  0.446  0.160   0.873   1.000   False
              status[T.moderate]     1.906  0.259  7.371   0.000   0.000    True
              status[T.severe]       2.935  0.304  9.649   0.000   0.000    True
              age                   -0.022  0.013 -1.697   0.090   0.538   False
    F3        Intercept             -1.690  0.570 -2.967   0.003   0.021    True
              status[T.moderate]     1.010  0.578  1.747   0.081   0.565   False
              status[T.severe]       0.717  0.517  1.388   0.165   1.000   False
              age                    0.063  0.022  2.871   0.004   0.033    True
    F4        Intercept              0.439  0.483  0.910   0.363   1.000   False
              status[T.moderate]    -0.046  0.405 -0.113   0.910   1.000   False
              status[T.severe]      -0.244  0.536 -0.455   0.649   1.000   False
              age                   -0.004  0.019 -0.199   0.842   1.000   False
    F5        Intercept             -1.227  0.393 -3.125   0.002   0.014    True
              status[T.moderate]     0.274  0.293  0.934   0.350   1.000   False
              status[T.severe]      -0.323  0.320 -1.011   0.312   1.000   False
              age                    0.027  0.014  2.001   0.045   0.318   False
    F6        Intercept             -0.439  0.440 -0.999   0.318   1.000   False
              status[T.moderate]     0.114  0.432  0.264   0.792   1.000   False
              status[T.severe]       0.375  0.544  0.691   0.490   1.000   False
              age                    0.008  0.016  0.464   0.643   1.000   False
    F7        Intercept              0.201  0.472  0.426   0.670   1.000   False
              status[T.moderate]     0.082  0.302  0.270   0.787   1.000   False
              status[T.severe]      -0.264  0.401 -0.658   0.510   1.000   False
              age                    0.004  0.016  0.256   0.798   1.000   False
    F8        Intercept             -0.168  0.526 -0.319   0.750   1.000   False
              status[T.moderate]    -0.402  0.584 -0.688   0.491   1.000   False
              status[T.severe]      -0.299  0.728 -0.411   0.681   1.000   False
              age                    0.022  0.018  1.222   0.222   1.000   False

    We found that feature "F2" is significantly differentially (more) abundant in
    "moderate" and "severe" groups compared with "mild", which serves as the reference
    group. Additionally, feature "F3" is separately correlated with age.

    >>> res_main.query('Covariate != "Intercept" & Signif == True').round(3)
                                  Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F2        status[T.moderate]     1.906  0.259  7.371   0.000   0.000    True
              status[T.severe]       2.935  0.304  9.649   0.000   0.000    True
    F3        age                    0.063  0.022  2.871   0.004   0.033    True

    **Global test**

    Next, we will perform a *global test* to identify features that are differentially
    abundant between at least two status.

    >>> res_global = res.global_test()
    >>> res_global.round(3)
                    W  pvalue  qvalue  Signif
    FeatureID
    F1          5.934   0.103   0.617   False
    F2         95.473   0.000   0.000    True
    F3          3.210   0.402   1.000   False
    F4          0.503   0.445   1.000   False
    F5          6.839   0.065   0.458   False
    F6          0.631   0.541   1.000   False
    F7          1.044   0.813   1.000   False
    F8          0.548   0.479   1.000   False

    The global test result suggests that "F2" is differentially abundant between two of
    the three groups (though it doesn't tell you which groups).

    **Structural zero test**

    The structural zero test identifies features that are systematically absent from
    certain sample groups. This test is an option of the R command ``ancombc``. In
    scikit-bio, :func:`struc_zero` is a standalone function, as it is generally useful
    with or without ANCOM-BC.

    >>> from skbio.stats.composition import struc_zero
    >>> res_zero = struc_zero(table, metadata, 'status')
    >>> res_zero
         mild  moderate  severe
    F1  False     False   False
    F2  False     False   False
    F3  False     False   False
    F4  False     False   False
    F5  False     False    True
    F6  False     False   False
    F7  False     False   False
    F8  False     False   False

    The result reveals that feature "F5" is a structural zero in the "severe" groups,
    as all or most of its values are zero. Although the ANCOM-BC test itself didn't
    identify "F5", we should consider it as differentially (less) abundant than in the
    other two groups.

    We can use this additional information to update the global test result. The rule
    is that a feature should be considered as globally differentially abundant if it is
    a structural zero in at least one but not all groups.

    >>> signif_global = res_global['Signif'] | (
    ...     ~res_zero.all(axis=1) & res_zero.any(axis=1))
    >>> signif_global
    FeatureID
    F1    False
    F2     True
    F3    False
    F4    False
    F5     True
    F6    False
    F7    False
    F8    False
    dtype: bool

    The main ANCOM-BC result can also be updated with structural zeros.

    >>> signif_main = res_main.query(
    ...     'Covariate.str.startswith("status[T.")')['Signif'].unstack()
    >>> signif_main.columns = signif_main.columns.str.removeprefix(
    ...     'status[T.').str.removesuffix(']')
    >>> signif_zero = res_zero.loc[signif_main.index, signif_main.columns]
    >>> signif_main |= signif_zero
    >>> signif_main.columns.name = None
    >>> signif_main
               moderate  severe
    FeatureID
    F1            False   False
    F2             True    True
    F3            False   False
    F4            False   False
    F5            False    True
    F6            False   False
    F7            False   False
    F8            False   False

    """
    return _ancombc_core(
        table=table,
        metadata=metadata,
        v2=False,
        formula=formula,
        grouping=grouping,
        max_iter=max_iter,
        tol=tol,
        alpha=alpha,
        p_adjust=p_adjust,
    )


def ancombc2(
    table,
    metadata,
    formula,
    grouping=None,
    pseudocount=0,
    aggregator=None,
    s0_perc=0.05,
    max_iter=100,
    tol=1e-5,
    p_adjust="holm",
    alpha=0.05,
):
    r"""Perform differential abundance test using ANCOM-BC2.

    Analysis of Compositions of Microbiomes with Bias Correction 2
    (ANCOM-BC2) [1]_ is a differential abundance testing method that extends
    ANCOM-BC with several improvements:

    - **SAM-like fudge factor** (``s0_perc``): Adds a small constant to
      the denominator of the test statistic to stabilize inference for
      rare taxa with extremely small standard errors.
    - **Pseudo-count sensitivity analysis**: Assesses robustness of results
      to the choice of pseudo-count added to zero counts.
    - **Multi-group comparisons**: Supports global test, pairwise
      directional test, Dunnett's type of test, and trend test.
    - **Mixed directional FDR** (mdFDR): Controls false discovery rate
      in multi-group settings using a bootstrap-based approach.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        Matrix containing count or proportional abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    metadata : pd.DataFrame or 2-D array_like
        Metadata of the samples. Rows correspond to samples and columns correspond
        to covariates (attributes). Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    formula : str or generic Formula object
        Formula defining the model using factors included in the metadata columns.
        Refer to `Patsy's documentation
        <https://patsy.readthedocs.io/en/latest/formulas.html>`_ on how to specify
        a formula.
    grouping : str, optional
        Metadata column defining sample groups for post-hoc analyses. Must be a term
        in ``formula`` and contain at least three groups. Other terms are treated as
        adjustment covariates. It has no impact on the main esult. If None (default),
        post-hoc analyses will be unavailable.
    pseudocount : int or float, optional
        Pseudocount to add to all abundance data. Default is 0.
    aggregator : callable, mapping, or 1-D array_like, optional
        Rule for aggregating features before final regression. Can be a function or a
        dictionary that maps each feature ID to an aggregate ID, or a plain list or
        array of aggregate ID per feature in table order. By default, no aggregation
        is performed.
    s0_perc : float, optional
        SAM-like fudge factor percentile. Default is 0.05.
    max_iter : int, optional
        Maximum number of iterations for the bias estimation process. Default is 100.
    tol : float, optional
        Absolute convergence tolerance for the bias estimation process. Default is
        1e-5.
    alpha : float, optional
        Significance level for the statistical tests. Must be in the range of (0, 1).
        Default is 0.05.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are
        Holm-Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-Hochberg
        ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported by
        statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.

    Returns
    -------
    :class:`ANCOMBCResult`
        Result object with primary results and post-hoc analysis methods.

    See Also
    --------
    ancombc : ANCOM-BC without sampling fraction correction and multi-group testing.
    struc_zero : Standalone structural zero detection.

    References
    ----------
    .. [1] Lin, H., & Peddada, S. D. (2024). Multigroup analysis of compositions of
       microbiomes with covariate adjustments and repeated measures. Nature Methods,
       21(1), 83-91.

    Examples
    --------
    >>> from skbio.stats.composition import ancombc2
    >>> import pandas as pd

    >>> samples = [f'S{i}' for i in range(1, 16)]
    >>> features = [f'F{i}' for i in range(1, 9)]
    >>> data = [[ 2,  0,  7,  0,  0,  2,  3,  2],
    ...         [ 0,  2,  0,  1,  1,  2,  0,  0],
    ...         [ 3,  1,  0,  9,  0,  1,  1,  0],
    ...         [ 2,  0,  1,  8,  1,  0, 11, 46],
    ...         [ 1,  0,  1,  1,  0,  0,  2,  2],
    ...         [ 0,  3, 22,  1,  1,  1,  3,  0],
    ...         [ 1,  7, 16,  1,  0,  0,  2,  2],
    ...         [ 0,  5,  6,  1,  2,  1,  0,  1],
    ...         [ 1,  7,  0,  2,  1,  0,  3,  2],
    ...         [ 0,  6,  4,  2,  0,  2,  2,  1],
    ...         [ 3, 13,  7,  0,  0,  0,  3,  9],
    ...         [ 1,  8,  5,  1,  0,  0,  0,  0],
    ...         [ 0,  5, 14,  1,  0,  1,  0,  1],
    ...         [ 5, 26,  3,  2,  0,  3,  1,  3],
    ...         [ 0, 18,  7,  0,  0,  3,  1,  0]]
    >>> table = pd.DataFrame(data, index=samples, columns=features)
    >>> status = ['mild'] * 5 + ['moderate'] * 5 + ['severe'] * 5
    >>> age = [39, 19, 20, 31, 15, 37, 27, 47, 26, 23, 39, 48, 46, 33, 36]
    >>> metadata = pd.DataFrame({'status': status, 'age': age}, index=samples)

    >>> res = ancombc2(table + 1, metadata, formula='status + age')
    >>> res_main = res.res
    >>> res_main.round(3)
                                  Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F1        Intercept              0.342  0.723  0.473   0.645   1.000   False
              status[T.moderate]    -0.337  0.517 -0.652   0.528   1.000   False
              status[T.severe]       0.278  0.673  0.413   0.688   1.000   False
              age                    0.001  0.024  0.054   0.958   1.000   False
    F2        Intercept             -0.541  0.796 -0.680   0.511   1.000   False
              status[T.moderate]     1.906  0.527  3.618   0.004   0.032    True
              status[T.severe]       2.935  0.639  4.590   0.001   0.006    True
              age                   -0.022  0.025 -0.877   0.399   1.000   False
    F3        Intercept             -2.242  0.893 -2.509   0.029   0.232   False
              status[T.moderate]     1.010  0.788  1.282   0.226   1.000   False
              status[T.severe]       0.717  0.803  0.893   0.391   1.000   False
              age                    0.063  0.032  1.955   0.077   0.612   False
    F4        Intercept              0.580  0.824  0.703   0.497   1.000   False
              status[T.moderate]    -0.046  0.639 -0.071   0.945   1.000   False
              status[T.severe]      -0.244  0.819 -0.298   0.771   1.000   False
              age                   -0.004  0.029 -0.126   0.902   1.000   False
    F5        Intercept             -0.502  0.757 -0.663   0.521   1.000   False
              status[T.moderate]     0.274  0.552  0.496   0.630   1.000   False
              status[T.severe]      -0.323  0.651 -0.497   0.629   1.000   False
              age                    0.027  0.025  1.069   0.308   1.000   False
    F6        Intercept             -0.045  0.791 -0.057   0.956   1.000   False
              status[T.moderate]     0.114  0.662  0.172   0.866   1.000   False
              status[T.severe]       0.375  0.825  0.455   0.658   1.000   False
              age                    0.008  0.028  0.275   0.788   1.000   False
    F7        Intercept              0.291  0.816  0.357   0.728   1.000   False
              status[T.moderate]     0.082  0.559  0.146   0.887   1.000   False
              status[T.severe]      -0.264  0.711 -0.371   0.718   1.000   False
              age                    0.004  0.027  0.150   0.883   1.000   False
    F8        Intercept             -0.118  0.858 -0.138   0.893   1.000   False
              status[T.moderate]    -0.402  0.793 -0.507   0.622   1.000   False
              status[T.severe]      -0.299  0.984 -0.304   0.767   1.000   False
              age                    0.022  0.029  0.763   0.462   1.000   False

    """
    return _ancombc_core(
        table=table,
        metadata=metadata,
        v2=True,
        formula=formula,
        grouping=grouping,
        max_iter=max_iter,
        tol=tol,
        alpha=alpha,
        p_adjust=p_adjust,
        pseudo=pseudocount,
        s0_perc=s0_perc,
        aggregator=aggregator,
    )


def _ancombc_core(
    table,
    metadata,
    v2=False,
    formula=None,
    grouping=None,
    aggregator=None,
    p_adjust="holm",
    pseudo=0,
    s0_perc=0.05,
    alpha=0.05,
    max_iter=100,
    tol=1e-5,
):
    """ANCOM-BC/BC2 core function."""
    matrix, samples, features = _ingest_table(table)

    # ANCOM-BC does not handle zeros in the input table. The user should have added a
    # pseudocount. ANCOM-BC2 can handle zeros.
    _check_composition(np, matrix, nozero=not v2)

    # Validate metadata and cast to numbers where applicable.
    metadata = _check_metadata(metadata, matrix, samples)
    metadata = _type_cast_to_float(metadata)

    # Create a design matrix based on metadata and formula.
    dmat = dmatrix(formula, metadata)

    # Obtain a list of covariates by selecting the relevant columns.
    covars = dmat.design_info.column_names
    n_covars = len(covars)

    # Validate sample grouping for post-hoc analysis
    # When provided, only this coefficient covariance submatrix is retained.
    groups = _validate_grouping(metadata, dmat, grouping)

    # Validate parameters
    if not 0 < alpha < 1:
        raise ValueError(f"`alpha`={alpha} is not within 0 and 1.")
    if not pseudo >= 0:
        raise ValueError(f"Pseudocount must be a non-negative number.")

    # Transform data
    data, missing = _transform_data(matrix, pseudo, center=v2)

    # Estimate initial model parameters.
    # ANCOM-BC2 always re-estimates parameters after sampling-fraction correction, so
    # its initial EM fit needs variances only. For ANCOM-BC this is the final fit, so
    # retain the grouping covariance submatrix when post-hoc analysis was requested.
    # Data matrix will be overwritten during dense parameter estimation. However it is
    # needed for sampling fraction estimation in the ANCOM-BC2 route. Therefore `v2`
    # instructs the function to make a copy.
    var_hat, beta, _, vcov_hat = _estimate_params(
        data, dmat, None if v2 else groups, missing, v2
    )

    # Estimate and correct for sampling bias via expectation-maximization (EM).
    # beta: (n_covariates, n_features); iterate over covariates (rows).
    bias = np.empty((n_covars, 3))
    for i in range(n_covars):
        bias[i] = _estimate_bias_em(beta[i], var_hat[:, i], tol=tol, max_iter=max_iter)

    delta_em, delta_wls, var_delta = bias[:, 0], bias[:, 1], bias[:, 2]

    # ANCOM-BC (original)
    if not v2:
        # Correct coefficients (logFC) according to estimated bias.
        beta_hat = beta.T - delta_em

        # Skip degree of freedom calculation.
        dof = None

    # ANCOM-BC2
    else:
        # Estimate sampling fractions
        theta_hat = _sample_fractions(data, dmat, beta, delta_em, missing=missing)

        # Aggregate data
        if aggregator is not None:
            matrix, features = _aggregate_features(matrix, aggregator, features)
            data, missing = _transform_data(matrix, pseudo, center=v2)

        # Correct data for sampling fractions
        data -= theta_hat[:, None]

        # Re-estimate parameters.
        # Since this is the final fit, retain the grouping covariance submatrix if
        # requested.
        var_hat, beta_hat, _, vcov_hat = _estimate_params(
            data, dmat, groups, missing, False
        )
        beta_hat = beta_hat.T

        # Adjust variances
        _adjust_variances(var_hat, vcov_hat, var_delta, s0_perc, groups)

        # Compute per-feature degree of freedom (observed samples - covariates)
        if missing is not None:
            n_valid = np.sum(~missing, axis=0)
            dof = np.where(n_valid > n_covars, n_valid - n_covars, np.nan)
        else:
            n_samps = matrix.shape[0]
            dof = n_samps - n_covars if n_samps > n_covars else np.nan

    # Output primary results. Compute statistics and populate the DataFrame
    # incrementally to minimize the number of feature-by-covariate arrays alive at
    # once. This is more memory-efficient than calculating all statistics first and
    # constructing the DataFrame from repeated Python label lists.
    if features is None:
        features = np.arange(matrix.shape[1])
    res = _format_results(beta_hat, var_hat, features, covars, alpha, p_adjust, dof)

    method = "ANCOM-BC" if not v2 else "ANCOM-BC2"

    return ANCOMBCResult(
        res=res,
        method=method,
        _dmat=dmat,
        _beta_hat=beta_hat,
        _var_hat=var_hat,
        _vcov_hat=vcov_hat,
        _grouping=grouping,
        _group_indices=groups,
        _dof=dof,
        _features=features,
        _covariates=covars,
        _p_adjust=p_adjust,
        _alpha=alpha,
    )


def _validate_grouping(metadata, dmat, grouping):
    """Validate a post-hoc grouping and return its design-matrix column indices."""
    if grouping is None:
        return None
    if not isinstance(grouping, str):
        raise TypeError("`grouping` must be a metadata column name or None.")
    if grouping not in metadata.columns:
        raise ValueError(f"`grouping`={grouping!r} is not a metadata column.")

    term_slices = dmat.design_info.term_name_slices
    if grouping not in term_slices:
        terms = [name for name in term_slices if name != "Intercept"]
        raise ValueError(
            f"`grouping`={grouping!r} must be a term in `formula`. "
            f"Available terms are: {terms}."
        )

    n_groups = metadata[grouping].nunique(dropna=True)
    if n_groups < 3:
        raise ValueError(
            f"`grouping`={grouping!r} must contain at least three observed groups "
            f"for post-hoc analysis; found {n_groups}."
        )

    s = term_slices[grouping]
    indices = np.arange(s.start, s.stop, dtype=int)
    if indices.size < 2:
        raise ValueError(
            f"`grouping`={grouping!r} does not produce at least two group "
            "coefficients in the design matrix. Ensure it is modeled as a "
            "categorical factor."
        )
    return indices


def _transform_data(data, pseudo=None, center=False):
    """Transform abundance data table.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Data table.
    pseudo : float or int, optional
        Pseudocount.
    center : bool, optional
        If True, center data at 0 for each feature.

    Returns
    -------
    data : ndarray of shape (n_samples, n_features)
        Transformed data table.
    missing : ndarray of shape (n_samples, n_features) or None
        Boolean mask of zero values in the table, or None if the table is zero-free or
        a pseudocount is added.

    Notes
    -----
    For ANCOM-BC, a simple log-transform is applied. For ANCOM-BC2, the log-transformed
    data is centered at 0 per feature. This is equivalent to `rclr(data, axis=0)`, but
    more efficient.

    The returned `data` is a new data table in any float type. The original data table
    is kept intact.

    """
    ### Step 1: Identify zeros values and/or add pseudocount. ###

    # Add pseudocount if provided. Output must be in float type.
    if pseudo:
        data = data + float(pseudo)
        missing = None

    # Otherwise, identify zero values and make a float copy of data.
    else:
        missing = data == 0
        if not np.any(missing):
            missing = None
        if np.issubdtype(data.dtype, np.floating):
            data = data.copy()
        else:
            data = data.astype(float)

    ### Step 2: Log-transform observed data and optionally center ###

    # Fill missing values with 1 for safe log (they become 0)
    if missing is not None:
        data[missing] = 1.0

    # Log-transform in place
    np.log(data, out=data)

    # Center at 0 (subtract by feature mean) in place
    if center:
        if missing is None:
            mean_ = np.mean(data, axis=0, keepdims=True)
        else:
            n_obs = data.shape[0] - np.sum(missing, axis=0, keepdims=True)
            # n_obs = np.where(n_obs > 0, n_obs, 1)
            mean_ = np.sum(data, axis=0, keepdims=True)
            mean_ /= n_obs
        data -= mean_

    # Fill missing values with NaN
    # NOTE: This is not needed except for `_sample_fractions`
    if missing is not None:
        data[missing] = np.nan

    return data, missing


def _estimate_params(data, dmat, groups, missing, keep_data=True):
    """Estimate model parameters.

    In ANCOM-BC, this function is used for initial estimation of model parameters
    (coefficients, variances and mean residuals) based on the observed data.

    In ANCOM-BC2, this function is additionally used for final parameter estimation on
    data corrected for sampling fractions.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Transformed data table.
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.
    groups : ndarray of shape (n_groups - 1,), bool, or None
        Indices of coefficients in the design matrix whose covariance submatrix will
        be calculated. If True, the entire covariance matrix will be calculated. If
        None, no covariance will be returned.
    missing : ndarray of shape (n_samples, n_features), or None
        Boolean mask of zero values in the pre-transformation data table, or None if
        the table is zero-free or a pseudocount was applied.
    keep_data : bool, optional
        If True, the original data table will be kept intact. Relevant when `missing`
        is not provided (dense route). Sparse route always keeps the data table.

    Returns
    -------
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances of regression coefficients.
    beta : ndarray of shape (n_covariates, n_features)
        Estimated coefficients (log-fold changes).
    theta : ndarray of shape (n_samples,)
        Per-sample mean residuals of estimated data.
    covmat : ndarray or None
        Estimated covariance matrices of coefficients.
        - When `groups=None`: None.
        - When `groups=True`: Full covariance matrices.
          Shape: (n_features, n_covariates, n_covariates)
        - When `groups` is ndarray of indices: Covariance submatrices of the indices.
          Shape: (n_features, n_groups - 1, n_groups - 1)

    Notes
    -----
    Returned `var_hat` is F-contiguous, which should (presumably) benefit the
    downstream EM optimization process.
    """
    dmat = np.asarray(dmat, dtype=data.dtype)
    if missing is not None:
        return _estimate_params_sparse(data, dmat, missing, groups)
    elif keep_data:
        return _estimate_params_dense(data.copy(), dmat, groups)
    else:
        return _estimate_params_dense(data, dmat, groups)


def _estimate_params_dense(data, dmat, groups=True):
    """Estimate model parameters from a dense matrix (no missing values).

    The original R code performs iterative maximum likelihood estimation to calculate
    the coefficients, and calls MASS:ginv to calculate the pseudo inverse Gram matrix
    using the Moore-Penrose method. We noticed that the former can be replaced with
    ordinary least squares by NumPy's `lstsq`; and the latter matches the method by
    NumPy's `pinv`. We further noticed that both functions perform singular value
    decomposition (SVD). Therefore, the following code only performs SVD once, and
    uses the intermediates to calculate coefficients and inverse Gram matrix.

    NOTE: This function overwrites the data table.
    """
    # Calculate coefficients using least-squares fitting
    dmat_inv, _ = _lstsq_fit(dmat)
    beta = dmat_inv @ data

    # Calculate residuals = data - fitted data
    # data -= dmat @ beta
    # The process below is equivalent but saves memory by avoiding materializing the
    # entire fitted matrix. TODO: Revisit.
    _calc_residual(data, dmat, beta)

    # Calculate per-sample mean residuals (theta)
    theta = np.mean(data, axis=1, keepdims=True)

    # Center and square residuals
    data -= theta  # eps
    np.square(data, out=data)  # eps2

    # Calculate variances and covariance matrix of coefficients
    var_hat, covmat = _calc_grouped_var_cov(data, dmat_inv.T, groups)

    return var_hat, beta, theta.reshape(-1), covmat


def _estimate_params_sparse(
    data,
    dmat,
    missing,
    groups=True,
    tol=1e-2,
    max_iter=20,
    direct=False,
    batched=True,
    batch=None,
):
    """Estimate model parameters from a sparse matrix (with missing values).

    This mirrors the fixed-effects ``.iter_mle`` path used by ANCOM-BC2 when
    ``theta`` is initially ``NULL``. The input ``data`` is expected to already
    be log-transformed and centered, with zeros represented by ``NaN``.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Zero-handled, log-transformed and centered data. Missing/zero entries
        must be represented by NaN.
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix. ANCOM-BC2 requires the fixed-effect design to be fully
        observed for this path.
    missing : ndarray of shape (n_samples, n_features), optional
        Boolean mask of zero/missing values in ``data`` before transformation.
    tol : float, default=1e-2
        Iteration tolerance. The ANCOM-BC2 default is 1e-2.
    max_iter : int, default=20
        Maximum number of iterations. The ANCOM-BC2 default is 20.
    direct : bool, default=False
        Solve the coupled sample effects directly instead of reproducing the
        ANCOM-BC2 iteration.
    groups : 1-D array_like of int, optional
        Coefficient indices whose covariance submatrix should be calculated when
        ``full_covariance=False``. All coefficient variances are always calculated.
        If None, no covariance matrix is returned.
    batched : bool, optional
        Missing-response least-squares implementation. False materializes the
        complete feature-by-sample pseudoinverse tensor and is retained for numerical
        comparison. True (default) performs the masked SVD in feature blocks and stores
        only compact per-feature spectral operators.
    batch : int, optional
        Number of features per masked-SVD block for ``solver="batched"``. By default a
        size is selected to keep the principal SVD input/output workspaces near 32 MiB.

    Returns
    -------
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances of regression coefficients.
    beta : ndarray of shape (n_covariates, n_features)
        Estimated regression coefficients.
    theta : ndarray of shape (n_samples,)
        Per-sample mean residuals / initial sampling-fraction estimates.
    covmat : ndarray or None
        Estimated covariance matrices of coefficients.

    """
    func = _lstsq_nan_batch if batched else _lstsq_nan_full
    theta, beta = func(data, dmat, missing, tol, max_iter, direct, batch)

    # Residuals used by the sandwich covariance estimator. Missing residuals are set to
    # zero here and handled separately below to preserve ANCOM-BC2/R's 0.1 correction.
    intm = np.empty_like(data)
    np.matmul(dmat, beta.T, out=intm)
    intm *= -1.0
    intm += data
    intm -= theta[:, None]
    np.square(intm, out=intm)
    intm[missing] = 0.0

    # Important R detail: .iter_mle uses ONE global ginv(X.T @ X), not a
    # feature-specific inverse based on each feature's observed rows. Its transpose is
    # the influence matrix H used by the GEMM covariance calculation below.
    dmat_inv, gram_sum = _lstsq_fit(dmat, gram=True)
    n_missing = missing.sum(axis=0)

    # Calculate variances and covariance matrix of coefficients
    var_hat, covmat = _calc_grouped_var_cov(intm, dmat_inv.T, groups)

    # This code matches R, which replaces every NaN element of resid^2 * x x^T with
    # 0.1. Since dmat is finite, each missing residual contributes a matrix filled with
    # 0.1 before transformation.
    if np.any(n_missing):
        var_hat += (0.1 * n_missing)[:, None] * (gram_sum * gram_sum)[None, :]
        if covmat is not None:
            if groups is not True:
                gram_sum = gram_sum[groups]
            missing_cov = np.outer(gram_sum, gram_sum)
            covmat += (0.1 * n_missing)[:, None, None] * missing_cov

    # Ensure the retained covariance diagonal is exactly the same array of variances
    # used by the primary analysis, including the R missing-value correction.
    if covmat is not None:
        if groups is True:
            groups = np.arange(var_hat.shape[1])
        diag = np.arange(groups.size)
        covmat[:, diag, diag] = var_hat[:, groups]

    return var_hat, beta.T.copy(), theta, covmat


def _lstsq_nan_full(data, dmat, missing, tol, max_iter, direct, batch=None):
    """Fit missing-response models using the original full batched pseudoinverse.

    This helper intentionally preserves the pre-optimization route so the compact
    implementation can be regression-tested against it. Its dominant arrays have shape
    `(n_features, n_samples, n_covariates)`.
    """
    n_samps = data.shape[0]

    W = 1.0 - missing.T
    X_w = dmat[None, :, :] * W[:, :, None]
    U, S, Vh = np.linalg.svd(X_w, full_matrices=False)
    S_inv = _invert_singular(S)
    V = np.swapaxes(Vh, -1, -2)
    dmat_inv = np.einsum("fpk,fk,fsk->fps", V, S_inv, U, optimize=True)

    intm = np.zeros_like(data)
    np.copyto(intm, data, where=~missing)
    theta = np.zeros(n_samps, dtype=data.dtype)
    beta = np.einsum("fps,sf->fp", dmat_inv, intm, optimize=True)

    if direct:
        theta, beta = _estimate_params_nan_direct(
            data, dmat, missing, W, dmat_inv, intm, beta
        )
    else:
        theta, beta = _estimate_params_nan_iter(
            data, dmat, missing, dmat_inv, intm, theta, beta, tol, max_iter
        )
    return theta, beta


def _lstsq_nan_batch(
    data,
    dmat,
    missing,
    tol,
    max_iter,
    direct,
    batch=None,
    compact_cond_max=1e4,
):
    """Fit missing-response models with chunked SVD and compact spectral operators.

    Each feature still receives an SVD of its masked design matrix, using the exact same
    singular-value cutoff as the legacy route. Only ``Vh`` and inverse singular values
    are retained. The large ``U`` and masked-design tensors exist for one feature block
    at a time, reducing persistent storage from O(F*N*P) to O(F*P^2).
    """
    n_samples, n_features = data.shape
    n_covariates = dmat.shape[1]
    n_components = min(n_samples, n_covariates)

    # The legacy route's `1.0 - missing` promotes masked designs to float64. Preserve
    # that behavior even when `dmat` itself is float32.
    operator_dtype = np.result_type(dmat.dtype, np.dtype(float))
    Vh_all = np.empty((n_features, n_components, n_covariates), dtype=operator_dtype)
    S_inv_all = np.empty((n_features, n_components), dtype=operator_dtype)

    if batch is None:
        # X_w and U are both approximately batch*N*P floats, while W contributes another
        # batch*N floats. Keep these principal workspaces near 32 MiB; LAPACK may use
        # additional internal workspace, so this is intentionally conservative.
        bytes_per_feature = (
            operator_dtype.itemsize * n_samples * (n_covariates + n_components + 1)
        )
        batch = max(1, (32 << 20) // max(bytes_per_feature, 1))
        batch = min(batch, n_features)
    elif batch < 1:
        raise ValueError("`svd_batch_size` must be a positive integer.")

    dmat_work = np.asarray(dmat, dtype=operator_dtype)
    fallback_blocks = {}
    for start in range(0, n_features, batch):
        stop = min(start + batch, n_features)
        # Preserve legacy dtype and SVD input exactly, but only for this feature block.
        W_block = 1.0 - missing[:, start:stop].T
        X_w = dmat_work[None, :, :] * W_block[:, :, None]
        U, S, Vh = np.linalg.svd(X_w, full_matrices=False)
        S_inv = _invert_singular(S)
        Vh_all[start:stop] = Vh
        S_inv_all[start:stop] = S_inv

        # Applying a compact operator through X.T @ y is algebraically equivalent to
        # X^+ @ y, but it behaves like a normal-equation calculation numerically and
        # can lose accuracy for severely ill-conditioned masked designs. Retain the
        # full, stable SVD pseudoinverse only for those exceptional features. This keeps
        # the common case O(F*P^2) while preserving legacy behavior where it matters.
        cutoff = 1e-15 * np.max(S, axis=1, keepdims=True)
        retained = S > cutoff
        min_retained = np.min(np.where(retained, S, np.inf), axis=1)
        max_s = np.max(S, axis=1)
        cond = np.divide(
            max_s,
            min_retained,
            out=np.ones_like(max_s),
            where=np.isfinite(min_retained) & (min_retained > 0),
        )
        unstable = cond > compact_cond_max
        if np.any(unstable):
            local = np.flatnonzero(unstable)
            V = np.swapaxes(Vh[local], -1, -2)
            pinv = np.einsum(
                "fpk,fk,fsk->fps", V, S_inv[local], U[local], optimize=True
            )
            fallback_blocks[start] = (local, pinv)

        # U and the masked design are the large block-local arrays.
        del U, X_w, W_block

    # Preserve the legacy route's response-workspace precision. In particular, a
    # float32 transformed response keeps float32 residual/theta updates even though the
    # masked design SVD itself is float64. The compact RHS/beta calculations still use
    # the promoted operator dtype, matching dmat_inv @ intm in the legacy route.
    intm = np.zeros_like(data)
    np.copyto(intm, data, where=~missing)
    theta = np.zeros(n_samples, dtype=data.dtype)

    solve_dtype = np.result_type(data.dtype, operator_dtype)
    rhs = np.empty((n_covariates, n_features), dtype=solve_dtype)
    beta = np.empty((n_features, n_covariates), dtype=solve_dtype)
    _apply_masked_spectral_operator(
        dmat_work,
        intm,
        Vh_all,
        S_inv_all,
        rhs=rhs,
        out=beta,
        fallback_blocks=fallback_blocks,
    )

    if direct:
        theta, beta = _estimate_params_nan_direct_batched(
            data,
            dmat_work,
            missing,
            Vh_all,
            S_inv_all,
            intm,
            beta,
            rhs,
            batch,
            fallback_blocks,
        )
    else:
        theta, beta = _estimate_params_nan_iter_batched(
            data,
            dmat_work,
            missing,
            Vh_all,
            S_inv_all,
            intm,
            theta,
            beta,
            rhs,
            tol,
            max_iter,
            fallback_blocks,
        )
    return theta, beta


def _apply_masked_spectral_operator(
    dmat,
    adjusted,
    Vh,
    S_inv,
    rhs=None,
    out=None,
    fallback_blocks=None,
):
    """Apply compact operators, with stable SVD fallbacks when needed."""
    n_features = S_inv.shape[0]
    n_covariates = Vh.shape[2]
    dtype = np.result_type(dmat.dtype, adjusted.dtype, Vh.dtype, S_inv.dtype)
    if rhs is None:
        rhs = np.empty((n_covariates, n_features), dtype=dtype)
    if out is None:
        out = np.empty((n_features, n_covariates), dtype=dtype)

    # Missing adjusted responses are already zero, so X.T @ adjusted is simultaneously
    # the RHS for every feature-specific masked regression.
    np.matmul(dmat.T, adjusted, out=rhs)

    # X_f = U S V.T and rhs = X_f.T y = V S U.T y. Therefore
    # X_f^+ y = V S^-2 V.T rhs. Retaining V and S instead of pinv(X.T X) preserves the
    # legacy SVD's rank decision and avoids forming normal equations.
    tmp = np.einsum("fkp,pf->fk", Vh, rhs, optimize=True)
    tmp *= S_inv * S_inv
    np.einsum("fpk,fk->fp", np.swapaxes(Vh, -1, -2), tmp, out=out, optimize=True)

    if fallback_blocks:
        for start, (local, pinv) in fallback_blocks.items():
            indices = start + local
            out[indices] = np.einsum(
                "kps,sk->kp", pinv, adjusted[:, indices], optimize=True
            )
    return out


def _estimate_params_nan_iter_batched(
    data,
    dmat,
    missing,
    Vh,
    S_inv,
    intm,
    theta,
    beta,
    rhs,
    tol,
    max_iter,
    fallback_blocks,
):
    """Reproduce ANCOM-BC2 iterations using compact feature operators."""
    epsilon = 100.0
    it = 0
    beta_new = np.empty_like(beta)
    while epsilon > tol and it < max_iter:
        np.subtract(data, theta[:, None], out=intm)
        np.copyto(intm, 0.0, where=missing)
        _apply_masked_spectral_operator(
            dmat,
            intm,
            Vh,
            S_inv,
            rhs=rhs,
            out=beta_new,
            fallback_blocks=fallback_blocks,
        )

        np.matmul(dmat, beta_new.T, out=intm)
        np.copyto(intm, 0.0, where=missing)
        np.subtract(data, intm, out=intm)
        theta_new = np.nanmean(intm, axis=1)

        epsilon = np.sqrt(
            np.nansum((beta_new - beta) ** 2) + np.nansum((theta_new - theta) ** 2)
        )
        beta, beta_new = beta_new, beta
        theta = theta_new
        it += 1
    return theta, beta


def _estimate_params_nan_direct_batched(
    data,
    dmat,
    missing,
    Vh,
    S_inv,
    intm,
    beta,
    rhs,
    svd_batch_size,
    fallback_blocks,
):
    """Direct fixed-point solve using compact operators and bounded reconstruction."""
    n_samples, n_features = data.shape
    n_covariates = dmat.shape[1]

    np.matmul(dmat, beta.T, out=intm)
    np.copyto(intm, 0.0, where=missing)
    np.subtract(data, intm, out=intm)
    np.nan_to_num(intm, copy=False, nan=0.0)
    residual_sum = intm.sum(axis=1)
    observed_counts = (~missing).sum(axis=1)

    system = np.diag(observed_counts.astype(np.result_type(dmat, float), copy=False))
    V = np.swapaxes(Vh, -1, -2)
    for start in range(0, n_features, svd_batch_size):
        stop = min(start + svd_batch_size, n_features)
        W_block = 1.0 - missing[:, start:stop].T

        # Reconstruct only this block of X_f^+ = V S^-2 V.T X_f.T. This is needed by
        # the direct sample-effect system but is discarded immediately afterwards.
        Vh_block = Vh[start:stop]
        V_block = V[start:stop]
        scaled = V_block * (S_inv[start:stop] ** 2)[:, None, :]
        gram_inv = np.matmul(scaled, Vh_block)
        xtw = dmat.T[None, :, :] * W_block[:, None, :]
        dmat_inv_block = np.matmul(gram_inv, xtw)
        if start in fallback_blocks:
            local, pinv = fallback_blocks[start]
            dmat_inv_block[local] = pinv
        system -= np.einsum(
            "bs,sp,bpt->st", W_block, dmat, dmat_inv_block, optimize=True
        )

    constraint = dmat.T * observed_counts
    augmented = np.block(
        [[system, constraint.T], [constraint, np.zeros((n_covariates,) * 2)]]
    )
    theta = np.linalg.lstsq(
        augmented,
        np.concatenate((residual_sum, np.zeros(n_covariates))),
        rcond=None,
    )[0][:n_samples]

    np.subtract(data, theta[:, None], out=intm)
    np.copyto(intm, 0.0, where=missing)
    _apply_masked_spectral_operator(
        dmat,
        intm,
        Vh,
        S_inv,
        rhs=rhs,
        out=beta,
        fallback_blocks=fallback_blocks,
    )
    return theta, beta


def _estimate_params_nan_iter(
    data, dmat, missing, dmat_inv, intm, theta, beta, tol, max_iter
):
    """Reproduce ANCOM-BC2's alternating beta/theta updates."""
    # Missing responses couple features through theta, so these are not
    # independent OLS fits.
    epsilon = 100.0
    it = 0
    while epsilon > tol and it < max_iter:
        np.subtract(data, theta[:, None], out=intm)
        np.copyto(intm, 0.0, where=missing)  # adjusted
        beta_new = np.einsum("fps,sf->fp", dmat_inv, intm, optimize=True)

        # R initializes each fitted vector with zero and writes fitted values
        # only at response rows used by lm().
        np.matmul(dmat, beta_new.T, out=intm)
        np.copyto(intm, 0.0, where=missing)
        np.subtract(data, intm, out=intm)  # fitted
        theta_new = np.nanmean(intm, axis=1)

        epsilon = np.sqrt(
            np.nansum((beta_new - beta) ** 2) + np.nansum((theta_new - theta) ** 2)
        )
        beta = beta_new
        theta = theta_new
        it += 1
    return theta, beta


def _estimate_params_nan_direct(data, dmat, missing, W, dmat_inv, intm, beta):
    """Solve the fully converged fixed point for zero-containing data.

    This exact method should produce numerically better result than the iterative
    optimization method. But it costs more compute due to solving a cubic system.

    """
    # The theta update is affine: theta_new = D^-1 (r + S @ theta).
    # Solve its fixed point while imposing the constraint selected by the
    # zero-initialized iteration, X.T @ D @ theta = 0.
    fitted = np.matmul(dmat, beta.T)
    np.copyto(fitted, 0.0, where=missing)
    np.subtract(data, fitted, out=intm)
    np.nan_to_num(intm, copy=False, nan=0.0)
    residual_sum = intm.sum(axis=1)
    observed_counts = (~missing).sum(axis=1)
    system = np.diag(observed_counts) - np.einsum(
        "fs,sp,fpt->st", W, dmat, dmat_inv, optimize=True
    )
    constraint = dmat.T * observed_counts
    augmented = np.block(
        [[system, constraint.T], [constraint, np.zeros((dmat.shape[1],) * 2)]]
    )
    theta = np.linalg.lstsq(
        augmented,
        np.concatenate((residual_sum, np.zeros(dmat.shape[1]))),
        rcond=None,
    )[0][: data.shape[0]]
    np.subtract(data, theta[:, None], out=intm)
    np.copyto(intm, 0.0, where=missing)
    beta = np.einsum("fps,sf->fp", dmat_inv, intm, optimize=True)
    return theta, beta


def _calc_residual(data, dmat, beta, target_bytes=8 << 20):
    """Calculate regression residuals in-place.

    ``data -= dmat @ beta`` would allocate a complete fitted matrix. Multiplying
    feature blocks into an approximately 8 MiB temporary bounds that extra memory while
    retaining efficient BLAS matrix multiplication.
    """
    n_samps, n_feats = data.shape
    if n_feats == 0:
        return

    cols = max(1, target_bytes // max(data.itemsize * n_samps, 1))
    cols = min(cols, n_feats)
    fitted = np.empty((n_samps, cols), dtype=data.dtype)

    for start in range(0, n_feats, cols):
        stop = min(start + cols, n_feats)
        block = fitted[:, : stop - start]
        np.matmul(dmat, beta[:, start:stop], out=block)
        data[:, start:stop] -= block


def _lstsq_fit(dmat, gram=False):
    """Calculate a least-squares operator and optional inverse-Gram sum vector.

    Both operators are derived from one singular value decomposition.
    """
    U, S, Vh = np.linalg.svd(dmat, full_matrices=False)
    S_inv = _invert_singular(S)
    V = Vh.T
    dmat_inv = (V * S_inv) @ U.T
    if gram:
        # Calculate a sum vector of inverse Gram matrix
        gram_sum = V @ ((S_inv**2) * np.sum(Vh, axis=1))
        # If one needs the full inverse Gram matrix, it is:
        # gmat_inv = (V * S_inv**2) @ Vh
    else:
        gram_sum = None
    return dmat_inv, gram_sum


def _invert_singular(S):
    """Invert significant singular values.

    The Gram matrix may be singular, if there are colinear covariates in the design
    matrix. In such cases, a direct `inv` will raise an error, and `pinv` (Moore-Penrose
    pseudoinverse)is the robust choice.

    Small singular values are set to zero prior to inversion. The process and threshold
    (1e-15) are consistent with the underlying algorithm of `pinv`.
    """
    cutoff = 1e-15 * np.max(S, axis=-1, keepdims=True)
    return np.divide(1.0, S, out=np.zeros_like(S), where=S > cutoff)


def _calc_grouped_var_cov(res2, hmat, groups):
    """Calculate variances and covariance matrix of coefficients.

    Parameters
    ----------
    res2 : ndarray of shape (n_samples, n_features)
        Squared, centered residuals. Any float type. Missing values must already be
        replaced with zero.
    hmat : ndarray of shape (n_samples, n_covariates)
        Influence matrix H. The same float type as `res2`.
    groups : ndarray of shape (n_groups - 1,), bool, or None
        Indices of coefficients representing groups.

    Returns
    -------
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances of regression coefficients.
    covmat : ndarray or None
        Estimated covariance matrices of coefficients.

    """
    if groups is True:
        covmat = _calc_covariance(res2, hmat)
        var_hat = np.diagonal(covmat, axis1=1, axis2=2)
        # Make a writable F-contiguous copy of the read-only diagonal view.
        var_hat = np.asfortranarray(var_hat)
    elif groups is not None:
        var_hat, covmat = _calc_var_cov(res2, hmat, groups)
    else:
        var_hat = _calc_variance(res2, hmat)
        covmat = None
    return var_hat, covmat


def _calc_variance(res2, hmat):
    """Calculate coefficient variances without forming full covariance matrices.

    For feature i, there is:

        diag(Cov(beta_i)) = sum_j resid_sq[j, i] * influence[j]**2

    This is a single matrix multiplication with an (n_features, n_covariates) output.

    This function outputs an F-contiguous array. Otherwise, it can be replaced with:

        `return resid_sq.T @ influence**2`
    """
    var_hat = np.empty((res2.shape[1], hmat.shape[1]), dtype=res2.dtype, order="F")
    np.matmul(res2.T, np.square(hmat), out=var_hat)
    return var_hat


def _calc_covariance(res2, hmat):
    """Calculate full covariance matrices.

    For feature i, the covariance is:

        Cov(beta_i) = sum_j resid2_ji * h_j h_j.T

    where H = X @ (X.T @ X)^+ is the influence matrix (^+ means Moore-Penrose
    pseudoinverse). This H is given by `dmat_inv.T` calculated upstream.

    Flattening h_j h_j.T turns the contraction over samples into one matrix
    multiplication across all features.
    """
    n_samps, n_covars = hmat.shape
    prods = (hmat[:, :, None] * hmat[:, None, :]).reshape(n_samps, n_covars * n_covars)
    return (res2.T @ prods).reshape(res2.shape[1], n_covars, n_covars)


def _calc_var_cov(res2, hmat, groups):
    """Calculate all variances and one symmetric covariance submatrix.

    All coefficient variances are calculated through `_calc_variance`, exactly the same
    path used when no grouping is requested. This guarantees that enabling post-hoc
    analysis does not perturb the primary result.

    A second, smaller matrix multiplication calculates only the k * (k - 1) / 2 unique
    off-diagonal covariances needed by the selected grouping, avoiding calculating the
    full covariance matrix.
    """
    n_samps, n_covars = hmat.shape
    n_feats = res2.shape[1]
    k = groups.size

    # Keep the primary variances bit-for-bit on the same calculation path regardless
    # of whether a grouping was requested.
    var_hat = _calc_variance(res2, hmat)
    covmat = np.empty((n_feats, k, k), dtype=var_hat.dtype)
    diag_idx = np.arange(k)
    covmat[:, diag_idx, diag_idx] = var_hat[:, groups]

    tri_i, tri_j = np.triu_indices(k, 1)
    n_offdiag = tri_i.size
    if n_offdiag:
        dtype = np.result_type(res2, hmat)
        pair_products = np.empty((n_samps, n_offdiag), dtype=dtype)
        np.multiply(
            hmat[:, groups[tri_i]],
            hmat[:, groups[tri_j]],
            out=pair_products,
        )
        offdiag = np.empty((n_feats, n_offdiag), dtype=dtype)
        np.matmul(res2.T, pair_products, out=offdiag)
        covmat[:, tri_i, tri_j] = offdiag
        covmat[:, tri_j, tri_i] = offdiag

    return var_hat, covmat


# Numerical optimization setting for sampling bias estimation through an Expectation-
# maximization (EM) algorithm.
#
# The Nelder-Mead simplex algorithm is used, to be consistent with the R code. Note: It
# does not actually enforce bounds during optimization. Rather, it merely clips at the
# bounds. Noted for further consideration.
#
# Tighter convergence tolerance produce more precise estimates in each EM step, leading
# to smaller accumulated errors. SciPy's default is 1e-4, which would produce different
# results from the R code. The value 1e-8 was selected according to our benchmarks.
OPTARG = dict(
    method="Nelder-Mead",
    bounds=[(0, None)],
    options={"xatol": 1e-8, "fatol": 1e-8},
)


def _estimate_bias_em(beta, var_hat, tol=1e-5, max_iter=100):
    """Estimate sampling bias through an expectation-maximization (EM) algorithm.

    This function models the observed coefficients (log-fold changes) for a given
    covariate as a Gaussian mixture distribution with three components: 0) null,
    1) negative, 2) positive. It aims to estimate a global bias term (delta) that
    affects all features.

    Parameters
    ----------
    beta : ndarray of shape (n_features,)
        Estimated coefficients (log-fold change before correction).
    var_hat : ndarray of shape (n_features,)
        Estimated variances of regression coefficients.
    tol : float, optional
        Absolute tolerance of EM iteration. Default is 1e-5.
    max_iter : int
        Max number of iteration of EM iteration. Default is 100.

    Returns
    -------
    delta_em : float
        EM estimator of bias.
    delta_wls : float
        WLS estimator of bias.
    var_delta : float
        Estimated variances of bias.

    """
    # This function involves careful memory optimization to avoid creating temporary
    # arrays during iteration. Technically, the arrays could have been pre-allocated
    # in the outer function `_ancombc_core` and re-used across covariates, which could
    # further enhance memory efficiency. It is left as-is for modularity.

    # The original R code has `na.rm = TRUE` in many commands. This is not necessary
    # in the current implementation, because the pre-correction coefficients (beta)
    # is guaranteed to not contain NaN values.

    # Mask NaN values (deemed unnecessary; left here for future examination).
    # beta = beta[~np.isnan(beta)]

    # There might be a chance that numeric optimization produces (near) zero weights
    # (pi) or variances (kappa), which can cause numerical stability issues in the
    # EM process. To safe-guard, one may use a small number `eps` as the floor of
    # those parameters. The original R code doesn't have this mechanism. Therefore, it
    # is currently disabled.
    # eps = 1e-12

    # Initial model parameters
    pi0, pi1, pi2 = 0.75, 0.125, 0.125  # weights of components (pi)
    delta, l1, l2, kappa1, kappa2 = _init_bias_params(beta)
    params = np.array([pi0, pi1, pi2, delta, l1, l2, kappa1, kappa2])
    updated = np.empty(8)

    # Pre-allocate memory for intermediates. Each array has three rows, representing
    # the three components (0, 1, 2), and columns representing individual features.
    # Let delta_i = mu_i1_hat - mu_i2_hat
    # The distribution of delta_i is modeled by Gaussian mixture:
    # f(delta_i) = pi0 * phi((delta_i - delta) / nu_i0) +
    #              pi1 * phi((delta_i - (delta + l1)) / nu_i1) +
    #              pi2 * phi((delta_i - (delta + l2)) / nu_i2)
    # where phi is the normal density function,
    # (delta + l1) and (delta + l2) are means for delta_i | C1 and delta_i | C2
    # nu_i0, nu_i1, and nu_i2 are variances of delta_i | C0, delta_i | C1, and
    # delta_i | C2 respectively.
    # We assume nu_i1 = nu_i0 + kappa1 and nu_i2 = nu_i0 + kappa2 for computational
    # simplicity.
    n_feats = beta.shape[0]
    shape = (3, n_feats)
    nu_inv = np.empty(shape)  # inverse of variances
    inv_stdevs = np.empty(shape)  # inverse standard deviations
    ratios = np.empty(shape)  # coefficients / variances

    # Mean coefficients
    means = np.empty(3)

    # Posterior probabilities of feature-component assignments (EM's responsibilities)
    resp = np.empty(shape)

    # Just a 2-row array to store random data
    intm = np.empty((2, n_feats))
    intm0, intm1 = intm

    # Initialize intermediates. The 1st row is constant, representing pre-correction
    # estimates, whereas the 2nd and 3rd rows are to be modified during iteration.
    # NOTE: Making `var_hat` C-contiguous can further improve performance, but this is
    # marginal compared with the EM process.
    np.reciprocal(var_hat, out=nu_inv[0])
    np.sqrt(nu_inv[0], out=inv_stdevs[0])
    np.divide(beta, var_hat, out=ratios[0])

    # Objective function for optimization of variance estimation. For fixed `loc` and
    # responsibilities, the negative weighted Gaussian log-likelihood is, up to an
    # additive constant:
    #
    #   0.5 * sum(resp * ((beta - loc)**2 / (var_hat + x)
    #                     + log(var_hat + x)))
    #
    # The omitted 0.5 * log(2 * pi) * sum(resp) term is independent of `x` and
    # therefore cannot affect the minimizer. The squared residual `(beta - loc)**2`
    # is also independent of `x`; it is computed once per component below and passed
    # into this callback, instead of being recomputed at every Nelder-Mead evaluation.
    # The factor 0.5 is retained because `OPTARG` uses an absolute function tolerance
    # (`fatol`); dropping the factor would rescale the stopping criterion.
    def func(x, sq_diff, resp):
        np.add(var_hat, x[0], out=intm0)
        np.log(intm0, out=intm1)
        np.divide(sq_diff, intm0, out=intm0)
        np.add(intm0, intm1, out=intm0)
        return 0.5 * np.dot(resp, intm0)

    # Expectation-maximization (E-M) iterations
    loss, epoch = np.inf, 0
    while loss > tol and epoch < max_iter:
        # Update intermediates (2nd and 3rd rows only)
        np.add(var_hat, params[6:8, None], out=intm)  # kappa1, kappa2
        np.reciprocal(intm, out=nu_inv[1:])
        np.sqrt(nu_inv[1:], out=inv_stdevs[1:])
        np.subtract(beta, params[4:6, None], out=ratios[1:])  # means (l)
        ratios[1:] *= nu_inv[1:]

        ### E-step ###
        # Mean coefficients
        delta = means[0] = params[3]  # global bias (delta)
        np.add(params[4:6], delta, means[1:])

        # Posterior probabilities = Gaussian densities weighted by component fractions.
        # The normalizing constant 1 / sqrt(2 * pi) is shared by all components and
        # cancels when responsibilities are normalized, so it is omitted. Compute the
        # remaining density directly in pre-allocated memory:
        #
        #   exp(-0.5 * (beta - mean_r)**2 / nu_r) / sqrt(nu_r)
        #
        # where `nu_inv` is 1 / nu_r and `inv_stdevs` is 1 / sqrt(nu_r).
        np.subtract(beta, means[:, None], out=resp)
        np.square(resp, out=resp)
        resp *= nu_inv
        resp *= -0.5
        np.exp(resp, out=resp)
        resp *= inv_stdevs
        resp *= params[:3, None]  # weights (pi)
        resp /= np.sum(resp, axis=0, keepdims=True)

        ### M-step ###
        # Weights of components (pi)
        # pi_r_new = mean(pi_r * pdf_r / (pi0 * pdf0 + pi1 * pdf1 + pi2 * pdf2)),
        # where r = 0, 1, 2
        np.mean(resp, axis=1, out=updated[:3])

        # Avoid zero weights.
        # updated[:3] = np.maximum(updated[:3], eps)
        # updated[:3] /= updated[:3].sum()

        # Gaussian mixture modeling of global bias (delta)
        # The following code produces the same result as:
        #   updated[3] = np.sum(resp * ratios) / np.sum(resp * nu_inv)
        # But it avoids creating intermediate arrays.
        # delta_new = sum(r_0i * beta / nu0 +
        #                 r_1i * (beta - l1) / (nu0 + kappa1) +
        #                 r_2i * (beta - l2) / (nu0 + kappa2)) /
        #             sum(r0i / nu0 + r1i / (nu0 + kappa1) + r2i / (nu0 + kappa2))
        updated[3] = np.vdot(resp, ratios) / np.vdot(resp, nu_inv)

        # Negative and positive components relative to delta (l)
        np.multiply(resp[1:], nu_inv[1:], out=intm)
        denom = np.sum(intm, axis=1)
        np.subtract(beta, delta, out=resp[0])  # reuse as intermediate
        intm *= resp[0]
        numer = np.sum(intm, axis=1)
        l1, l2 = numer / denom
        # l1_new = min(sum(r1i * (beta - delta) / (nu0 + kappa1)) /
        #              sum(r1i / (nu0 + kappa1)), 0)
        updated[4] = np.minimum(l1, 0)
        # l2_new = min(sum(r2i * (beta - delta) / (nu0 + kappa2)) /
        #              sum(r2i / (nu0 + kappa2)), 0)
        updated[5] = np.maximum(l2, 0)

        # Perform numeric optimization to minimize variances of negative and positive
        # components (kappa). The component means and responsibilities stay fixed
        # throughout each M-step optimization, so precompute squared residuals once.
        # `ratios[1:]` is no longer needed after the delta update above and is reused
        # here as workspace.
        np.subtract(beta, means[1:, None], out=ratios[1:])
        np.square(ratios[1:], out=ratios[1:])

        # TODO: Consider scenarios where optimization doesn't converge.
        updated[6] = minimize(func, params[6], args=(ratios[1], resp[1]), **OPTARG).x[0]
        updated[7] = minimize(func, params[7], args=(ratios[2], resp[2]), **OPTARG).x[0]

        # Avoid zero variances.
        # updated[6] = max(updated[6], eps)
        # updated[7] = max(updated[7], eps)

        # Loss (epsilon)
        # epsilon = sqrt((pi0_new - pi0)^2 + (pi1_new - pi1)^2 + (pi2_new - pi2)^2 +
        #                (delta_new - delta)^2 + (l1_new - l1)^2 + (l2_new - l2)^2 +
        #                (kappa1_new - kappa1)^2 + (kappa2_new - kappa2)^2)
        loss = np.linalg.norm(updated - params)

        params[:] = updated
        epoch += 1

    return _estimate_bias_var(beta, var_hat, params)


def _init_bias_params(beta):
    """Initialize parameters for iterative bias estimation.

    Parameters
    ----------
    beta : ndarray of shape (n_features,)
        Estimated coefficients before correction.

    Returns
    -------
    delta : float
        Initial global bias term.
    l1, l2 : float
        Initial means of negative and positive components, respectively.
    kappa1, kappa2 : float
        Initial variances of negative and positive components, respectively.

    """
    edges = np.quantile(beta, (0.125, 0.25, 0.75, 0.875))

    # Estimate delta (mean of values between q1 and q3)
    if np.any(mask := (beta >= edges[1]) & (beta <= edges[2])):
        delta = np.mean(beta[mask])
    else:
        delta = np.mean(beta)

    # Estimate l1
    if np.any(mask := beta < edges[0]):
        l1 = np.mean(beta_ := beta[mask])
        if beta_.size > 1:
            kappa1 = np.var(beta_, ddof=1, mean=l1) or 1.0
        else:
            kappa1 = 1.0
    else:
        l1 = np.min(beta)
        kappa1 = 1.0

    # Estimate l2
    if np.any(mask := beta > edges[3]):
        l2 = np.mean(beta_ := beta[mask])
        if beta_.size > 1:
            kappa2 = np.var(beta_, ddof=1, mean=l2) or 1.0
        else:
            kappa2 = 1.0
    else:
        l2 = np.max(beta)
        kappa2 = 1.0

    return delta, l1, l2, kappa1, kappa2


def _estimate_bias_var(beta, var_hat, params):
    """Estimate variance of sampling bias according to EM-optimized parameters.

    This function calculates the weighted least squares (WLS) estimator of bias.

    Parameters
    ----------
    beta : ndarray of shape (n_features,)
        Estimated coefficients (log-fold change before correction).
    var_hat : ndarray of shape (n_features,)
        Estimated variances of regression coefficients.
    params : ndarray of shape (8,)
        Model parameters optimized by EM.

    Returns
    -------
    delta_em : float
        EM estimator of bias.
    delta_wls : float
        WLS estimator of bias.
    var_delta : float
        Estimated variance of bias.

    """
    pi0, pi1, pi2, delta, l1, l2, kappa1, kappa2 = params

    # The EM estimator of bias
    delta_em = delta

    # Assign features to Gaussian components (C)
    q1, q2 = np.quantile(beta, (pi1, 1.0 - pi2))
    C1 = np.flatnonzero(beta < q1)
    C2 = np.flatnonzero(beta >= q2)

    # Numerator of the WLS estimator
    nu_inv = var_hat.copy()
    nu_inv[C1] += kappa1
    nu_inv[C2] += kappa2
    nu_inv[:] = 1.0 / nu_inv
    wls_denom = np.sum(nu_inv)

    # Denominator of the WLS estimator
    beta_ = beta.copy()
    beta_[C1] -= l1
    beta_[C2] -= l2
    nu_inv *= beta_
    wls_numer = np.sum(nu_inv)

    # Estimate the variance of bias
    wls_denom_inv = 1.0 / wls_denom
    delta_wls = wls_numer * wls_denom_inv
    var_delta = np.nan_to_num(wls_denom_inv)

    return delta_em, delta_wls, var_delta


def _sample_fractions(data, dmat, beta, delta_em, missing=None):
    """Estimate sampling fractions.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Data table. Zero-handled. Transformed.
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.
    beta : ndarray of shape (n_covariates, n_features)
        Estimated coefficients.
    delta_em : ndarray of shape (n_covariates,)
        Estimated biases.
    missing : ndarray of shape (n_samples, n_features), optional
        Boolean mask of unobserved entries. When absent, use an algebraically equivalent
        dense calculation that avoids an n_samples x n_features residual array.

    Returns
    -------
    theta_hat : ndarray of shape (n_covariates,)
        Estimated sampling fractions.

    """
    # Compute sampling fractions. For dense data, linearity of the mean gives
    #
    # mean(data - X @ beta + X @ delta, axis=features)
    #   = mean(data) - X @ mean(beta) + X @ delta,
    #
    # avoiding another n_samples x n_features matrix while `data` must remain intact.
    if missing is None:
        theta_hat = np.mean(data, axis=1)
        theta_hat -= dmat @ np.mean(beta, axis=1)
        theta_hat += dmat @ delta_em
    else:
        # With feature-specific missingness the mean cannot be pulled through X @ beta,
        # so retain the existing NaN-aware calculation.
        intm = dmat @ beta
        intm -= (dmat @ delta_em)[:, None]
        intm -= data
        intm *= -1.0
        theta_hat = np.nanmean(intm, axis=1)

    # Handle NaN in theta (samples with all-NaN residuals)
    # TODO: This may not be necessary if empty samples are not allowed.
    nan_theta = np.isnan(theta_hat)
    if np.any(nan_theta):
        theta_hat[nan_theta] = 0.0

    return theta_hat


def _adjust_variances(var_hat, vcov_hat, var_delta, s0_perc, groups=None):
    """Adjust variances.

    Parameters
    ----------
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances.
    vcov_hat : ndarray or None
        Full coefficient covariance matrices or a retained covariance submatrix. May be
        None when only coefficient variances were requested.
    var_delta : ndarray of shape (n_covariates,)
        Estimated variances of biases.
    s0_perc : float
        SAM-like fudge factor.
    covariance_indices : 1-D array_like of int, optional
        Full-model coefficient indices represented by ``vcov_hat`` when it is a
        covariance submatrix.

    Notes
    -----
    This function updates ``var_hat`` in place and, when provided, updates the diagonal
    of ``vcov_hat`` to remain consistent.

    """
    # vars += delta + 2 * |vars * delta|^0.5
    var_delta_t = var_delta[None, :]
    var_prod = var_hat * var_delta_t
    np.abs(var_prod, out=var_prod)
    np.sqrt(var_prod, out=var_prod)
    var_prod *= 2.0
    var_hat += var_delta_t
    var_hat += var_prod

    # SAM-like fudge factor
    # TODO: Handle NaN values
    if s0_perc:
        s02 = np.nanquantile(var_hat, s0_perc, axis=0)
        var_hat += s02[None, :]
        # var_hat[np.isnan(beta_hat)] = np.nan

    # Keep a retained covariance tensor consistent with the adjusted variances.
    if vcov_hat is not None:
        if groups is None:
            groups = np.arange(var_hat.shape[1])
        diag_idx = np.arange(groups.size)
        vcov_hat[:, diag_idx, diag_idx] = var_hat[:, groups]


def _format_results(beta_hat, var_hat, features, covariates, alpha, p_adjust, dof=None):
    """Format primary ANCOM-BC/BC2 statistics as a DataFrame.

    This function is specialized for the public result-construction path. Unlike
    :func:`_calc_statistics`, which returns every intermediate statistic as a NumPy
    array, this routine inserts each statistic into the DataFrame as soon as it is
    calculated and then reuses its temporary workspace. This reduces peak memory for
    large feature-by-covariate result tables.

    Parameters
    ----------
    beta_hat : ndarray of shape (n_features, n_covariates)
        Estimated coefficients post correction.
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances.
    features : 1-D array_like of length n_features
        Feature identifiers.
    covariates : 1-D array_like of length n_covariates
        Covariate identifiers.
    alpha : float
        Significance level.
    p_adjust : str
        FDR correction method.
    dof : float or ndarray of shape (n_features,), optional
        Degrees of freedom.

    Returns
    -------
    pd.DataFrame
        Primary results with a (FeatureID, Covariate) MultiIndex and columns
        Log2(FC), SE, W, pvalue, qvalue and Signif.

    """
    beta_hat = np.asarray(beta_hat)
    var_hat = np.asarray(var_hat)
    if beta_hat.shape != var_hat.shape:
        raise ValueError("`beta_hat` and `var_hat` must have matching shapes.")

    n_feats, n_covars = beta_hat.shape
    if len(features) != n_feats or len(covariates) != n_covars:
        raise ValueError(
            "Feature and covariate identifiers must match the result dimensions."
        )

    # Construct the MultiIndex directly. This avoids materializing two Python lists
    # of length n_features * n_covariates before converting them into an index.
    index = pd.MultiIndex.from_product(
        (features, covariates), names=("FeatureID", "Covariate")
    )
    res = pd.DataFrame(index=index)
    res["Log2(FC)"] = beta_hat.ravel()

    # A single feature-by-covariate workspace is enough for SE and W because pandas
    # copies each column on assignment.
    work = np.sqrt(var_hat)
    res["SE"] = work.ravel()
    np.divide(beta_hat, work, out=work)
    res["W"] = work.ravel()

    # scipy.stats allocates the p-value array. Once p-values have been copied into
    # the DataFrame, reuse that same array for adjusted p-values.
    np.abs(work, out=work)
    pval = _calc_pvalues_abs(work, dof)
    res["pvalue"] = pval.ravel()
    del work

    _adjust_pvalues(pval, p_adjust, out=pval)
    res["qvalue"] = pval.ravel()
    res["Signif"] = pd.array(pval.ravel() <= alpha, dtype="boolean")

    return res


def _calc_statistics(beta_hat, var_hat, alpha, p_adjust, dof=None):
    """Calculate statistical significance while correcting for multiple testing.

    Parameters
    ----------
    beta_hat : ndarray of shape (n_features, n_covariates)
        Estimated coefficients post correction.
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances.
    alpha : float
        Significance level.
    p_adjust : str
        FDR correction method.
    dof : float or ndarray of shape (n_features,), optional
        Degrees of freedom.

    Returns
    -------
    se_hat : ndarray of shape (n_features, n_covariates)
        Estimated standard errors.
    W : ndarray of shape (n_features, n_covariates)
        Test statistics.
    pval : ndarray of shape (n_features, n_covariates)
        p-values.
    qval : ndarray of shape (n_features, n_covariates)
        Adjusted p-values.
    reject : ndarray of shape (n_features, n_covariates)
        Significant differential abundance (reject null hypothesis).

    """
    # Estimate standard error
    se_hat = np.sqrt(var_hat)
    # NOTE: Code was below. But I think this is unnecessary: First, var_hat cannot
    # contain negative values. Second, np.sqrt on NaN values will return NaN without
    # warning.
    # se_hat = np.sqrt(np.maximum(var_hat, 0))
    # se_hat[np.isnan(var_hat)] = np.nan

    # Calculate test statistic (W)
    W = beta_hat / se_hat
    # NOTE: Original code is below. But I think there is no need to handle NaN here.
    # W = np.where(np.isnan(se_hat), np.nan, beta_hat / se_hat)

    # Calculate p-values
    pval = _calc_pvalues(W, dof)

    # FDR correction of p-values
    qval = _adjust_pvalues(pval, p_adjust)

    # Reject null hypothesis
    reject = qval <= alpha

    return se_hat, W, pval, qval, reject


def _calc_pvalues(W, dof=None):
    """Calculate p-values from test statistics (W).

    p-values are calculated using t-test (with degrees of freedom (dof) provided, as in
    ANCOM-BC2) or Z-test (without dof, as in ANCOM-BC).

    Parameters
    ----------
    W : ndarray of shape (n_features, n_covariates)
        Test statistics.
    dof : float or ndarray of shape (n_features,), optional
        Degrees of freedom.

    Returns
    -------
    pval : ndarray of shape (n_features, n_covariates)
        p-values.

    """
    return _calc_pvalues_abs(np.abs(W), dof)


def _calc_pvalues_abs(W_abs, dof=None):
    """Calculate two-sided p-values from absolute test statistics."""
    if dof is not None:  # t-test with dof
        if not np.isscalar(dof):
            dof = np.asarray(dof)[:, None]  # broadcast to 2D
        pval = t.sf(W_abs, dof)
    else:  # Z-test
        pval = norm.sf(W_abs)
    pval *= 2.0  # two-sided test
    return pval


def _adjust_pvalues(pval, method, out=None):
    """Adjust p-values for multiple-testing correction.

    This function applies FDR correction to non-NaN entries. This behavior matches R's
    `p.adjust`, whereas statsmodels' `multipletests` has inconsistent behavior in some
    methods.

    NaN p-values could emerge when there are NaN in dof, which in turn could happen if
    there are less observed samples than covariates in some features. This issue is
    independent from the zero values in the input data.

    Parameters
    ----------
    pval : ndarray of shape (n_features, n_covariates)
        p-values.
    out : ndarray, optional
        Array in which to store adjusted p-values. May be the same array as ``pval``.
        By default, allocate a new array.

    Returns
    -------
    qval : ndarray of shape (n_features, n_covariates)
        q-values.

    """
    pval = np.asarray(pval)
    valid = ~np.isnan(pval)
    if out is None:
        qval = np.full_like(pval, np.nan)
    else:
        qval = np.asarray(out)
        if qval.shape != pval.shape:
            raise ValueError("`out` must have the same shape as `pval`.")
        if not np.can_cast(pval.dtype, qval.dtype, casting="same_kind"):
            raise TypeError("`out` has an incompatible dtype.")
        # If output doesn't alias the input, initialize invalid entries to NaN. When
        # it does alias, those entries are already NaN and clearing the array first
        # would destroy the p-values before adjustment.
        if not np.shares_memory(qval, pval):
            qval.fill(np.nan)
    # Holm and Benjamini-Hochberg are the documented/common choices. Their formulas
    # are simple enough to apply directly, avoiding repeated statsmodels dispatch and
    # validation overhead for every covariate. Processing one column at a time keeps
    # scratch space O(n_features), rather than allocating another full result matrix.
    key = None if method is None else str(method).lower()
    holm = key in {"holm", "holm-bonferroni"}
    bh = key in {"bh", "fdr_bh", "benjamini-hochberg"}

    if key is None:
        if not np.shares_memory(qval, pval):
            np.copyto(qval, pval)
        return qval

    if holm or bh:
        cols = (None,) if pval.ndim == 1 else range(pval.shape[1])
        for col in cols:
            pcol = pval if col is None else pval[:, col]
            qcol = qval if col is None else qval[:, col]
            valid_idx = np.flatnonzero(~np.isnan(pcol))
            n = valid_idx.size
            if not n:
                continue

            values = pcol[valid_idx]
            order = np.argsort(values)
            adjusted = values[order].copy()
            if holm:
                adjusted *= n - np.arange(n)
                np.maximum.accumulate(adjusted, out=adjusted)
            else:
                adjusted /= np.arange(1, n + 1) / float(n)
                adjusted[:] = np.minimum.accumulate(adjusted[::-1])[::-1]
            np.minimum(adjusted, 1.0, out=adjusted)
            qcol[valid_idx[order]] = adjusted
        return qval

    # Preserve support for every other method accepted by `_check_p_adjust`.
    func = _check_p_adjust(method)
    if pval.ndim == 1:
        qval[valid] = func(pval[valid])
    else:
        for col in range(pval.shape[1]):
            valid_col = valid[:, col]
            qval[valid_col, col] = func(pval[valid_col, col])

    return qval


class ANCOMBCResult:
    """Results for ANCOM-BC and ANCOM-BC2 analyses.

    This class contains the primary differential abundance results. Post-hoc analyses
    are available as methods that compute on-demand when a ``grouping`` was specified
    in the upstream :func:`ancombc` or :func:`ancombc2` call. Only the covariance
    submatrices required for that grouping are retained.

    Attributes
    ----------
    res : pd.DataFrame
        Primary results with (FeatureID, Covariate) multi-index. Columns are: Log2(FC),
        SE, W, pvalue, qvalue, Signif.
    method : {"ANCOM-BC", "ANCOM-BC2"}
        Differential abundance method used for the analysis.
    grouping : str or None
        Grouping term selected for post-hoc analysis, or None when post-hoc analyses
        were not enabled.
    has_covariance : bool
        Whether a grouping covariance submatrix was retained.

    Methods
    -------
    global_test
        Global test for differential abundance across >= 3 groups.
    dunnett_test
        Dunnett's test: each group vs. reference, with mdFDR correction.
    pairwise_test
        Pairwise directional test between all group pairs, with mdFDR.
    trend_test
        Trend test for ordered patterns in group effects.
    sensitivity_analysis
        Pseudo-count sensitivity analysis for robustness assessment.

    See Also
    --------
    ancombc : ANCOM-BC function.
    ancombc2 : ANCOM-BC2 function.
    """

    _private_defaults = {
        "_dmat": None,
        "_beta_hat": None,
        "_var_hat": None,
        "_vcov_hat": None,
        "_grouping": None,
        "_group_indices": None,
        "_dof": None,
        "_features": None,
        "_covariates": None,
        "_alpha": 0.05,
        "_p_adjust": "holm",
        "_s0_perc": 0.05,
        "_max_iter": 100,
        "_tol": 1e-5,
        "_pseudo": 0,
    }

    def __init__(
        self,
        res: pd.DataFrame,
        method: str,
        **kwargs,
    ):
        unexpected = set(kwargs).difference(self._private_defaults)
        if unexpected:
            names = ", ".join(sorted(unexpected))
            raise TypeError(f"Unexpected ANCOMBCResult argument(s): {names}")
        if method not in {"ANCOM-BC", "ANCOM-BC2"}:
            raise ValueError("`method` must be either 'ANCOM-BC' or 'ANCOM-BC2'.")

        self.res = res
        self._method = method
        for name, default in self._private_defaults.items():
            setattr(self, name, kwargs.get(name, default))

    @property
    def res(self) -> pd.DataFrame:
        return self._res

    @res.setter
    def res(self, value: pd.DataFrame):
        self._res = value

    @property
    def method(self) -> str:
        return self._method

    @property
    def grouping(self) -> str | None:
        """Grouping term selected for post-hoc analysis, if any."""
        return self._grouping

    @property
    def has_covariance(self) -> bool:
        """Whether a grouping covariance submatrix was retained."""
        return self._vcov_hat is not None

    def _require_grouping(self, method):
        """Require an upstream grouping for post-hoc analyses."""
        if self._grouping is None:
            func = "ancombc2" if self.method == "ANCOM-BC2" else "ancombc"
            raise ValueError(
                f"`{method}` requires a post-hoc grouping. Rerun "
                f"`{func}(..., grouping=<metadata column>)` to enable post-hoc "
                "analysis."
            )

    def __getitem__(self, key):
        return getattr(self, key)

    def keys(self):
        """Return a list of attribute names that are not private and are not None."""
        return [
            name
            for name in (
                "res",
                "method",
            )
            if getattr(self, name) is not None
        ]

    # formated output
    def __repr__(self):
        n_feats = len(self.res.index.get_level_values("FeatureID").unique())
        n_covars = len(self.res.index.get_level_values("Covariate").unique())
        n_signif = int(self.res["Signif"].sum())
        return (
            f"ANCOMBCResult(method={self.method!r}, "
            f"n_taxa={n_feats}, n_covariates={n_covars}, "
            f"n_signif={n_signif})"
        )

    def _stat_params(self, alpha, p_adjust):
        """Get FDR-correction method and significance level from upstream."""
        if alpha == "inherit":
            alpha = self._alpha
        if p_adjust == "inherit":
            p_adjust = self._p_adjust
        return alpha, p_adjust

    def global_test(
        self, alpha: float | str = "inherit", p_adjust: str = "inherit"
    ) -> pd.DataFrame:
        """Perform global test for differential abundance across groups.

        The global test identifies features that are differentially abundant between at
        least two groups across three or more groups.

        Parameters
        ----------
        alpha : float or "inherit", optional
            Significance level, or the value used by :func:`ancombc` or
            :func:`ancombc2`. Default is "inherit".
        p_adjust : str, optional
            Multiple-testing correction method, or "inherit" to use the value
            supplied upstream. Default is "inherit".

        Returns
        -------
        res_global : pd.DataFrame
            Global test result. Columns are:

            - ``FeatureID``: Feature identifier, i.e., dependent variable.

            - ``W``: *W*-statistic quantifying the overall evidence against null
              hypothesis (mean abundance of the feature is the same across all groups).

            - ``pvalue``: *p*-value of the *W*-statistic.

            - ``qvalue``: Corrected *p*-value.

            - ``Signif``: Whether at least one group mean is different from others.

        """
        self._require_grouping("global_test")
        alpha, p_adjust = self._stat_params(alpha, p_adjust)
        W, pval, qval, reject = _global_test(
            dmat=self._dmat,
            grouping=self._grouping,
            beta_hat=self._beta_hat,
            vcov_hat=self._vcov_hat,
            group_indices=self._group_indices,
            p_adjust=p_adjust,
            alpha=alpha,
            dof=self._dof,
        )
        result = pd.DataFrame(
            {"W": W, "pvalue": pval, "qvalue": qval, "Signif": reject},
            index=self._features,
        )
        result.index.name = "FeatureID"
        return result

    def pairwise_test(
        self,
        alpha: float | str = "inherit",
        p_adjust: str = "inherit",
    ) -> pd.DataFrame:
        """Perform pairwise directional test between all group pairs.

        Uses mixed directional FDR (mdFDR) correction via bootstrap.

        Parameters
        ----------
        alpha : float or "inherit", optional
            Significance level, or the value supplied upstream. Default is
            "inherit".
        p_adjust : str, optional
            Family wise error (FWER) controlling method. Default is "inhert", which
            will use the *p*-value correction method supplied upstream.

        Returns
        -------
        pd.DataFrame
            DataFrame with (FeatureID, Comparison) multi-index and columns:
            Log2(FC), SE, W, pvalue, qvalue, Signif.

        """
        self._require_grouping("pairwise_test")
        alpha, p_adjust = self._stat_params(alpha, p_adjust)

        raw = _pairwise_test(
            dmat=self._dmat,
            grouping=self._grouping,
            beta_hat=self._beta_hat,
            var_hat=self._var_hat,
            vcov_hat=self._vcov_hat,
            group_indices=self._group_indices,
            dof=self._dof,
            p_adjust=p_adjust,
            alpha=alpha,
        )
        comp_names = raw["comp_names"]
        n_comp = len(comp_names)
        n_feats = len(self._features)
        result = pd.DataFrame(
            {
                "FeatureID": [x for x in self._features for _ in range(n_comp)],
                "Comparison": comp_names * n_feats,
                "Log2(FC)": raw["beta"].ravel(),
                "SE": raw["se"].ravel(),
                "W": raw["W"].ravel(),
                "pvalue": raw["p_val"].ravel(),
                "qvalue": raw["q_val"].ravel(),
                "Signif": raw["reject"].ravel(),
            }
        )
        result.set_index(["FeatureID", "Comparison"], inplace=True)
        return result

    def dunnett_test(
        self,
        alpha: float | str = "inherit",
        p_adjust: str = "inherit",
        bootstraps: int = 100,
        seed: SeedLike | None = None,
    ) -> pd.DataFrame:
        """Perform Dunnett's test (each group vs. reference) with mdFDR.

        Parameters
        ----------
        alpha : float or "inherit", optional
            Significance level, or the value supplied upstream. Default is
            "inherit".
        p_adjust : str, optional
            Family wise error (FWER) controlling method. Default in "inhert", which
            will use the *p*-value correction method supplied upstream.
        bootstraps : int, optional
            Number of bootstrap iterations. Default is 100.
        seed : int, Generator, or RandomState, optional
            A user-provided random seed or generator for bootstrap samples. See
            :func:`details <skbio.util.get_rng>`.

        Returns
        -------
        pd.DataFrame
            DataFrame with (FeatureID, Comparison) multi-index and columns:
            Log2(FC), SE, W, pvalue, qvalue, Signif.

        """
        self._require_grouping("dunnett_test")
        rng = get_rng(seed)
        alpha, p_adjust = self._stat_params(alpha, p_adjust)

        raw = _dunnett_test(
            dmat=self._dmat,
            grouping=self._grouping,
            beta_hat=self._beta_hat,
            group_indices=self._group_indices,
            var_hat=self._var_hat,
            dof=self._dof,
            bootstraps=bootstraps,
            rng=rng,
            p_adjust=p_adjust,
            alpha=alpha,
        )
        comp_names = raw["comp_names"]
        n_comp = len(comp_names)
        n_feats = len(self._features)
        result = pd.DataFrame(
            {
                "FeatureID": [x for x in self._features for _ in range(n_comp)],
                "Comparison": comp_names * n_feats,
                "Log2(FC)": raw["beta"].ravel(),
                "SE": raw["se"].ravel(),
                "W": raw["W"].ravel(),
                "pvalue": raw["p_val"].ravel(),
                "qvalue": raw["q_val"].ravel(),
                "Signif": raw["reject"].ravel(),
            }
        )
        result.set_index(["FeatureID", "Comparison"], inplace=True)
        return result

    def trend_test(
        self,
        alpha: float | str = "inherit",
        p_adjust: str = "inherit",
        trend_contrast: dict | None = None,
        trend_node: dict | None = None,
        bootstraps: int = 100,
        seed: SeedLike | None = None,
    ) -> pd.DataFrame:
        r"""Perform trend test for ordered patterns in group effects.

        Uses constrained optimization to test monotone increasing/decreasing
        patterns in group-level effects.

        Parameters
        ----------
        alpha : float or "inherit", optional
            Significance level, or the value supplied upstream. Default is
            "inherit".
        p_adjust : str, optional
            Multiple-testing correction method, or "inherit" to use the value
            supplied upstream. Default is "inherit".
        trend_contrast, trend_node : dict, optional
            Trend-test contrast matrices and their node indices.
        bootstraps : int, optional
            Number of bootstrap iterations. Default is 100.
        seed : int, Generator, or RandomState, optional
            A user-provided random seed or generator for bootstrap samples. See
            :func:`details <skbio.util.get_rng>`.

        Returns
        -------
        pd.DataFrame
            DataFrame indexed by FeatureID with columns: W, pvalue, qvalue, Signif.

        Notes
        -----
        The trend test is highly stochastic and requires a large number of bootstraps
        (e.g., 10,000) to stabilize the estimated *p*- and *q*-values. This is because
        every *p*-value is estimated through independent bootstrapping. In comparsion,
        the Dunnett's test (:meth:`dunnett_test`) is more stable.

        """
        self._require_grouping("trend_test")
        rng = get_rng(seed)
        alpha, p_adjust = self._stat_params(alpha, p_adjust)

        raw = _trend_test(
            dmat=self._dmat,
            grouping=self._grouping,
            beta_hat=self._beta_hat,
            var_hat=self._var_hat,
            vcov_hat=self._vcov_hat,
            group_indices=self._group_indices,
            p_adjust=p_adjust,
            alpha=alpha,
            trend_contrast=trend_contrast,
            trend_node=trend_node,
            bootstraps=bootstraps,
            rng=rng,
        )
        result = pd.DataFrame(
            {
                "W": raw["W"],
                "pvalue": raw["p_val"],
                "qvalue": raw["q_val"],
                "Signif": raw["reject"],
            },
            index=self._features,
        )
        result.index.name = "FeatureID"
        return result


def _group_indices_from_design(dmat, grouping):
    """Return the design-matrix column indices for a named model term."""
    s = dmat.design_info.term_name_slices[grouping]
    return np.arange(s.start, s.stop, dtype=int)


def _select_group_covariance(vcov_hat, group_indices):
    """Return a grouping covariance submatrix from full or already-subset storage."""
    group_indices = np.asarray(group_indices, dtype=int)
    k = group_indices.size
    if vcov_hat.shape[1:] == (k, k):
        return vcov_hat
    return vcov_hat[:, group_indices][:, :, group_indices]


def _global_test(
    dmat,
    grouping,
    beta_hat,
    vcov_hat,
    alpha=0.05,
    p_adjust="holm",
    dof=None,
    group_indices=None,
):
    """Perform ANCOM-BC global test.

    The global test is to determine features that are differentially abundant between
    at least 2 sample groups across 3 or more groups.

    Parameters
    ----------
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.
    grouping : str
        The group variable of interests in metadata.
    beta_hat : ndarray of shape (n_features, n_covariates)
        Corrected coefficients.
    vcov_hat : ndarray
        Full coefficient covariance matrices or the covariance submatrices for the
        grouping term.
    alpha : float, optional
        Significance level for the statistical tests. Must be in the range of (0, 1).
        Default is 0.05.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are Holm-
        Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-
        Hochberg ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported
        by statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.
    dof : float or ndarray of shape (n_features,), optional
        Degrees of freedom. When provided, calculate p-values using F distribution (as
        in ANCOM-BC2), otherwise use chi-square distribution (as in ANCOM-BC).

    Returns
    -------
    W_global : ndarray of shape (n_features,)
        W-statistics of global test.
    pval : ndarray of shape (n_features,)
        p-values of global test.
    qval : ndarray of shape (n_features,)
        Adjusted p-values of global test.
    reject : ndarray of shape (n_features,)
        If the variable is differentially abundant.

    """
    if group_indices is None:
        group_ind = _group_indices_from_design(dmat, grouping)
    else:
        group_ind = np.asarray(group_indices, dtype=int)
    n_groups = group_ind.size

    # `vcov_hat` may be the full model covariance (legacy/internal calls) or the
    # grouping-only covariance retained by the public API.
    beta_hat_sub = beta_hat[:, group_ind]
    vcov_hat_sub = _select_group_covariance(vcov_hat, group_ind)

    # Inverse the subset of vcov_hat
    vcov_hat_sub_inv = np.linalg.pinv(vcov_hat_sub)

    # NOTE: The R code uses an identity matrix A, which is omitted here since it does
    # not change anything in multiplication. The following math is more efficient.
    W_global = np.einsum(
        "ni,nij,nj->n", beta_hat_sub, vcov_hat_sub_inv, beta_hat_sub, optimize=True
    )
    # NOTE: A more performant math is as follows. But it may be unsafe because `solve`
    # requires a non-singular square matrix.
    # intm = np.linalg.solve(vcov_hat_sub, beta_hat_sub[..., None])[..., 0]
    # W_global = np.einsum("ni,ni->n", beta_hat_sub, intm)

    # Calculate p-values. ANCOM-BC uses chi-square; ANCOM-BC2 uses F with per-feature
    # residual degrees of freedom.
    if dof is None:
        p_lower = chi2.cdf(W_global, n_groups)
        p_upper = chi2.sf(W_global, n_groups)
    else:
        p_lower = f.cdf(W_global, n_groups, dof)
        p_upper = f.sf(W_global, n_groups, dof)
    pval = 2 * np.minimum(p_lower, p_upper)

    # R's p.adjust excludes NA values; ANCOM-BC2 then treats invalid global
    # tests as nonsignificant.
    qval = _adjust_pvalues(pval, p_adjust)
    qval = np.where(np.isnan(qval), 1.0, qval)
    reject = qval <= alpha

    return W_global, pval, qval, reject


def _pairwise_test(
    dmat,
    grouping,
    beta_hat,
    var_hat,
    vcov_hat,
    dof,
    p_adjust="holm",
    alpha=0.05,
    group_indices=None,
):
    """ANCOM-BC2 pairwise directional test.

    For each pair of group levels, compute the difference in coefficients
    and its variance, then apply mdFDR correction.
    """
    covariates = dmat.design_info.column_names
    if group_indices is None:
        group_ind = _group_indices_from_design(dmat, grouping)
    else:
        group_ind = np.asarray(group_indices, dtype=int)

    beta_hat_sub = beta_hat[:, group_ind]
    group_covars = [covariates[i] for i in group_ind]
    vcov_group = _select_group_covariance(vcov_hat, group_ind)

    # Compute pairwise differences and their variances
    n_tax = beta_hat.shape[0]
    n_group = group_ind.size

    # Generate all pairwise comparisons
    pair_names = []
    pair_indices = []
    for i, j in combinations(range(n_group), 2):
        pair_names.append(f"{group_covars[j]}_{group_covars[i]}")
        pair_indices.append((j, i))

    # Also include individual group effects
    all_names = list(group_covars) + pair_names
    n_comp = len(all_names)

    beta_pair = np.zeros((n_tax, n_comp))
    var_pair = np.zeros((n_tax, n_comp))

    # Individual effects
    for k in range(n_group):
        beta_pair[:, k] = beta_hat_sub[:, k]
        var_pair[:, k] = var_hat[:, group_ind[k]]

    # Pairwise differences
    for k, (j, i) in enumerate(pair_indices):
        beta_pair[:, n_group + k] = beta_hat_sub[:, j] - beta_hat_sub[:, i]
        var_pair[:, n_group + k] = (
            vcov_group[:, j, j] + vcov_group[:, i, i] - 2.0 * vcov_group[:, j, i]
        )

    se_pair = np.sqrt(np.maximum(var_pair, 0))
    W_pair = beta_pair / se_pair

    # Apply mdFDR correction
    pval, qval = _mdfdr_pairwise(
        W=W_pair,
        dof=dof,
        fwer_ctrl=p_adjust,
        dmat=dmat,
        group=grouping,
        beta_hat=beta_hat,
        vcov_hat=vcov_group,
        alpha=alpha,
        dof_global=dof,
        group_indices=group_ind,
    )
    reject = qval <= alpha

    return {
        "beta": beta_pair,
        "se": se_pair,
        "W": W_pair,
        "p_val": pval,
        "q_val": qval,
        "reject": reject,
        "comp_names": all_names,
    }


def _mdfdr_pairwise(
    W,
    dof,
    fwer_ctrl,
    dmat,
    group,
    beta_hat,
    vcov_hat,
    alpha,
    dof_global=None,
    group_indices=None,
):
    """Perform mixed directional FDR (mdFDR) correction for pairwise tests.

    1. Run global test to screen significant features (count R).
    2. Only consider R significant features for pairwise p-values.
    3. Adjust at level R * alpha / d.

    """
    n_feats, n_comps = W.shape

    # Calculate p-values
    p_val = _calc_pvalues(W, dof)

    # Screen for significant comparisons using global test
    _, _, _, signif = _global_test(
        dmat=dmat,
        grouping=group,
        beta_hat=beta_hat,
        vcov_hat=vcov_hat,
        p_adjust="BH",  # TODO: Question: Is "BH" hard-coded? Not inherit?
        alpha=alpha,
        dof=dof_global,
        group_indices=group_indices,
    )
    n_signs = signif.sum().item()  # R

    # Zero out p-values for features not significant in global test
    p_val = p_val * signif[:, None]
    p_val[p_val == 0] = 1.0
    p_val = np.where(np.isnan(p_val), 1.0, p_val)

    # Apply mdFDR correction. `n_rejs` adjusts each feature's comparisons as one
    # family, with its size inflated by the global-screening factor
    # n_feats / n_rejs.
    func = _check_p_adjust(fwer_ctrl)
    if n_signs:
        n_tests = n_comps * n_feats // n_signs
        n_padding = n_tests - n_comps

        def adjust(pvals):
            # R's p.adjust(..., n=n_tests) is equivalent to adding n_tests -
            # len(pvals) unit p-values before correction.
            if n_padding:
                pvals = np.pad(pvals, (0, n_padding), constant_values=1.0)
            return func(pvals)[:n_comps]

        q_val = np.apply_along_axis(adjust, 1, p_val)
    else:
        q_val = np.ones_like(p_val)

    return p_val, q_val


def _var_diff(vcov_sub):
    """Compute variance of difference for a pair from a covariance sub-matrix.

    For a 2x2 covariance matrix [[var_a, cov_ab], [cov_ab, var_b]],
    var(a - b) = var_a + var_b - 2*cov_ab = trace - 2*off_diag_sum.

    For a general kxk matrix, var(sum of signed differences) =
    sum of all elements with appropriate signs. For the simple pairwise
    case with contrast [1, -1], this is trace - 2*sum(off-diagonal).
    """
    vcov_sub = np.asarray(vcov_sub)
    # For the pairwise case: var(a-b) = var(a) + var(b) - 2*cov(a,b)
    # = sum(diag) - 2*sum(off-diagonal upper/lower)
    # = sum(all elements) - 2*sum(off-diagonal)
    # More precisely: c' V c where c = [1, -1] for 2x2
    if vcov_sub.shape == (2, 2):
        return vcov_sub[0, 0] + vcov_sub[1, 1] - 2 * vcov_sub[0, 1]
    else:
        # General case: for contrast vector c, var = c' V c
        # Default contrast is [1, -1, 0, ...] etc.
        off_diag = vcov_sub.sum() - np.trace(vcov_sub)
        return np.trace(vcov_sub) - off_diag


def _dunnett_test(
    dmat,
    grouping,
    beta_hat,
    var_hat,
    dof,
    bootstraps,
    rng,
    p_adjust,
    alpha,
    group_indices=None,
):
    """ANCOM-BC2 Dunnett's type of test.

    Compare each group to the reference group with mdFDR correction.
    """
    covariates = dmat.design_info.column_names
    if group_indices is None:
        group_ind = _group_indices_from_design(dmat, grouping)
    else:
        group_ind = np.asarray(group_indices, dtype=int)

    beta_hat_dunn = beta_hat[:, group_ind]
    var_hat_dunn = var_hat[:, group_ind]
    se_hat_dunn = np.sqrt(np.maximum(var_hat_dunn, 0))
    W_dunn = beta_hat_dunn / se_hat_dunn

    # mdFDR correction
    p_val, q_val = _mdfdr_dunnett(
        W=W_dunn,
        dof=dof,
        fwer_ctrl=p_adjust,
        dmat=dmat,
        group=grouping,
        bootstraps=bootstraps,
        rng=rng,
        alpha=alpha,
    )

    return {
        "beta": beta_hat_dunn,
        "se": se_hat_dunn,
        "W": W_dunn,
        "p_val": p_val,
        "q_val": q_val,
        "reject": q_val <= alpha,
        "comp_names": [covariates[i] for i in group_ind],
    }


def _mdfdr_dunnett(W, dof, fwer_ctrl, dmat, group, bootstraps, alpha, rng):
    """mdFDR correction for Dunnett's test."""
    n_feats, n_comps = W.shape

    # Step 1: Global test screening via bootstrap
    res_screen = _dunn_global(
        dmat=dmat,
        group=group,
        W=W,
        bootstraps=bootstraps,
        dof=dof,
        p_adjust="BH",
        alpha=alpha,
        rng=rng,
    )
    # TODO: Likewise, should p_adjust be hard-coded as "BH"?
    n_signs = int(res_screen["reject"].sum())
    screen_ind = res_screen["reject"].values

    # Step 2: Compute p-values
    if dof is not None:
        dof = np.asarray(dof)
        p_val = 2 * t.sf(np.abs(W), df=dof if dof.ndim == 0 else dof[:, None])
    else:
        p_val = 2 * norm.sf(np.abs(W))

    # Zero out for non-significant taxa
    p_val = p_val * screen_ind[:, np.newaxis]
    p_val[p_val == 0] = 1.0
    p_val = np.where(np.isnan(p_val), 1.0, p_val)

    # Step 3: Adjust each feature's comparisons. R's p.adjust ``n`` parameter
    # is equivalent to padding each family with unit p-values.
    func = _check_p_adjust(fwer_ctrl)
    if n_signs:
        n_tests = n_comps * n_feats // n_signs
        n_padding = n_tests - n_comps

        def adjust(pvals):
            if n_padding:
                pvals = np.pad(pvals, (0, n_padding), constant_values=1.0)
            return func(pvals)[:n_comps]

        q_val = np.apply_along_axis(adjust, 1, p_val)
    else:
        q_val = np.ones_like(p_val)

    return p_val, q_val


def _dunn_global(dmat, group, W, bootstraps, dof, p_adjust, alpha, rng):
    """Dunnett's global test for mdFDR.

    Bootstrap-based: generate null W from t-distribution, take max |W|.
    """
    n_tax = W.shape[0]
    covariates = dmat.design_info.column_names
    group_ind = np.array([group in c and ":" not in c for c in covariates])
    n_group = int(np.sum(group_ind))

    # Observed global statistic: max |W| per taxon
    W_global = np.max(np.abs(W), axis=1)

    # Bootstrap null distribution
    W_global_null = np.zeros((n_tax, bootstraps))

    for b in range(bootstraps):
        # Generate null W from the per-feature t-distribution.
        if dof is not None:
            dof_null = np.nan_to_num(dof, nan=999.0)
            W_null = rng.standard_t(
                dof_null if np.ndim(dof_null) == 0 else dof_null[:, None],
                size=W.shape,
            )
        else:
            W_null = rng.standard_normal(size=W.shape)

        W_global_null[:, b] = np.max(np.abs(W_null), axis=1)

    # P-values from bootstrap
    p_global = np.mean(W_global_null > W_global[:, np.newaxis], axis=1)

    # R's p.adjust excludes NA values; ANCOM-BC2 then treats invalid global
    # tests as nonsignificant.
    q_global = _adjust_pvalues(p_global, p_adjust)
    q_global = np.where(np.isnan(q_global), 1.0, q_global)
    reject = q_global <= alpha

    return pd.DataFrame(
        {
            "W": W_global,
            "p_val": p_global,
            "q_val": q_global,
            "reject": reject,
        }
    )


def _trend_test(
    dmat,
    grouping,
    beta_hat,
    var_hat,
    vcov_hat,
    p_adjust="holm",
    alpha=0.05,
    trend_contrast=None,
    trend_node=None,
    bootstraps=100,
    rng=None,
    group_indices=None,
):
    """ANCOM-BC2 trend test (pattern analysis).

    Uses constrained optimization to test ordered patterns in group effects.
    """
    n_feats = beta_hat.shape[0]
    if group_indices is None:
        group_ind = _group_indices_from_design(dmat, grouping)
    else:
        group_ind = np.asarray(group_indices, dtype=int)
    n_group = group_ind.size

    beta_hat_sub = beta_hat[:, group_ind]
    var_hat_sub = var_hat[:, group_ind]
    vcov_hat_sub = _select_group_covariance(vcov_hat, group_ind)

    # Test both monotone directions. Each contrast operates on the non-reference
    # group coefficients, with the final coefficient as the trend node.
    if trend_contrast is None:
        increasing = np.eye(n_group)
        increasing[1:, :-1] -= np.eye(n_group - 1)
        trend_contrast = {
            "increasing": increasing,
            "decreasing": -increasing,
        }
        trend_node = {name: n_group - 1 for name in trend_contrast}

    n_trend = len(trend_contrast)
    trend_names = list(trend_contrast.keys())

    # Constrained estimation for each taxon and each trend pattern
    beta_hat_opt_all = np.zeros((n_feats, n_group * n_trend))
    l_vals = np.zeros((n_feats, n_trend))

    for i in range(n_feats):
        for t_idx, (tname, contrast) in enumerate(trend_contrast.items()):
            C = np.asarray(contrast)
            beta_opt = _constrain_est(beta_hat_sub[i], vcov_hat_sub[i], C)
            start_col = t_idx * n_group
            beta_hat_opt_all[i, start_col : start_col + n_group] = beta_opt

            # Compute l_infinity norm
            node = trend_node[tname]
            l_vals[i, t_idx] = _l_infty(beta_opt, node)

    # W_trend = max l_infinity across patterns
    W_trend = np.max(l_vals, axis=1)
    opt_trend_idx = np.argmax(l_vals, axis=1)

    # Select the optimal trend's beta for each taxon
    beta_hat_trend = np.zeros((n_feats, n_group))
    for i in range(n_feats):
        t_idx = opt_trend_idx[i]
        start_col = t_idx * n_group
        beta_hat_trend[i] = beta_hat_opt_all[i, start_col : start_col + n_group]

    # Bootstrap null distribution
    if rng is None:
        rng = np.random.default_rng()
    W_trend_null = np.zeros((n_feats, bootstraps))

    var_hat_sub_dup = np.nan_to_num(var_hat_sub, nan=1.0)

    for b in range(bootstraps):
        # Generate null beta from N(0, I)
        beta_null = rng.standard_normal(size=(n_feats, n_group))

        # Constrained estimation under null
        l_null = np.zeros((n_feats, n_trend))
        for t_idx, (tname, contrast) in enumerate(trend_contrast.items()):
            beta_null_opt = _constrain_est_identity(beta_null, contrast)
            beta_null_opt *= np.sqrt(np.maximum(var_hat_sub_dup, 0))
            node = trend_node[tname]
            l_null[:, t_idx] = np.maximum(
                np.abs(beta_null_opt[:, node]),
                np.abs(beta_null_opt[:, node] - beta_null_opt[:, -1]),
            )

        W_trend_null[:, b] = np.max(l_null, axis=1)

    # P-values from bootstrap
    p_trend = np.mean(W_trend_null > W_trend[:, np.newaxis], axis=1)

    # R's p.adjust excludes NA values; ANCOM-BC2 then treats invalid trend
    # tests as nonsignificant.
    q_trend = _adjust_pvalues(p_trend, p_adjust)
    q_trend = np.where(np.isnan(q_trend), 1.0, q_trend)
    diff_trend = q_trend <= alpha

    return {
        "beta": beta_hat_trend,
        "se": np.sqrt(np.maximum(var_hat_sub, 0)),
        "W": W_trend,
        "p_val": p_trend,
        "q_val": q_trend,
        "reject": diff_trend,
    }


def _constrain_est(beta_hat, vcov_hat, contrast):
    """Constrained estimation via quadratic programming.

    Minimize (beta - beta_opt)' P (beta - beta_opt) subject to
    contrast @ beta_opt >= 0 (each row of contrast defines one inequality).

    Uses scipy.optimize.minimize with SLSQP method, replacing R's quadprog::solve.QP.

    """
    P = _safe_inverse_spd(vcov_hat)
    n = len(beta_hat)
    contrast = np.asarray(contrast)

    # Objective: 0.5 * x' D x - d' x  (quadratic form)
    Dmat = 2 * P
    Dmat = (Dmat + Dmat.T) / 2
    # Add small ridge for numerical stability (matches R code: diag(Dmat) += 1e-10)
    Dmat += 1e-10 * np.eye(n)
    dvec = 2 * P @ beta_hat

    def objective(x):
        return 0.5 * x @ Dmat @ x - dvec @ x

    def grad(x):
        return Dmat @ x - dvec

    # Each row of contrast defines one inequality: contrast[i] @ x >= 0
    # SLSQP expects a single constraint function returning a vector
    n_constraints = contrast.shape[0]
    constraints = {
        "type": "ineq",
        "fun": lambda x: contrast @ x,  # Returns vector; each element must be >= 0
        "jac": lambda x: contrast,
    }

    try:
        result = minimize(
            objective,
            x0=beta_hat,
            jac=grad,
            method="SLSQP",
            constraints=constraints,
            options={"ftol": 1e-12, "maxiter": 1000},
        )
        return result.x
    except Exception:
        return np.zeros(n)


def _constrain_est_identity(beta_hat, contrast):
    """Project coefficient vectors onto linear inequalities under identity covariance.

    This is the identity-covariance specialization of :func:`_constrain_est` used by
    the trend bootstrap. It enumerates active constraint sets and evaluates all
    feature vectors simultaneously, avoiding repeated SLSQP calls.
    """
    beta_hat = np.asarray(beta_hat)
    contrast = np.asarray(contrast)
    n_features, n_coefficients = beta_hat.shape
    n_constraints = contrast.shape[0]

    # Active-set enumeration is efficient for the small ordinal contrasts used
    # by trend tests, but grows exponentially for larger custom systems.
    if n_constraints > 10:
        identity = np.eye(n_coefficients)
        return np.array([_constrain_est(beta, identity, contrast) for beta in beta_hat])

    candidates = []
    valid = []

    for size in range(n_constraints + 1):
        for active in combinations(range(n_constraints), size):
            if active:
                active_contrast = contrast[list(active)]
                gram = active_contrast @ active_contrast.T
                multipliers = -(beta_hat @ active_contrast.T @ np.linalg.pinv(gram))
                candidate = beta_hat + multipliers @ active_contrast
                active_valid = np.all(multipliers >= -1e-10, axis=1)
            else:
                candidate = beta_hat
                active_valid = np.ones(n_features, dtype=bool)

            candidates.append(candidate)
            valid.append(
                active_valid & np.all(contrast @ candidate.T >= -1e-10, axis=0)
            )

    candidates = np.asarray(candidates)
    valid = np.asarray(valid)
    distances = np.sum((candidates - beta_hat) ** 2, axis=2)
    distances[~valid] = np.inf
    selected = np.argmin(distances, axis=0)
    return candidates[selected, np.arange(n_features)]


def _safe_inverse_spd(A, ridge=1e-8):
    """Safe inverse of a symmetric positive-definite matrix."""
    A = (A + A.T) / 2
    try:
        L = np.linalg.cholesky(A)
        return np.linalg.solve(A, np.eye(A.shape[0]))
    except np.linalg.LinAlgError:
        return np.linalg.inv(A + ridge * np.eye(A.shape[0]))


def _l_infty(beta_opt, node):
    """Compute l_infinity norm for a pattern at specified node.

    Matches R's .l_infty function:
    l = max(abs(beta_opt[node]), abs(beta_opt[node] - beta_opt[length(beta_opt)]))

    Note: node is 0-based in Python (1-based in R). The last element
    beta_opt[-1] corresponds to beta_opt[length(beta_opt)] in R.
    """
    beta_opt = np.asarray(beta_opt)
    return max(
        abs(beta_opt[node]),
        abs(beta_opt[node] - beta_opt[-1]),
    )


def struc_zero(table, metadata, grouping, neg_lb=False):
    r"""Identify features with structural zeros.

    .. versionadded:: 0.7.1

    Structural zeros refer to features that are systematically absent from certain
    sample groups. Consequently, the observed feature frequencies are all zeros, or
    mostly zeros, due to variability in technical factors. This function tests
    whether the proportion of observed zeros is close to zero, which suggests the
    absence of a feature in a given sample group.

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing count or proportional abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    metadata : pd.DataFrame or 2-D array_like
        Metadata of the samples. Rows correspond to samples and columns correspond
        to covariates (attributes). Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    grouping : str
        A metadata column name indicating the assignment of samples to groups.
    neg_lb : bool, optional
        Determine whether to use negative lower bound when calculating sample
        proportions. Default is False. Generally, it is recommended to set it as True
        when the sample size per group is relatively large.

    Returns
    -------
    pd.DataFrame of bool of shape (n_features, n_groups)
        A table indicating whether each feature (row) is a structural zero in each
        group (column) (True: structural zero, False: not structural zero).

    Notes
    -----
    The structural zero test was initially proposed and implemented in the ANCOM-II
    method [1]_. It was adopted to the ANCOM-BC method [2]_ as a recommended method to
    complement test results. See :func:`ancombc` for how to use this function along
    with the ANCOM-BC test. Nevertheless, this function is generally useful with or
    without explicit statistical tests of feature abundances.

    A feature found to be a structural zero in a group should be automatically
    considered as differentially (less) abundant compared with other groups in which
    this feature is not a structural zero. Meanwhile, this feature should be excluded
    from subsequent analyses that involves this group. If a feature is identified as a
    structural zero in all groups, this feature should be removed entirely from
    downstream analyses.

    Note that the structural zero test should be applied to the original table before
    adding a pseudocount, which will otherwise mask all zeros and invalidate this test.

    References
    ----------
    .. [1] Kaul, A., Mandal, S., Davidov, O., & Peddada, S. D. (2017). Analysis of
       microbiome data in the presence of excess zeros. Frontiers in Microbiology, 8,
       2114.

    .. [2] Lin, H. and Peddada, S.D., 2020. Analysis of compositions of microbiomes
       with bias correction. Nature Communications, 11(1), p.3514.

    Examples
    --------
    >>> from skbio.stats.composition import struc_zero
    >>> import pandas as pd

    Generate a DataFrame with 10 samples and 6 features with 0's in specific groups:

    >>> table = pd.DataFrame([[ 7,  1,  0, 11,  3,  1],
    ...                       [ 1,  1,  0, 13, 13,  0],
    ...                       [11,  5,  0,  1,  4,  1],
    ...                       [ 2,  2,  0, 16,  4,  0],
    ...                       [ 0,  1,  0,  0,  6,  0],
    ...                       [14,  8,  7,  9,  0,  5],
    ...                       [ 0,  7,  4,  1,  0, 26],
    ...                       [ 8,  1,  4, 28,  0, 10],
    ...                       [ 2,  2,  2,  4,  0,  5],
    ...                       [ 6,  4, 10,  1,  0,  9]],
    ...                      index=[f's{i}' for i in range(10)],
    ...                      columns=[f'f{i}' for i in range(6)])

    Then create a grouping vector. In this example, there is a treatment group
    and a placebo group.

    >>> metadata = pd.DataFrame(
    ...     {'grouping': ['treatment'] * 5 + ['placebo'] * 5},
    ...     index=[f's{i}' for i in range(10)])

    The ``struc_zero`` function will identify features with structural zeros. Features
    that are identified as structural zeros in given groups should not be used in
    further analyses such as ``ancombc`` and  ``dirmult_ttest``.

    Setting ``neg_lb=True`` declares that the true prevalence of a feature in a group
    is not significantly different from zero.

    >>> result = struc_zero(table, metadata, grouping='grouping', neg_lb=True)
    >>> result
        placebo  treatment
    f0    False      False
    f1    False      False
    f2    False       True
    f3    False      False
    f4     True      False
    f5    False       True

    """
    # Validate feature table and metadata
    matrix, samples, features = _ingest_table(table)
    metadata = _check_metadata(metadata, matrix, samples)
    metadata = _type_cast_to_float(metadata)

    unique_groups, group_indices, group_counts = np.unique(
        metadata[grouping], return_inverse=True, return_counts=True
    )

    # Create a boolean matrix to indicate whether the value in table is 0 or not
    tmp = np.nan_to_num(matrix) != 0

    n_groups = len(unique_groups)
    n_features = tmp.shape[1]

    # Initialize group sum matrix
    group_sums = np.zeros((n_groups, n_features))
    np.add.at(group_sums, group_indices, tmp.astype(int))

    # Calculate sample sizes of the groups for each feature
    sample_size = group_counts[:, np.newaxis]

    # Calculate sample proportions of the groups for each feature
    p_hat = group_sums / sample_size

    # Calculate the lower bound of a 95% confidence interval for sample proportion
    if neg_lb:
        p_hat = p_hat - 1.96 * (p_hat * (1 - p_hat) / sample_size) ** 0.5

    zero_idx = p_hat <= 0

    # Output structural zero as a DataFrame
    return pd.DataFrame(zero_idx.T, index=features, columns=unique_groups)
