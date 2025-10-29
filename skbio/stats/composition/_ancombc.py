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

import numpy as np
import pandas as pd
from scipy.stats import norm, chi2
from scipy.optimize import minimize
from patsy import dmatrix

from skbio.table._tabular import _ingest_table
from ._base import _check_composition
from ._utils import _check_metadata, _check_p_adjust, _type_cast_to_float


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

    Parameters
    ----------
    table : table_like of shape (n_samples, n_features)
        A matrix containing strictly positive count or proportional abundance data of
        the samples. See :ref:`supported formats <table_like>`.

        .. note::
            If the table contains zero values, one should add a pseudocount or apply
            :func:`multi_replace` to convert all values into positive numbers.

    metadata : pd.DataFrame or 2-D array_like
        Metadata of the samples. Rows correspond to samples and columns correspond
        to covariates (attributes). Must be a pandas DataFrame or convertible to a
        pandas DataFrame.
    formula : str or generic Formula object
        A formula defining the model using factors included in the metadata columns.
        Refer to `Patsy's documentation
        <https://patsy.readthedocs.io/en/latest/formulas.html>`_ on how to specify
        a formula.
    grouping : str, optional
        A metadata column name of interest for *global test*, which identifies features
        that are differentially abundant between at least two groups across three or
        more groups in that column. Must be one of the factors in ``formula``. Default
        is None, which skips global test.
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
    res_main : pd.DataFrame
        A table of features and covariates, their log-fold changes and other relevant
        statistics of the primary ANCOM-BC analysis.

        - ``FeatureID``: Feature identifier, i.e., dependent variable.

        - ``Covariate``: Covariate name, i.e., independent variable.

        - ``Log2(FC)``: Expected log2-fold change of abundance from the reference
          category to the covariate category defined in the formula. The value is
          expressed in the center log ratio (see :func:`clr`) transformed coordinates.

        - ``SE``: Standard error of the estimated Log2(FC).

        - ``W``: *W*-statistic, or the number of features that the current feature is
          tested to be significantly different against.

        - ``pvalue``: *p*-value of the ANCOM-BC test.

        - ``qvalue``: Corrected *p*-value of the ANCOM-BC test for multiple
          comparisons.

        - ``Signif``: Whether the covariate category is significantly differentially
          abundant from the reference category. A feature-covariate pair is marked as
          "True" if the *q*-value is less than or equal to the significance level
          (``alpha``).

    res_global : pd.DataFrame, optional
        A table of features and statistics from the global test (when ``grouping`` is
        set).

        - ``FeatureID``: Feature identifier, i.e., dependent variable.

        - ``W``: *W*-statistic of Chi-squared statistic that quantifies the overall
          evidence against null hypothesis (the mean abundance of the feature is the
          same across all groups).

        - ``pvalue``: *p*-value of the *W*-statistic in global test.

        - ``qvalue``: Corrected *p*-value in global test.

        - ``Signif``: Whether at least one group mean is different from other
          in grouping variable.

    See Also
    --------
    ancom
    struc_zero
    multi_replace

    Notes
    -----
    This function is a Python re-implementation of the ANCOM-BC method [1]_, which was
    originally implemented in the R package ``ANCOMBC``. This function provides an
    efficient and scalable algorithm, with a simple interface consistent with other
    scikit-bio components. The output of this function should match that of the R
    package.

    Comparing with the R command ``ancombc``, this function defers the flexibility of
    data preprocessing to the user. Most importantly, if the input table contains zero
    values (which is common), one needs to remove them by, e.g., adding a pseudocount
    of 1 (``pseudo=1`` in the R command's parameters) (assuming ``table`` is a pandas
    DataFrame):

    .. code-block:: python

       table += 1

    See also :func:`multi_replace` for additional information on zero handling.

    Some other pre-processing options provided by the R command can be performed with:

    To aggregate data at a given taxonomic level (``tax_level="Family"``):

    .. code-block:: python

       table = table.T.groupby(feature_to_family_dict).sum().T

    To discard features with prevalence < 10% among samples (``prv_cut=0.1``):

    .. code-block:: python

       table = table.loc[:, (table > 0).mean() >= 0.1]

    To discard samples with a total abundance < 1 million (``lib_cut=1e6``):

    .. code-block:: python

       table = table.loc[table.sum(axis=1) >= 1e6]

    Categorical columns in metadata are sorted alphabetically, and the reference level
    for each column is automatically set to be the first category. If this behavior is
    not intended, you will need to change the order manually, like:

    .. code-block:: python

       metadata['bmi'] = pd.Categorical(
           metadata['bmi'], categories=['lean', 'overweight', 'obese'], ordered=True)

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

    >>> result = ancombc(table + 1, metadata, formula='status')
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
    while controlling for the other. Additionally, we instruct the function to perform
    a *global test* on "status", which identifies features that are differentially
    abundant between at least two status.

    >>> res_main, res_global = ancombc(
    ...     table + 1, metadata, formula='status + age', grouping='status')
    >>> res_main.round(3)
                                  Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F1        Intercept             -0.006  0.346 -0.016   0.987   1.000   False
              status[T.moderate]    -0.337  0.245 -1.376   0.169   1.000   False
              status[T.severe]       0.278  0.350  0.792   0.428   1.000   False
              age                    0.002  0.011  0.141   0.888   1.000   False
    F2        Intercept              0.072  0.446  0.160   0.873   1.000   False
              status[T.moderate]     1.906  0.259  7.371   0.000   0.000    True
              status[T.severe]       2.935  0.304  9.649   0.000   0.000    True
              age                   -0.022  0.013 -1.676   0.094   0.562   False
    F3        Intercept             -1.690  0.570 -2.967   0.003   0.021    True
              status[T.moderate]     1.010  0.578  1.747   0.081   0.565   False
              status[T.severe]       0.717  0.517  1.388   0.165   1.000   False
              age                    0.063  0.022  2.883   0.004   0.031    True
    F4        Intercept              0.439  0.483  0.910   0.363   1.000   False
              status[T.moderate]    -0.046  0.405 -0.113   0.910   1.000   False
              status[T.severe]      -0.244  0.536 -0.455   0.649   1.000   False
              age                   -0.003  0.019 -0.185   0.854   1.000   False
    F5        Intercept             -1.227  0.393 -3.125   0.002   0.014    True
              status[T.moderate]     0.274  0.293  0.934   0.350   1.000   False
              status[T.severe]      -0.323  0.320 -1.011   0.312   1.000   False
              age                    0.027  0.014  2.021   0.043   0.303   False
    F6        Intercept             -0.439  0.440 -0.999   0.318   1.000   False
              status[T.moderate]     0.114  0.432  0.264   0.792   1.000   False
              status[T.severe]       0.375  0.544  0.691   0.490   1.000   False
              age                    0.008  0.016  0.480   0.631   1.000   False
    F7        Intercept              0.201  0.472  0.426   0.670   1.000   False
              status[T.moderate]     0.082  0.302  0.270   0.787   1.000   False
              status[T.severe]      -0.264  0.401 -0.658   0.510   1.000   False
              age                    0.004  0.016  0.272   0.785   1.000   False
    F8        Intercept             -0.168  0.526 -0.319   0.750   1.000   False
              status[T.moderate]    -0.402  0.584 -0.688   0.491   1.000   False
              status[T.severe]      -0.299  0.728 -0.411   0.681   1.000   False
              age                    0.022  0.018  1.237   0.216   1.000   False

    We found that feature "F2" is significantly differentially (more) abundant in
    "moderate" and "severe" groups compared with "mild", which serves as the reference
    group. Additionally, feature "F3" is separately correlated with age.

    >>> res_main.query('Covariate != "Intercept" & Signif == True').round(3)
                                  Log2(FC)     SE      W  pvalue  qvalue  Signif
    FeatureID Covariate
    F2        status[T.moderate]     1.906  0.259  7.371   0.000   0.000    True
              status[T.severe]       2.935  0.304  9.649   0.000   0.000    True
    F3        age                    0.063  0.022  2.883   0.004   0.031    True

    The global test result suggests that "F2" and "F4" are differentially abundant
    between two of the three groups (though it doesn't tell which groups).

    >>> res_global.round(3)
                    W  pvalue  qvalue  Signif
    FeatureID
    F1          1.855   0.791   1.000   False
    F2         80.771   0.000   0.000    True
    F3          2.925   0.463   1.000   False
    F4         -0.093   0.000   0.000    True
    F5          1.068   0.827   1.000   False
    F6          0.121   0.117   0.704   False
    F7          0.220   0.208   1.000   False
    F8          0.485   0.430   1.000   False

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
    F4     True
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
    # Note: A pseudocount should have been added to the table by the user prior to
    # calling this function.
    matrix, samples, features = _ingest_table(table)
    _check_composition(np, matrix, nozero=True)
    n_feats = matrix.shape[1]
    if features is None:
        features = np.arange(n_feats)

    # Log-transform count matrix.
    matrix = np.log(matrix)

    # Validate metadata and cast to numbers where applicable.
    metadata = _check_metadata(metadata, matrix, samples)
    metadata = _type_cast_to_float(metadata)

    if grouping is not None and len(metadata[grouping].unique()) < 3:
        raise ValueError("Global test cannot be performed on less than three groups.")

    # Create a design matrix based on metadata and formula.
    dmat = dmatrix(formula, metadata)

    # Obtain a list of covariates by selecting the relevant columns.
    covars = dmat.design_info.column_names
    n_covars = len(covars)

    # validate parameters
    if not 0 < alpha < 1:
        raise ValueError("`alpha`=%f is not within 0 and 1." % alpha)

    # Estimate initial model parameters.
    var_hat, beta, _, vcov_hat = _estimate_params(matrix, dmat)

    # Estimate and correct for sampling bias via expectation-maximization (EM).
    bias = np.empty((n_covars, 3))
    for i in range(n_covars):
        bias[i] = _estimate_bias_em(beta[i], var_hat[:, i], tol=tol, max_iter=max_iter)
    delta_em = bias[:, 0]

    # Correct coefficients (logFC) according to estimated bias.
    beta_hat = beta.T - delta_em

    # Calculate statistics while performing multiple testing correction.
    se_hat, W, pval, qval = _calc_statistics(beta_hat, var_hat, method=p_adjust)

    # Identify significantly differentially abundance feature-covariate pairs.
    reject = qval <= alpha

    # Output the primary results.
    res = pd.DataFrame.from_dict(
        {
            "FeatureID": [x for x in features for _ in range(n_covars)],
            "Covariate": list(covars) * n_feats,
            "Log2(FC)": beta_hat.ravel(),
            "SE": se_hat.ravel(),
            "W": W.ravel(),
            "pvalue": pval.ravel(),
            "qvalue": qval.ravel(),
        }
    )
    # Pandas' nullable boolean type
    res["Signif"] = pd.Series(reject.ravel(), dtype="boolean")
    res.set_index(["FeatureID", "Covariate"], inplace=True)

    # Output global test results
    if grouping is not None:
        res_global = _global_test(dmat, grouping, beta_hat, vcov_hat, alpha, p_adjust)
        res_global_df = pd.DataFrame.from_dict(
            {
                "FeatureID": features,
                "W": res_global[0],
                "pvalue": res_global[1],
                "qvalue": res_global[2],
                "Signif": res_global[3],
            }
        )
        res_global_df.set_index("FeatureID", inplace=True)
        return res, res_global_df
    return res  # , None


def _estimate_params(data, dmat):
    """Estimate initial model parameters.

    Perform initial estimation of model parameters (coefficients, variances and
    mean residuals) based on the observed data prior to bias correction.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Data table. Must be zero-handled and log-transformed.
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.

    Returns
    -------
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances of regression coefficients.
    beta : ndarray of shape (n_features, n_covariates)
        Estimated coefficients (log-fold changes before correction).
    theta : ndarray of shape (n_samples,)
        Per-sample mean residuals of estimated data.
    beta_covmat : ndarray of shape (n_features, n_covariates, n_covariates)
        Estimated covariance matrices.

    """
    # The original R code performs iterative maximum likelihood estimation to calculate
    # the coefficients, and calls MASS:ginv to calculate the pseudo inverse Gram matrix
    # using the Moore-Penrose method. We noticed that the former can be replaced with
    # ordinary least squares by NumPy's `lstsq`; and the latter matches the method by
    # NumPy's `pinv`. We further noticed that both functions perform singular value
    # decomposition (SVD). Therefore, the following code only performs SVD once, and
    # uses the intermediates to calculate coefficients and inverse Gram matrix.

    # Perform thin SVD and set small singular values to zero. The process and threshold
    # (1e-15) are consistent with the underlying algorithm of `pinv`.
    # Note: the Gram matrix may be singular, if there are colinear covariates in the
    # design matrix. In such cases, a direct `inv` will raise an error, and `pinv` is
    # the robust choice.
    U, S, Vh = np.linalg.svd(dmat, full_matrices=False)
    S_inv = np.where(S > 1e-15 * np.max(S), 1.0 / S, 0.0)

    # Regression coefficients
    V = Vh.T
    dmat_inv = (V * S_inv) @ U.T
    beta = dmat_inv @ data

    # Inverse Gram matrix
    gmat_inv = (V * S_inv**2) @ Vh

    # Per-sample mean residuals (theta)
    diff = data - dmat @ beta
    theta = np.mean(diff, axis=1, keepdims=True)

    # Centered residuals
    eps = diff - theta

    # Calculate the covariance matrix of the coefficients. The estimated variances are
    # the diagonal of the covariance matrix.
    # Note: The original R code uses nested `for` loops over samples and features. The
    # current implementation is fully vectorized.
    # Note: The original R code patches NaN with 0.1 when calculating covariances. This
    # is not needed in the current implementation, as the input data are guaranteed to
    # contain only finite real numbers.
    intm_mat = np.einsum("ji,jp,jq->ipq", eps**2, dmat, dmat)
    beta_covmat = (gmat_inv @ intm_mat) @ gmat_inv
    var_hat = np.diagonal(beta_covmat, axis1=1, axis2=2)

    # Note: Residuals are needed for ANCOM-BC2.
    return var_hat, beta, theta.reshape(-1), beta_covmat


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
    # in the outer function `ancombc` and re-used across covariates, which could
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
    stdevs = np.empty(shape)  # standard deviations
    ratios = np.empty(shape)  # coefficients / variances

    # Mean coefficients
    means = np.empty(3)

    # Posterior probabilities of feature-component assignments (EM's responsibilities)
    resp = np.empty(shape)

    # Just a 2-row array to store random data
    intm = np.empty((2, n_feats))

    # Initialize intermediates. The 1st row is constant, representing pre-correction
    # estimates, whereas the 2nd and 3rd rows are to be modified during iteration.
    np.reciprocal(var_hat, out=nu_inv[0])
    np.sqrt(var_hat, out=stdevs[0])
    np.divide(beta, var_hat, out=ratios[0])

    # Objective function for numeric optimization of variance estimation
    # Note: `norm.logpdf` doesn't have an `out` parameter. To further optimize this,
    # one needs to manually implement the under-the-hood algorithm.
    def func(x, loc, resp):
        log_pdf = norm.logpdf(beta, loc=loc, scale=(var_hat + x) ** 0.5)
        return -np.dot(resp, log_pdf)

    # Optimizer arguments. The Nelder-Mead simplex algorithm is used, which is
    # consistent with the original R code.
    # Note: The Nelder–Mead method doesn't actually enforce bounds during optimization.
    # It merely clips at the bounds. For further consideration.
    args = dict(method="Nelder-Mead", bounds=[(0, None)])

    # Expectation-maximization (E-M) iterations
    loss, epoch = np.inf, 0
    while loss > tol and epoch < max_iter:
        # Update intermediates (2nd and 3rd rows only)
        np.add(var_hat, params[6:8, None], out=intm)  # kappa1, kappa2
        np.reciprocal(intm, out=nu_inv[1:])
        np.sqrt(intm, out=stdevs[1:])
        np.subtract(beta, params[4:6, None], out=ratios[1:])  # means (l)
        ratios[1:] *= nu_inv[1:]

        ### E-step ###
        # Mean coefficients
        delta = means[0] = params[3]  # global bias (delta)
        np.add(params[4:6], delta, means[1:])

        # Posterior probabilities = mean probability density functions weighted by
        # component fractions
        # Note: `norm.pdf` doesn't have an `out` parameter. To further optimize this,
        # one needs to manually implement the under-the-hood algorithm.
        # p_r,i = (pi_r * phi(delta_i - (delta + l_r) / nu_ir)) /
        #         sum_r(pi_r * phi((delta_i - (dleta + l_r)) / nu_ir)),
        # where r = 0, 1, 2; i = 1, ..., n_features
        resp[:] = norm.pdf(beta, means[:, None], stdevs)
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
        # components (kappa).
        # TODO: Consider scenarios where optimization doesn't converge.
        updated[6] = minimize(func, params[6], args=(means[1], resp[1]), **args).x[0]
        updated[7] = minimize(func, params[7], args=(means[2], resp[2]), **args).x[0]

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

    # TODO: var_delta will be used if conserve=True to account for the variance of
    # delta_hat
    return delta_em, delta_wls, var_delta


def _sample_fractions(data, dmat, beta_hat):
    """Estimate sampling fractions.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        Data table. Zero-handled. Log-transformed.
    dmat : ndarray of shape (n_samples, n_covariates)
        Design matrix.
    beta_hat : ndarray of shape (n_features, n_covariates)
        Corrected coefficients.

    Returns
    -------
    theta_hat : ndarray of shape (n_samples,)
        Sampling fractions.

    """
    return np.mean(data - dmat @ beta_hat.T, axis=1)


def _calc_statistics(beta_hat, var_hat, method="holm"):
    """Calculate statistical significance while correcting for multiple testing.

    Parameters
    ----------
    beta_hat : ndarray of shape (n_features, n_covariates)
        Estimated coefficients post correction.
    var_hat : ndarray of shape (n_features, n_covariates)
        Estimated variances.

    Returns
    -------
    se_hat : ndarray of shape (n_features, n_covariates)
        Estimated standard errors.
    W : ndarray of shape (n_features, n_covariates)
        Test statistics.
    p : ndarray of shape (n_features, n_covariates)
        p-values.
    q : ndarray of shape (n_features, n_covariates)
        Adjusted p-values.

    """
    se_hat = var_hat**0.5
    W = beta_hat / se_hat
    pval = 2.0 * norm.sf(abs(W), loc=0, scale=1)
    func = _check_p_adjust(method)
    qval = np.apply_along_axis(func, 0, pval)
    return se_hat, W, pval, qval


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
    adding a pseudocount (see :func:`multi_replace`), which will otherwise mask all
    zeros and invalidate this test.

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


def _global_test(dmat, grouping, beta_hat, vcov_hat, alpha=0.05, p_adjust="holm"):
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
    vcov_hat : ndarray of shape (n_features, n_covariates, n_covariates)
        Estimated covariance matrices.
    alpha : float, optional
        Significance level for the statistical tests. Must be in the range of (0, 1).
        Default is 0.05.
    p_adjust : str, optional
        Method to correct *p*-values for multiple comparisons. Options are Holm-
        Boniferroni ("holm" or "holm-bonferroni") (default), Benjamini-
        Hochberg ("bh", "fdr_bh" or "benjamini-hochberg"), or any method supported
        by statsmodels' :func:`~statsmodels.stats.multitest.multipletests` function.
        Case-insensitive. If None, no correction will be performed.

    Returns
    -------
    W_global : ndarray of float of shape (n_features,)
        W-statistics of global test.
    pval : ndarray of float of shape (n_features,)
        p-values of global test.
    qval : ndarray of float of shape (n_features,)
        Adjusted p-values of global test.
    reject : ndarray of bool of shape (n_features,)
        If the variable is differentially abundant.

    """
    # Slices of columns in the dmat that the terms in the grouping is mapped to
    s = dmat.design_info.term_name_slices[grouping]

    # Get the index of the terms in the grouping
    group_ind = np.array(range(*s.indices(s.stop)))

    # Subset beta_hat and vcov_hat by grouping indices
    beta_hat_sub = beta_hat[:, group_ind]
    vcov_hat_sub = vcov_hat[:, group_ind][:, :, group_ind]

    # Inverse the subset of vcov_hat
    vcov_hat_sub_inv = np.linalg.pinv(vcov_hat_sub)

    dof = group_ind.size
    A = np.identity(dof)

    # for each feature, calculate test statistics W by the following formula:
    # W = (A @ beta_hat_sub).T @ inv(A @ vcov_hat_sub @ A.T) @ (A @ beta_hat_sub)
    term = np.einsum("ik,jk->ji", A, beta_hat_sub)
    W_global = np.einsum("ni,nij,ni->n", term, vcov_hat_sub_inv, term)

    # Derive p-values from W statistics
    p_lower = chi2.cdf(W_global, dof)
    p_upper = chi2.sf(W_global, dof)
    pval = 2 * np.minimum(p_lower, p_upper)

    # Correct p-values
    func = _check_p_adjust(p_adjust)
    qval = np.apply_along_axis(func, 0, pval)
    reject = qval <= alpha

    return W_global, pval, qval, reject
