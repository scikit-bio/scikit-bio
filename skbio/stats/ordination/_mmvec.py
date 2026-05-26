# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""MMvec - Microbe-Metabolite Vectors.

This module implements MMvec for learning joint embeddings of two omics modalities from
their co-occurrence patterns.

This implementation was adopted, modified and improved from the original mmvec package,
licensed under BSD-3-Clause:

- https://github.com/biocore/mmvec

"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scipy.sparse import coo_array

from skbio._base import SkbioObject
from skbio.stats.composition import clr_inv as softmax
from skbio.stats.composition import ilr_inv
from skbio.util import get_rng
from skbio.table._tabular import _ingest_table, _create_table, _create_table_1d

if TYPE_CHECKING:  # pragma: no cover
    from skbio.util._typing import SeedLike, TableLike


def mmvec(
    X: TableLike,
    Y: TableLike,
    dimensions: int = 3,
    optimizer: str = "lbfgs",
    max_iter: int = 1000,
    x_prior_mean: float = 0.0,
    x_prior_scale: float = 1.0,
    y_prior_mean: float = 0.0,
    y_prior_scale: float = 1.0,
    learning_rate: float = 1e-3,
    batch_size: int = 50,
    beta_1: float = 0.9,
    beta_2: float = 0.95,
    clipnorm: float = 10.0,
    batch_norm: str = "unbiased",
    seed: SeedLike | None = None,
    verbose: bool = False,
    output_format: str | None = None,
) -> MMvecResult:
    r"""Perform multiomics Microbe-Metabolite Vectors (MMvec) analysis.

    MMvec learns joint embeddings of two feature sets from their co-occurrence patterns
    using a multinomial likelihood model [1]_. The model learns:

    .. math::

        P(Y_j | X_i) = \text{softmax}(\hat{X}_i \cdot \hat{Y}_j +
        b_{\hat{X}_i} + b_{\hat{Y}_j})

    where :math:`\hat{X}_i` is the learned embedding vector of feature :math:`i` in the
    X modality (e.g., microbiome), :math:`\hat{Y}_j` is the learned embedding vector of
    feature :math:`j` in the Y modality (e.g., metabolome), and :math:`b_{\hat{X}_i}`,
    :math:`b_{\hat{Y}_j}` are learned bias terms.

    While MMvec was originally developed to analyze microbe-metabolite co-occurrence,
    this method is **generic** and can be applied to any two modalities that can be
    represented as compositional (count-based) data and share the same samples.

    .. versionadded:: 0.7.3

    Parameters
    ----------
    X : table_like of shape (n_samples, n_features_x)
        Abundance counts for the first modality (e.g., microbes). This modality is
        treated as the "conditioning" variable. See
        :ref:`supported formats <table_like>`.

    Y : table_like of shape (n_samples, n_features_y)
        Abundance counts for the second modality (e.g., metabolites). This modality is
        treated as the "conditioned" variable. See above. Must have the same samples as
        ``X``.

    dimensions : int, optional
        Number of latent dimensions for embeddings. Default is 3.

    optimizer : {'lbfgs', 'adam'}, optional
        Optimization algorithm to use. Default is 'lbfgs'.

        - 'lbfgs': L-BFGS-B quasi-Newton method. Recommended for most cases. Typically
          converges in 50-200 iterations. Deterministic.
        - 'adam': Stochastic gradient descent with Adam. Use for very large datasets or
          when stochastic behavior is desired.
    max_iter : int, optional
        Maximum number of iterations. Default is 1000. For 'lbfgs', this is the max
        number of L-BFGS iterations. For 'adam', this is the number of epochs.
    x_prior_mean : float, optional
        Mean of Gaussian prior on first modality embeddings. Default is 0.0.
    x_prior_scale : float, optional
        Scale (std) of Gaussian prior on first modality embeddings. Default is 1.0.
        Smaller values increase regularization.
    y_prior_mean : float, optional
        Mean of Gaussian prior on second modality embeddings. Default is 0.0.
    y_prior_scale : float, optional
        Scale (std) of Gaussian prior on second modality embeddings. Default is 1.0.
        Smaller values increase regularization.

    learning_rate : float, optional
        Adam learning rate. Ignored for 'lbfgs'. Default is 1e-3.
    batch_size : int, optional
        Adam mini-batch size. Ignored for 'lbfgs'. Default is 50.
    beta_1 : float, optional
        Adam exponential decay rate for first moment. Ignored for 'lbfgs'.
        Default is 0.9.
    beta_2 : float, optional
        Adam exponential decay rate for second moment. Ignored for 'lbfgs'.
        Default is 0.95.
    clipnorm : float, optional
        Adam gradient clipping threshold (global L2 norm). Ignored for 'lbfgs'.
        Default is 10.0.
    batch_norm : {'unbiased', 'legacy'}, optional
        Method for scaling mini-batch likelihood in Adam. Ignored for 'lbfgs'.

        - 'unbiased' (default): Uses norm = sum(n_features_x) / batch_size.
        - 'legacy': Uses norm = n_samples / batch_size.

    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.
    verbose : bool, optional
        Print training progress. Default is False.
    output_format : str, optional
        Output table format. See :ref:`table_params` for details.

    Returns
    -------
    :class:`MMvecResult`
        Result of MMvec model fitting, including learned embeddings, conditional
        probabilities, and convergence history. It also provides methods for predicting
        target distributions, and scoring predictive performance.

    See Also
    --------
    MMvecResult

    References
    ----------
    .. [1] Morton, J. T., Aksenov, A. A., Nothias, L. F., Foulds, J. R., Quinn, R. A.,
       Badri, M. H., ... & Knight, R. (2019). Learning representations of
       microbe-metabolite interactions. Nature Methods, 16(12), 1306-1314.

    Examples
    --------
    >>> from skbio.stats.ordination import mmvec

    Create a toy example of paired microbiome and metabolome datasets.

    >>> import pandas as pd
    >>> samples = ['S1', 'S2', 'S3', 'S4', 'S5']
    >>> microbes = ['O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8']
    >>> microbiome = pd.DataFrame([
    ...     [1, 0, 0, 7, 0, 2, 1, 0],
    ...     [0, 1, 2, 9, 1, 1, 0, 1],
    ...     [2, 0, 1, 6, 0, 0, 2, 0],
    ...     [1, 2, 1, 8, 3, 1, 0, 2],
    ...     [0, 1, 3, 3, 1, 4, 1, 0],
    ... ], index=samples, columns=microbes)
    >>> metabolites = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6']
    >>> metabolome = pd.DataFrame([
    ...     [0, 1, 0, 1, 1, 2],
    ...     [1, 0, 2, 1, 5, 0],
    ...     [0, 3, 1, 0, 1, 1],
    ...     [3, 1, 0, 2, 2, 0],
    ...     [0, 0, 4, 1, 7, 3],
    ... ], index=samples, columns=metabolites)

    Fit the MMvec model to the data. This will perform numerical optimization for up to
    100 iterations to learn the embeddings and co-occurrence patterns.

    >>> res = mmvec(microbiome, metabolome, max_iter=100, seed=42)

    Features of the two datasets (i.e., microbes and metabolites) are projected to the
    same latent space. Their coordinates (i.e., embeddings) are stored in two matrices:
    ``x_embeddings`` and ``y_embeddings``, respectively. The last column of each matrix
    is the bias term.

    >>> res.x_embeddings.round(3)
          PC1    PC2    PC3   bias
    O1 -0.442 -0.358  0.701  0.036
    O2  0.136  0.378 -0.280 -0.304
    O3  0.226 -0.206 -0.668  0.376
    O4  0.036 -0.015  0.088  0.055
    O5  0.137  0.490 -0.068 -0.423
    O6 -0.268 -0.047 -0.899  0.622
    O7  0.009 -0.924  0.025  1.099
    O8  0.575  0.823  0.398 -0.593

    >>> res.y_embeddings.round(3)
          PC1    PC2    PC3   bias
    C1  0.000  0.000  0.000  0.000
    C2 -0.153 -0.861  0.936 -0.149
    C3  0.186 -0.475 -0.742  0.102
    C4 -0.079  0.429 -0.001  0.020
    C5  0.266 -0.074 -0.541  1.028
    C6 -0.743 -0.965 -0.525 -0.132

    One can generate a biplot visualizing the co-occurrence patterns of microbes and
    metabolites. The distance and angle between features indicate their association
    strength. Note that the origin of the plot alway lies on the first metabolite in
    the dataset, which serves as a reference.

    >>> import matplotlib.pyplot as plt
    >>> plt.scatter('PC1', 'PC2', data=res.x_embeddings)  # doctest: +SKIP
    >>> plt.scatter('PC1', 'PC2', data=res.y_embeddings)  # doctest: +SKIP
    >>> plt.grid()  # doctest: +SKIP

    The ``ranks`` matrix stores the log conditional probabilities of microbe-metabolite
    co-occurrences. Larger positive values indicate a stronger likelihood of
    co-occurrence. Low and negative values indicate no relationship, not necessarily a
    negative correlation. Values are row-centered. One can sort each row to reveal
    which metabolites are most strongly associated with each microbe.

    >>> res.ranks.round(3)
           C1     C2     C3     C4     C5     C6
    O1 -0.228  0.691 -0.522 -0.291  0.366 -0.017
    O2  0.202 -0.859  0.054  0.070  1.086 -0.553
    O3 -0.602 -0.858  0.511 -0.313  1.239  0.023
    O4 -0.180 -0.184 -0.074 -0.114  0.866 -0.314
    O5  0.368 -0.710 -0.109  0.165  1.011 -0.725
    O6 -0.833 -1.119  0.531 -0.189  1.236  0.373
    O7 -1.355  0.411  0.266 -0.634  0.829  0.483
    O8  0.725 -0.441 -0.347  0.459  1.036 -1.431

    To obtain actual conditional probabilities, apply inverse CLR transform (a.k.a.
    softmax) to the ranks matrix.

    >>> from skbio.stats.composition import clr_inv
    >>> probs = clr_inv(res.ranks)
    >>> probs.round(3)
    array([[ 0.121,  0.304,  0.09 ,  0.114,  0.22 ,  0.15 ],
           [ 0.167,  0.058,  0.144,  0.147,  0.405,  0.079],
           [ 0.07 ,  0.054,  0.212,  0.093,  0.44 ,  0.13 ],
           [ 0.127,  0.126,  0.141,  0.135,  0.361,  0.111],
           [ 0.2  ,  0.068,  0.124,  0.163,  0.379,  0.067],
           [ 0.053,  0.04 ,  0.208,  0.101,  0.421,  0.177],
           [ 0.034,  0.201,  0.174,  0.071,  0.305,  0.216],
           [ 0.256,  0.08 ,  0.088,  0.197,  0.35 ,  0.03 ]])

    The ``convergence`` vector stores the loss curve over iterations during model
    training. If the model has converged, the curve should decrease and stabilize.

    >>> plt.plot(res.convergence)  # doctest: +SKIP

    A fitted MMvec model can be used to predict metabolite abundances given new
    microbiome data using the ``predict`` method.

    >>> samples_test = ['S6', 'S7', 'S8']
    >>> microbiome_test = pd.DataFrame([
    ...     [3, 0, 2, 6, 1, 4, 4, 0],
    ...     [2, 1, 5, 0, 0, 3, 6, 1],
    ...     [8, 4, 1, 1, 5, 0, 5, 1],
    ... ], index=samples_test, columns=microbes)
    >>> metabolome_pred = res.predict(microbiome_test)
    >>> metabolome_pred.round(3)
           C1     C2     C3     C4     C5     C6
    S6  0.091  0.140  0.160  0.109  0.349  0.151
    S7  0.077  0.130  0.175  0.098  0.360  0.160
    S8  0.131  0.171  0.129  0.124  0.318  0.128

    To evaluate the model's generalizability, one can compute the :math:`Q^2` statistic
    on held-out test data using the ``score`` method. A positive :math:`Q^2` indicates
    better-than-mean predictive performance (with 1.0 indicating perfect prediction). A
    negative or close-to-zero :math:`Q^2` is a sign of model overfitting or poor
    generalization.

    >>> metabolome_test = pd.DataFrame([
    ...     [0, 5, 1, 1, 8, 0],
    ...     [1, 0, 5, 3, 2, 3],
    ...     [4, 1, 0, 1, 6, 3],
    ... ], index=samples_test, columns=metabolites)
    >>> q2 = res.score(microbiome_test, metabolome_test)
    >>> float(q2.round(5))
    0.02291

    """
    estimator = MMvec(
        n_components=dimensions,
        optimizer=optimizer,
        max_iter=max_iter,
        x_prior_mean=x_prior_mean,
        x_prior_scale=x_prior_scale,
        y_prior_mean=y_prior_mean,
        y_prior_scale=y_prior_scale,
        learning_rate=learning_rate,
        batch_size=batch_size,
        beta_1=beta_1,
        beta_2=beta_2,
        clipnorm=clipnorm,
        batch_norm=batch_norm,
        seed=seed,
        verbose=verbose,
        output_format=output_format,
    )
    fitted = estimator.fit(X, Y)
    return MMvecResult(fitted)


class MMvecResult:
    r"""Result of an MMvec analysis.

    This class contains the learned embeddings and co-occurrence patterns from fitting
    an MMvec model. It enables both interpretation (which conditioning (X) and
    conditioned (Y) features co-occur with each other) and prediction (expected Y
    composition given X composition).

    Attributes
    ----------
    x_embeddings : table_like of shape (n_features_x, n_dimensions + 1)
    y_embeddings : table_like of shape (n_features_y, n_dimensions + 1)
    ranks : table_like of shape (n_features_x, n_features_y)
    convergence : table_like of shape (n_iterations,)

    See Also
    --------
    mmvec

    Notes
    -----
    **Detecting overfitting with Q-squared**

    Overfitting occurs when the model memorizes training data rather than learning
    generalizable patterns. To detect overfitting:

    1. Split your data into training and test sets before fitting.
    2. Fit the model on training data only.
    3. Use :meth:`score` to compute :math:`Q^2` on held-out test data.

    Interpretation of :math:`Q^2` values:

    - **Close to 1**: Excellent predictive performance.
    - **Close to 0**: Model predicts no better than the mean.
    - **Negative**: Model performs worse than predicting the mean,
      indicating overfitting or model misspecification.

    If :math:`Q^2` is much lower than expected, try:

    - Reducing ``dimensions`` (fewer latent dimensions).
    - Increasing regularization via smaller ``x_prior_scale`` and ``y_prior_scale``
      values.
    - Collecting more training samples.

    **Embedding Interpretation**

    The embeddings place X and Y features in the same latent space. The inner product
    between an X embedding vector and a Y embedding vector (plus their bias terms)
    gives the log-odds of their co-occurrence:

    .. math::

        \log \frac{P(m_j | \mu_i)}{P(m_{\text{ref}} | \mu_i)} =
        X_i \cdot Y_j + b_{X_i} + b_{Y_j}

    This means that X and Y features pointing in similar directions associate with each
    other, and the angle between a pair of X and Y vectors indicates their association
    strength.

    """

    def __init__(self, estimator: MMvec):
        self._estimator = estimator

    def __str__(self) -> str:
        n_features_x = self.x_embeddings.shape[0]
        n_features_y = self.y_embeddings.shape[0]
        n_dimensions = self.x_embeddings.shape[1] - 1
        return (
            f"MMvecResult\n"
            f"  Features: {n_features_x}\n"
            f"  Targets: {n_features_y}\n"
            f"  Components: {n_dimensions}\n"
            f"  Iterations: {len(self.convergence)}"
        )

    @property
    def x_embeddings(self) -> TableLike:
        r"""Learned coordinates of conditioning features (X) in latent space.

        Each row is a vector representation of an X feature. Features with similar
        embeddings tend to co-occur with similar sets of Y features. The Euclidean
        distance or cosine similarity between embeddings can be used to identify
        feature relatedness. The last column (+1; "bias") captures the baseline
        tendency of each X feature to associate with Y features overall.

        """
        return self._estimator.x_embeddings_

    @property
    def y_embeddings(self) -> TableLike:
        r"""Learned coordinates of conditioned features (Y) in latent space.

        See ``x_embeddings``. The first row is all zeros and represents the reference
        Y feature used for identifiability. The last column (+1; "bias") captures the
        baseline abundance of each Y feature.

        """
        return self._estimator.y_embeddings_

    @property
    def ranks(self) -> TableLike:
        r"""Log conditional probability matrix of co-occurrence of X and Y features.

        Entry (i, j) represents the log-odds of observing :math:`Y_j` given
        :math:`X_i`, relative to the row mean. Higher values indicate stronger
        positive associations. This matrix is row-centered (each row sums to zero)
        for identifiability.

        To obtain actual conditional probabilities, transform the matrix using
        :func:`~skbio.stats.composition.clr_inv`.

        """
        return self._estimator.ranks_

    @property
    def convergence(self) -> TableLike:
        r"""Loss (negative log-posterior) over iterations during training.

        Use this to diagnose training issues. The loss should generally decrease and
        stabilize. If the loss is still decreasing at the final iteration, consider
        increasing ``max_iter``. If the loss oscillates (Adam optimizer), try reducing
        ``learning_rate``.

        """
        return self._estimator.loss_curve_

    def predict(self, X: TableLike) -> TableLike:
        r"""Predict conditioned feature compositions given conditioning features.

        The expected conditioned (Y) feature distribution for each sample is computed
        by marginalizing over the given conditioning (X) feature composition and the
        learned conditional probabilities:

        .. math::

            P(Y) = \sum_i P(X_i)\,P(Y \mid X_i)

        Parameters
        ----------
        X : table_like of shape (n_samples, n_features_x)
            Feature abundance table of the conditioning (X) modality. Columns must
            match the features used during training.

        Returns
        -------
        table_like of shape (n_samples, n_features_y)
            Predicted feature proportions of the conditioned (Y) modality for each
            sample. Each row sums to 1.

        """
        return self._estimator.predict(X)

    def score(self, X: TableLike, Y: TableLike) -> float:
        r"""Compute Q-squared (coefficient of prediction) on held-out data.

        :math:`Q^2` measures predictive performance on test data, analogous to
        :math:`R^2` but for cross-validation. Values range from -inf to 1, where 1
        indicates perfect prediction and 0 indicates prediction no better than the
        mean.

        .. math::

            Q^2 = 1 - \frac{SS_{res}}{SS_{tot}}
                = 1 - \frac{\sum(y - \hat{y})^2}{\sum(y - \bar{y}_j)^2}

        where :math:`\bar{y}_j` is the per-target mean across samples.

        Parameters
        ----------
        X : table_like of shape (n_samples, n_features_x)
            Feature abundance table of the conditioning (X) modality. Columns must
            match the features used during training.
        Y : table_like of shape (n_samples, n_features_y)
            Feature abundance table of the conditioned (Y) modality.

        Returns
        -------
        q2 : float
            Q-squared score. Higher is better, with 1.0 being perfect prediction.

        """
        return self._estimator.score(X, Y)


class MMvec(SkbioObject):
    r"""MMvec estimator with scikit-learn-style interface.

    MMvec models conditional target distributions given feature compositions. In sklearn
    terms, the first modality (``X``) are "features" and the second modality (``y``,
    vector or matrix) are "targets". In the original MMvec work, features typically
    correspond to microbes and targets to metabolites.

    .. versionadded:: 0.7.3

    Attributes
    ----------
    x_embeddings_ : table_like of shape (n_features, n_components + 1)
        Learned feature embeddings (+1 is the bias term).
    y_embeddings_ : table_like of shape (n_outputs, n_components + 1)
        Learned target embeddings (+1 is the bias term).
    ranks_ : table_like of shape (n_features, n_outputs)
        Row-centered log conditional probabilities.
    loss_curve_ : table_like of shape (n_iterations,)
        Optimization losses over training updates.
    n_iter_ : int
        Number of optimization iterations completed.
    n_features_in_ : int
        Number of features in `X` seen during :meth:`fit`.
    n_outputs_ : int
        Number of features in `y` seen during :meth:`fit`.
    x_feature_ids_ : tuple or None
        Feature IDs of `X`, if available.
    y_feature_ids_ : tuple or None
        Feature IDs of `y`, if available.
    """

    def __init__(
        self,
        n_components: int = 3,
        optimizer: str = "lbfgs",
        max_iter: int = 1000,
        x_prior_mean: float = 0.0,
        x_prior_scale: float = 1.0,
        y_prior_mean: float = 0.0,
        y_prior_scale: float = 1.0,
        learning_rate: float = 1e-3,
        batch_size: int = 50,
        beta_1: float = 0.9,
        beta_2: float = 0.95,
        clipnorm: float = 10.0,
        batch_norm: str = "unbiased",
        seed: SeedLike | None = None,
        verbose: bool = False,
        output_format: str | None = None,
    ) -> None:
        self.n_components = n_components
        self.optimizer = optimizer
        self.max_iter = max_iter
        self.x_prior_mean = x_prior_mean
        self.x_prior_scale = x_prior_scale
        self.y_prior_mean = y_prior_mean
        self.y_prior_scale = y_prior_scale
        self.learning_rate = learning_rate
        self.batch_size = batch_size
        self.beta_1 = beta_1
        self.beta_2 = beta_2
        self.clipnorm = clipnorm
        self.batch_norm = batch_norm
        self.seed = seed
        self.verbose = verbose
        self.output_format = output_format

    def fit(self, X: TableLike, y: TableLike) -> MMvec:
        """Fit MMvec model.

        Parameters
        ----------
        X : table_like of shape (n_samples, n_features)
            Conditioning feature matrix (e.g., microbes).
        y : table_like of shape (n_samples,) or (n_samples, n_outputs)
            Conditioned feature matrix or vector (e.g., metabolites).

        Returns
        -------
        MMvec
            Fitted estimator.

        """
        # Parse tabular inputs
        X_arr, _, x_feature_ids = _ingest_table(X)
        y_arr, _, y_feature_ids = _ingest_table(y)

        # Check for sample count consistency
        if X_arr.shape[0] != y_arr.shape[0]:
            raise ValueError(
                f"X and y must have the same number of samples. "
                f"Got {X_arr.shape[0]} and {y_arr.shape[0]}."
            )

        # Check for all-zero columns / rows
        for mod, name in (X_arr, "X"), (y_arr, "y"):
            for i, axis in enumerate(("columns", "rows")):
                if np.any(mod.sum(axis=i) == 0):
                    raise ValueError(
                        f"{name} contains all-zero {axis}. Remove them before calling."
                    )

        # Check for positive prior scales
        if self.x_prior_scale <= 0:
            raise ValueError("x_prior_scale must be positive.")
        if self.y_prior_scale <= 0:
            raise ValueError("y_prior_scale must be positive.")

        # Validate optimizer
        optimizer = self.optimizer.lower()
        if optimizer not in ("lbfgs", "adam"):
            raise ValueError("Optimizer must be 'lbfgs' or 'adam'.")

        # Create RNG
        rng = get_rng(self.seed)

        n_features_x = X_arr.shape[1]
        n_features_y = y_arr.shape[1]

        # Initialize model. Note the change of terminology.
        model = _MMvecModel(
            n_features_x=n_features_x,
            n_features_y=n_features_y,
            n_components=self.n_components,
            x_prior_mean=self.x_prior_mean,
            x_prior_scale=self.x_prior_scale,
            y_prior_mean=self.y_prior_mean,
            y_prior_scale=self.y_prior_scale,
            rng=rng,
        )

        # Convert to sparse COO format
        X_coo = coo_array(X_arr)

        # Train model
        if optimizer == "lbfgs":
            losses = _train_lbfgs(model, X_coo, y_arr, self.max_iter, self.verbose)
        else:
            losses = _train_adam(
                model=model,
                X_coo=X_coo,
                Y=y_arr,
                rng=rng,
                max_iter=self.max_iter,
                learning_rate=self.learning_rate,
                batch_size=self.batch_size,
                beta_1=self.beta_1,
                beta_2=self.beta_2,
                clipnorm=self.clipnorm,
                batch_norm=self.batch_norm,
                verbose=self.verbose,
            )

        pc_cols = [f"PC{i + 1}" for i in range(model.n_components)] + ["bias"]
        x_embed = np.hstack([model.x_main, model.x_bias])
        y_embed = np.vstack(
            [
                np.zeros((1, model.n_components + 1)),
                np.hstack([model.y_main.T, model.y_bias.T]),
            ]
        )
        self.x_embeddings_ = _create_table(
            x_embed, columns=pc_cols, index=x_feature_ids, backend=self.output_format
        )
        self.y_embeddings_ = _create_table(
            y_embed, columns=pc_cols, index=y_feature_ids, backend=self.output_format
        )
        self.loss_curve_ = _create_table_1d(losses)

        # Keep the learned low-rank parameters available for lazy post-fit views.
        # `ranks_` and `probs` are both dense (n_features_x, n_features_y) matrices,
        # so we avoid materializing either one until a user actually needs it.
        self._x_main = model.x_main
        self._x_bias = model.x_bias
        self._y_main = model.y_main
        self._y_bias = model.y_bias
        self._ranks = None
        self._probs = None

        self.x_feature_ids_ = (
            tuple(x_feature_ids) if x_feature_ids is not None else None
        )
        self.y_feature_ids_ = (
            tuple(y_feature_ids) if y_feature_ids is not None else None
        )
        self.n_features_in_ = X_arr.shape[1]
        self.n_outputs_ = y_arr.shape[1]
        self.n_iter_ = len(losses)
        self.is_fitted_ = True

        return self

    def _check_is_fitted(self) -> None:
        if not getattr(self, "is_fitted_", False):
            raise ValueError("MMvec estimator is not fitted. Call fit(X, y) first.")

    def __str__(self) -> str:
        """Return string representation of MMvec."""
        self._check_is_fitted()
        n_features_x = self.x_embeddings_.shape[0]
        n_features_y = self.y_embeddings_.shape[0]
        n_components = self.x_embeddings_.shape[1] - 1  # exclude bias
        n_iterations = len(self.loss_curve_)
        return (
            f"MMvec\n"
            f"  Features: {n_features_x}\n"
            f"  Targets: {n_features_y}\n"
            f"  Components: {n_components}\n"
            f"  Iterations: {n_iterations}"
        )

    def get_params(self, deep: bool = True) -> dict[str, object]:
        """Get estimator parameters for sklearn compatibility."""
        return {
            "n_components": self.n_components,
            "optimizer": self.optimizer,
            "max_iter": self.max_iter,
            "learning_rate": self.learning_rate,
            "batch_size": self.batch_size,
            "x_prior_mean": self.x_prior_mean,
            "x_prior_scale": self.x_prior_scale,
            "y_prior_mean": self.y_prior_mean,
            "y_prior_scale": self.y_prior_scale,
            "beta_1": self.beta_1,
            "beta_2": self.beta_2,
            "clipnorm": self.clipnorm,
            "batch_norm": self.batch_norm,
            "seed": self.seed,
            "verbose": self.verbose,
            "output_format": self.output_format,
        }

    def set_params(self, **params: object) -> MMvec:
        """Set estimator parameters for sklearn compatibility."""
        for key, value in params.items():
            if not hasattr(self, key):
                raise ValueError(f"Invalid parameter '{key}' for MMvec estimator.")
            setattr(self, key, value)
        return self

    @property
    def ranks_(self) -> TableLike:
        """Row-centered log conditional probabilities.

        This dense matrix is cached lazily because it can be very large for
        high-dimensional multiomics data and is primarily used for inspection.
        """
        self._check_is_fitted()
        if self._ranks is None:
            ranks = self._compute_logits()
            ranks -= np.mean(ranks, axis=1, keepdims=True)
            self._ranks = _create_table(
                ranks,
                columns=self.y_feature_ids_,
                index=self.x_feature_ids_,
                backend=self.output_format,
            )
        return self._ranks

    def _compute_logits(self) -> np.ndarray:
        """Construct full Y logits from the learned low-rank parameters.

        The first Y feature is the reference category and is fixed at zero logits.
        Remaining columns are computed directly from the fitted X/Y embeddings and
        bias terms.
        """
        logits = np.empty((self.n_features_in_, self.n_outputs_))
        logits[:, 0] = 0.0
        logits_nr = logits[:, 1:]
        np.matmul(self._x_main, self._y_main, out=logits_nr)
        logits_nr += self._x_bias
        logits_nr += self._y_bias
        return logits

    def _get_probs(self) -> np.ndarray:
        """Lazily cache the dense conditional probability matrix P(Y | X)."""
        if self._probs is None:
            self._probs = softmax(self._compute_logits(), validate=False)
        return self._probs

    def _predict(self, X: np.ndarray) -> np.ndarray:
        """Predict Y proportions given X, without table I/O."""
        # Normalize abundances to proportions before applying P(Y | X).
        row_sums = np.sum(X, axis=1, keepdims=True)
        if np.any(row_sums == 0):
            raise ValueError(
                "X contains samples with all-zero counts. Remove these samples before "
                "calling predict."
            )
        X_props = X / row_sums
        return X_props @ self._get_probs()

    def predict(self, X: TableLike) -> TableLike:
        """Predict target distributions given feature compositions.

        Parameters
        ----------
        X : table_like of shape (n_samples, n_features)
            Feature abundance/count table. Columns must match the features used during
            training.

        Returns
        -------
        y_pred : table_like of shape (n_samples, n_outputs)
            Predicted target (y) proportions for each sample. Each row sums to 1.

        """
        self._check_is_fitted()
        X, sample_ids, _ = _ingest_table(X)
        y_pred = self._predict(X)
        return _create_table(
            y_pred,
            columns=self.y_feature_ids_,
            index=sample_ids,
            backend=self.output_format,
        )

    def score(self, X: TableLike, y: TableLike) -> float:
        r"""Compute Q-squared (coefficient of prediction) on held-out data.

        Parameters
        ----------
        X : table_like of shape (n_samples, n_features)
            Test feature table.
        y : table_like of shape (n_samples,) or (n_samples, n_outputs)
            Test target table (or vector).

        Returns
        -------
        q2 : float
            Q-squared score. Higher is better, with 1.0 being perfect prediction.

        """
        self._check_is_fitted()
        X, _, _ = _ingest_table(X)
        Y, _, _ = _ingest_table(y)

        # Predict directly into an array to avoid a table round-trip.
        Y_pred = self._predict(X)

        # Normalize abundances to proportions
        row_sums = np.sum(Y, axis=1, keepdims=True)
        if np.any(row_sums == 0):
            raise ValueError(
                "y contains samples with all-zero counts. "
                "Remove these samples before calling score."
            )
        Y_true = Y / row_sums

        # Q^2 = 1 - SS_res / SS_tot
        ss_res = np.sum((Y_true - Y_pred) ** 2)
        ss_tot = np.sum((Y_true - Y_true.mean(axis=0)) ** 2)

        # All true values are the same - return 0 to avoid division by zero
        # NOTE: Mathematically, this should happen when every feature has a uniform
        # value across samples. But due to floating-point arithmetic issues, this
        # won't happen.
        if ss_tot == 0:  # pragma: no cover
            return 0.0

        return 1.0 - ss_res / ss_tot


class _MMvecModel:
    """Internal model class for MMvec optimization."""

    def __init__(
        self,
        n_features_x: int,
        n_features_y: int,
        n_components: int,
        x_prior_mean: float,
        x_prior_scale: float,
        y_prior_mean: float,
        y_prior_scale: float,
        rng: np.random.Generator,
    ) -> None:
        """Initialize MMvec model parameters.

        Parameters
        ----------
        n_features_x : int
            Number of features in the conditioning modality (X).
        n_features_y : int
            Number of features in the conditioned modality (Y).
        n_components : int
            Latent dimensionality (p).
        x_prior_mean : float
            Mean of Gaussian prior on X-side embeddings and biases.
        x_prior_scale : float
            Scale of Gaussian prior on X-side embeddings and biases.
        y_prior_mean : float
            Mean of Gaussian prior on Y-side embeddings and biases.
        y_prior_scale : float
            Scale of Gaussian prior on Y-side embeddings and biases.
        rng : numpy.random.Generator
            Random number generator.

        """
        self.n_features_x = n_features_x
        self.n_features_y = n_features_y
        self.n_components = n_components
        self.x_prior_mean = x_prior_mean
        self.x_prior_scale = x_prior_scale
        self.y_prior_mean = y_prior_mean
        self.y_prior_scale = y_prior_scale
        # Cache squared scales used repeatedly in objective/gradient evaluation.
        self._x_prior_scale_sq = x_prior_scale**2
        self._y_prior_scale_sq = y_prior_scale**2

        # Initialize parameters with random normal
        self.x_main = rng.standard_normal((n_features_x, n_components))
        self.x_bias = rng.standard_normal((n_features_x, 1))
        self.y_main = rng.standard_normal((n_components, n_features_y - 1))
        self.y_bias = rng.standard_normal((1, n_features_y - 1))

    def loss_and_gradients(
        self,
        X_coo: coo_array,
        Y: np.ndarray,
        size: int,
        norm: float,
        weights: np.ndarray,
        rng: np.random.Generator,
    ) -> tuple[
        float,
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    ]:
        """Compute loss and gradients for a mini-batch.

        Parameters
        ----------
        X_coo : coo_array of shape (n_samples, n_features_x)
            Microbe counts in COO format.
        Y : np.ndarray of shape (n_samples, n_features_y)
            Metabolite counts.
        size : int
            Mini-batch size.
        norm : float
            Batch normalization factor.
        weights : np.ndarray of shape (n_samples,)
            Sample weights for weighted batch sampling.
        rng : numpy.random.Generator
            Random number generator for batch sampling.

        Returns
        -------
        loss : float
            Negative log posterior.
        grads : tuple of np.ndarray
            Gradients in the fixed order
            ``(dx_main, dx_bias, dy_main, dy_bias)``.
        """
        batch_idx = rng.choice(len(X_coo.data), size=size, replace=True, p=weights)

        sample_ids = X_coo.row[batch_idx]
        X_ids = X_coo.col[batch_idx]

        # Build the non-reference logits directly for the sampled X features.
        Y_batch = Y[sample_ids, :]  # (B, d2)
        x_main_batch = self.x_main[X_ids, :]  # (B, p)
        x_bias_batch = self.x_bias[X_ids, :]  # (B, 1)
        logits_nr = x_main_batch @ self.y_main  # (B, d2-1)
        logits_nr += x_bias_batch
        logits_nr += self.y_bias

        # Stable log(1 + sum(exp(logits_nr))) with the reference Y category
        # kept implicit as a zero logit.
        row_max = np.max(logits_nr, axis=1, keepdims=True)
        np.maximum(row_max, 0.0, out=row_max)
        exp_shifted = np.exp(logits_nr - row_max)
        log_norm = row_max + np.log(
            np.exp(-row_max) + exp_shifted.sum(axis=1, keepdims=True)
        )

        # Non-reference probabilities and grouped residuals for the sampled batch.
        probs_nr = exp_shifted / np.exp(log_norm - row_max)
        totals = np.sum(Y_batch, axis=1, keepdims=True)
        delta = Y_batch[:, 1:] - totals * probs_nr

        # The reference category has logit 0, so only non-reference logits appear in
        # the inner-product term.
        loglik = np.sum(Y_batch[:, 1:] * logits_nr) - np.sum(totals * log_norm)

        # Gradient w.r.t. Y-side embedding parameters
        dy_main = -norm * (x_main_batch.T @ delta)
        dy_main += (self.y_main - self.y_prior_mean) / (self.y_prior_scale**2)

        # Y-side bias gradient is the row-wise sum over non-reference residuals.
        dy_bias = -norm * delta.sum(axis=0, keepdims=True)
        dy_bias += (self.y_bias - self.y_prior_mean) / (self.y_prior_scale**2)

        # Project residuals back through the Y embeddings for the X-side gradients.
        dx_main_batch = delta @ self.y_main.T  # (B, p)
        dx_bias_batch = delta.sum(axis=1)  # (B,)

        dx_main = np.zeros_like(self.x_main)
        dx_bias = np.zeros_like(self.x_bias)

        # Scatter-add the sampled X-side contributions back to full parameter arrays.
        np.add.at(dx_main, X_ids, -norm * dx_main_batch)
        np.add.at(dx_bias[:, 0], X_ids, -norm * dx_bias_batch)

        # Add prior gradients
        dx_main += (self.x_main - self.x_prior_mean) / (self.x_prior_scale**2)
        dx_bias += (self.x_bias - self.x_prior_mean) / (self.x_prior_scale**2)

        # Compute total loss (negative log posterior)
        prior_loss = 0.0
        prior_loss += (
            0.5
            * np.sum((self.x_main - self.x_prior_mean) ** 2)
            / self._x_prior_scale_sq
        )
        prior_loss += (
            0.5
            * np.sum((self.x_bias - self.x_prior_mean) ** 2)
            / self._x_prior_scale_sq
        )
        prior_loss += (
            0.5
            * np.sum((self.y_main - self.y_prior_mean) ** 2)
            / self._y_prior_scale_sq
        )
        prior_loss += (
            0.5
            * np.sum((self.y_bias - self.y_prior_mean) ** 2)
            / self._y_prior_scale_sq
        )

        loss = -norm * loglik + prior_loss

        grads = (dx_main, dx_bias, dy_main, dy_bias)

        return loss, grads

    def pack_params(self) -> np.ndarray:
        """Flatten all parameters into a single vector.

        Returns
        -------
        theta : np.ndarray
            Flattened parameters [x_main, x_bias, y_main, y_bias].
        """
        return np.concatenate(
            [
                self.x_main.ravel(),
                self.x_bias.ravel(),
                self.y_main.ravel(),
                self.y_bias.ravel(),
            ]
        )

    def unpack_params(self, theta: np.ndarray) -> None:
        """Restore parameter matrices from flat vector.

        Parameters
        ----------
        theta : np.ndarray
            Flattened parameters.
        """
        d1 = self.n_features_x
        d2 = self.n_features_y
        p = self.n_components

        idx = 0
        self.x_main = theta[idx : idx + d1 * p].reshape(d1, p)
        idx += d1 * p
        self.x_bias = theta[idx : idx + d1].reshape(d1, 1)
        idx += d1
        self.y_main = theta[idx : idx + p * (d2 - 1)].reshape(p, d2 - 1)
        idx += p * (d2 - 1)
        self.y_bias = theta[idx : idx + (d2 - 1)].reshape(1, d2 - 1)

    def full_batch_loss_and_gradient(
        self,
        y_sums: np.ndarray,
        n_sums: np.ndarray,
        logits: np.ndarray,
        resids: np.ndarray,
        row_max: np.ndarray,
        log_norm: np.ndarray,
        row_scale: np.ndarray,
    ) -> tuple[float, np.ndarray]:
        """Compute full-batch loss and gradient for L-BFGS.

        Parameters
        ----------
        y_sums : np.ndarray of shape (n_features_x, n_features_y)
            Weighted Y counts aggregated by X feature index.
            Each row stores the total Y counts contributed by all nonzero X entries
            for one X feature.
        n_sums : np.ndarray of shape (n_features_x,)
            Weighted Y sample totals aggregated by X feature index.
            Each entry stores the total Y count paired with one X feature.
        logits : np.ndarray of shape (n_features_x, n_features_y - 1)
            Reusable workspace for non-reference logits.
        resids : np.ndarray of shape (n_features_x, n_features_y - 1)
            Reusable workspace for grouped residuals and non-reference probabilities.
        row_max : np.ndarray of shape (n_features_x,)
            Reusable workspace for row-wise maxima.
        log_norm : np.ndarray of shape (n_features_x,)
            Reusable workspace for row-wise log normalizers.
        row_scale : np.ndarray of shape (n_features_x,)
            Reusable workspace for row-wise normalizing factors.

        Returns
        -------
        loss : float
            Negative log-posterior.
        grad : np.ndarray
            Flattened gradient vector.

        Notes
        -----
        This routine evaluates the full negative log-posterior for the L-BFGS solver.

        It first groups the data by X feature using two sufficient statistics:

        - y_sums: grouped Y counts
        - n_sums: grouped Y totals

        With those grouped arrays in hand, the dense part of the objective is computed
        from the non-reference logits:

            x_main @ y_main + x_bias + y_bias

        without materializing the reference Y column. The data term is then:

            sum(n_sums * log_norm) - sum(y_sums[:, 1:] * logits)

        followed by the Gaussian prior terms for the X-side and Y-side parameters.
        The same grouped quantities are reused to form the gradient.

        """
        # Symbols used in comments below:
        # d1 : number of X features
        # d2 : number of Y features
        # p : number of latent dimensions
        # nnz : number of non-zero entries in X

        # Compute the non-reference logits directly:
        #   logits = x_main @ y_main + x_bias + y_bias
        # This is algebraically the same as the old augmented-matrix product, but
        # avoids the reference column (1).
        np.matmul(self.x_main, self.y_main, out=logits)
        logits += self.x_bias
        logits += self.y_bias

        # Stable softmax for all X features (rows).
        # Compute log_norm = log(1 + \sum exp(logits)) in a numerically stable way that
        # resembles SciPy's logsumexp, but without allocating the reference column (0):
        #   log_norm = m + log(exp(-m) + \sum exp(logits - m))
        # In which m = max(0, max(logits))
        # And exp(logits - m) is stored in resids for reuse.
        np.max(logits, axis=1, out=row_max)
        np.maximum(row_max, 0.0, out=row_max)
        np.subtract(logits, row_max[:, None], out=resids)
        np.exp(resids, out=resids)
        np.sum(resids, axis=1, out=log_norm)
        np.negative(row_max, out=row_scale)
        np.exp(row_scale, out=row_scale)
        log_norm += row_scale
        np.log(log_norm, out=log_norm)
        log_norm += row_max

        # Compute the grouped negative log-likelihood (part of the loss).
        # For each X feature, we combine:
        # - the weighted observed Y totals (`y_sums`)
        # - the weighted Y sample totals (`n_sums`)
        # - the corresponding log-normalizer (`log_norm`)
        # In compact form:
        #   loss_data = sum(n_sums * log_norm) - sum(y_sums[:, 1:] * logits)
        # The reference category has logit 0, so only non-reference logits appear
        # in the inner-product term.
        loss = np.dot(n_sums, log_norm) - np.einsum("ij,ij->", y_sums[:, 1:], logits)

        # Recover non-reference probabilities from the shifted exponentials.
        # At this point resids holds exp(logits - row_max). Divide by
        # exp(log_norm - row_max) to recover the non-reference probabilities.
        np.subtract(log_norm, row_max, out=row_scale)
        np.exp(row_scale, out=row_scale)
        resids /= row_scale[:, None]

        # Turn probabilities into grouped residuals: observed weighted totals
        # minus expected weighted totals for each non-reference category.
        resids *= n_sums[:, None]
        resids *= -1.0
        resids += y_sums[:, 1:]

        # Gradient for the Y-side embedding matrix.
        dy_main = -self.x_main.T @ resids  # (p, d2-1)

        # Gradient for Y-side biases: sum residuals across X features.
        dy_bias = -resids.sum(axis=0, keepdims=True)  # (1, d2-1)

        # Gradient for the X-side embedding matrix.
        dx_main = -(resids @ self.y_main.T)  # (d1, p)

        # Gradient for X-side biases: sum residuals across Y categories.
        dx_bias = -resids.sum(axis=1, keepdims=True)  # (d1, 1)

        # Add prior terms (L2 regularization)
        # These are exactly the Gaussian-prior contributions to the negative
        # log-posterior and its gradient.
        x_main_diff = self.x_main - self.x_prior_mean
        x_bias_diff = self.x_bias - self.x_prior_mean
        y_main_diff = self.y_main - self.y_prior_mean
        y_bias_diff = self.y_bias - self.y_prior_mean

        loss += 0.5 * np.sum(x_main_diff**2) / self._x_prior_scale_sq
        loss += 0.5 * np.sum(x_bias_diff**2) / self._x_prior_scale_sq
        loss += 0.5 * np.sum(y_main_diff**2) / self._y_prior_scale_sq
        loss += 0.5 * np.sum(y_bias_diff**2) / self._y_prior_scale_sq

        dx_main += x_main_diff / self._x_prior_scale_sq
        dx_bias += x_bias_diff / self._x_prior_scale_sq
        dy_main += y_main_diff / self._y_prior_scale_sq
        dy_bias += y_bias_diff / self._y_prior_scale_sq

        # Pack gradients
        grad = np.concatenate(
            [
                dx_main.ravel(),
                dx_bias.ravel(),
                dy_main.ravel(),
                dy_bias.ravel(),
            ]
        )

        return loss, grad


def _train_lbfgs(
    model: _MMvecModel,
    X_coo: coo_array,
    Y: np.ndarray,
    max_iter: int,
    verbose: bool,
) -> list[float]:
    """Train MMvec model using L-BFGS-B optimization.

    Parameters
    ----------
    model : _MMvecModel
        Initialized model to train.
    X_coo : coo_array of shape (n_samples, n_features_x)
        X feature counts in sparse COO format.
    Y : np.ndarray of shape (n_samples, n_features_y)
        Y feature counts.
    max_iter : int
        Maximum number of L-BFGS iterations.
    verbose : bool
        Print training progress.

    Returns
    -------
    losses : list of float
        Loss per iteration.

    """
    from scipy.optimize import minimize

    rows, cols, data = X_coo.row, X_coo.col, X_coo.data
    d1 = model.n_features_x
    d2 = model.n_features_y

    # Infer the data type for computation. Currently it promotes to at least float64.
    # But we should consider float32 support in the future.
    dtype = np.result_type(Y.dtype, data.dtype, np.float64)

    # Precompute X feature-grouped sufficient statistics reused by every iteration.
    # y_sums groups Y counts by X feature. Each row is the total Y abundance associated
    # with one X feature across all nonzero X entries.
    prods_ = data[:, None] * Y[rows, :]  # (nnz, d2)
    y_sums = np.empty((d1, d2), dtype=dtype)
    for j in range(d2):
        y_sums[:, j] = np.bincount(cols, weights=prods_[:, j], minlength=d1)

    # n_sums stores the corresponding grouped Y totals. It is the row sum partner to
    # y_sums and tells us how much total Y abundance is attached to each X feature.
    sums_ = Y.sum(axis=1)[rows]  # (nnz,)
    n_sums = np.bincount(cols, weights=data * sums_, minlength=d1).astype(
        dtype, copy=False
    )

    # Reusable array workspaces for each iteration.
    # Logits for each X feature, excluding the reference (first column, all zeros).
    logits = np.empty((d1, d2 - 1), dtype=dtype)

    # Exponential-shifted logits reused later as grouped residuals.
    resids = np.empty((d1, d2 - 1), dtype=dtype)

    # Row-wise maximum over logits for numerically stable normalization.
    row_max = np.empty(d1, dtype=dtype)

    # Row-wise log normalizer, i.e., log(1 + sum(exp(logits))).
    log_norm = np.empty(d1, dtype=dtype)

    # Row-wise scaling factors derived from row_max/log_norm.
    row_scale = np.empty(d1, dtype=dtype)

    losses = []
    it = 0

    def func(theta: np.ndarray) -> tuple[float, np.ndarray]:
        nonlocal it
        model.unpack_params(theta)
        loss, grad = model.full_batch_loss_and_gradient(
            y_sums, n_sums, logits, resids, row_max, log_norm, row_scale
        )
        it += 1
        losses.append(loss)

        if verbose and it % 10 == 0:
            print(f"Iteration {it}, Loss: {loss:.4f}")

        return loss, grad

    # Initial parameters
    theta0 = model.pack_params()

    # Run L-BFGS-B
    result = minimize(
        func,
        theta0,
        method="L-BFGS-B",
        jac=True,
        options={
            "maxiter": max_iter,
            "ftol": 1e-9,
            "gtol": 1e-5,
        },
    )

    # Unpack final parameters
    model.unpack_params(result.x)

    if verbose:
        print(
            f"L-BFGS-B converged: {result.success}, "
            f"iterations: {result.nit}, message: {result.message}"
        )

    return losses


def _train_adam(
    model: _MMvecModel,
    X_coo: coo_array,
    Y: np.ndarray,
    rng: np.random.Generator,
    max_iter: int,
    learning_rate: float,
    batch_size: int,
    beta_1: float,
    beta_2: float,
    clipnorm: float,
    batch_norm: str,
    verbose: bool,
) -> list[float]:
    """Train MMvec model using Adam optimization.

    Parameters
    ----------
    model : _MMvecModel
        Initialized model to train.
    X_coo : coo_array of shape (n_samples, n_features_x)
        Microbe counts in sparse COO format.
    Y : np.ndarray of shape (n_samples, n_features_y)
        Metabolite counts.
    rng : np.random.Generator
        Random number generator.
    max_iter : int
        Number of epochs.
    learning_rate : float
        Adam learning rate.
    batch_size : int
        Mini-batch size.
    beta_1 : float
        Adam exponential decay rate for first moment.
    beta_2 : float
        Adam exponential decay rate for second moment.
    clipnorm : float
        Gradient clipping threshold.
    batch_norm : {'legacy', 'unbiased'}
        Batch normalization mode.
    verbose : bool
        Print training progress.

    Returns
    -------
    losses : list of dict
        Training metrics per iteration.

    """
    losses = []

    # Initialize Adam moments
    x_main_m = np.zeros_like(model.x_main)
    x_main_v = np.zeros_like(model.x_main)
    x_bias_m = np.zeros_like(model.x_bias)
    x_bias_v = np.zeros_like(model.x_bias)
    y_main_m = np.zeros_like(model.y_main)
    y_main_v = np.zeros_like(model.y_main)
    y_bias_m = np.zeros_like(model.y_bias)
    y_bias_v = np.zeros_like(model.y_bias)

    # Compute number of iterations per epoch
    data = X_coo.data
    nnz = len(data)
    iter_per_epoch = max(1, nnz // batch_size)

    # Compute normalization factor
    if batch_norm == "legacy":
        norm = X_coo.shape[0] / batch_size
    else:  # unbiased
        norm = data.sum() / batch_size

    # Sample batch weighted by abundance
    weights = data / data.sum()

    params = (learning_rate, beta_1, beta_2)
    it = 0
    for epoch in range(max_iter):
        for _ in range(iter_per_epoch):
            it += 1

            # Compute loss and gradients
            loss, grads = model.loss_and_gradients(
                X_coo, Y, batch_size, norm, weights, rng
            )

            # Gradient clipping
            dx_main, dx_bias, dy_main, dy_bias = _clip_gradients(grads, clipnorm)

            # Adam updates
            bc1 = 1 - beta_1**it
            bc2 = 1 - beta_2**it

            model.x_main, x_main_m, x_main_v = _adam_update(
                model.x_main, dx_main, x_main_m, x_main_v, bc1, bc2, *params
            )
            model.x_bias, x_bias_m, x_bias_v = _adam_update(
                model.x_bias, dx_bias, x_bias_m, x_bias_v, bc1, bc2, *params
            )
            model.y_main, y_main_m, y_main_v = _adam_update(
                model.y_main, dy_main, y_main_m, y_main_v, bc1, bc2, *params
            )
            model.y_bias, y_bias_m, y_bias_v = _adam_update(
                model.y_bias, dy_bias, y_bias_m, y_bias_v, bc1, bc2, *params
            )

            losses.append(loss)

        if verbose:
            print(f"Epoch {epoch + 1}/{max_iter}, Loss: {loss:.4f}")

    return losses


def _clip_gradients(
    grads: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    clipnorm: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Apply global norm gradient clipping.

    Parameters
    ----------
    grads : tuple of np.ndarray
        Gradients in the fixed order
        ``(dx_main, dx_bias, dy_main, dy_bias)``.
    clipnorm : float
        Maximum gradient norm.

    Returns
    -------
    grads : tuple of np.ndarray
        Clipped gradients in the same fixed order.
    """
    global_norm = np.sqrt(sum(np.sum(g**2) for g in grads))

    if global_norm > clipnorm:
        scale = clipnorm / global_norm
        grads = tuple(g * scale for g in grads)

    return grads


def _adam_update(
    param: np.ndarray,
    grad: np.ndarray,
    m: np.ndarray,
    v: np.ndarray,
    bc1: float,
    bc2: float,
    lr: float,
    beta_1: float,
    beta_2: float,
    eps: float = 1e-8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Perform single Adam update step.

    Parameters
    ----------
    param : np.ndarray
        Current parameter values.
    grad : np.ndarray
        Gradient.
    m : np.ndarray
        First moment estimate.
    v : np.ndarray
        Second moment estimate.
    bc1 : float
        First-moment bias-correction denominator.
    bc2 : float
        Second-moment bias-correction denominator.
    lr : float
        Learning rate.
    beta_1 : float
        Exponential decay rate for first moment.
    beta_2 : float
        Exponential decay rate for second moment.
    eps : float
        Small constant for numerical stability.

    Returns
    -------
    param : np.ndarray
        Updated parameter values.
    m : np.ndarray
        Updated first moment.
    v : np.ndarray
        Updated second moment.
    """
    m = beta_1 * m + (1 - beta_1) * grad
    v = beta_2 * v + (1 - beta_2) * grad**2

    # Bias correction
    m_hat = m / bc1
    v_hat = v / bc2
    param = param - lr * m_hat / (np.sqrt(v_hat) + eps)

    return param, m, v


def random_multimodal(
    n_features_x: int = 20,
    n_features_y: int = 100,
    n_samples: int = 100,
    n_components: int = 3,
    grad_low: float = -1,
    grad_high: float = 1,
    x_total: int = 10,
    y_total: int = 100,
    coef_mean: float = 0,
    coef_std: float = 2,
    x_noise_std: float = 0.1,
    x_mean: float = 0,
    x_std: float = 1,
    y_mean: float = 0,
    y_std: float = 1,
    seed: SeedLike | None = 0,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Generate synthetic co-occurrence data of two modalities.

    This helper function is currently kept private because it is specific to MMvec
    tests. In the future, it may be promoted to a public and generalized simulation
    utility if broader use cases emerge.

    Parameters
    ----------
    n_features_x : int
        Number of conditioning (X) features to simulate.
    n_features_y : int
        Number of conditioned (Y) features to simulate.
    n_samples : int
        Number of samples to generate.
    n_components : int
        Number of latent dimensions for the embeddings.
    grad_low : float
        Lower bound of the sample gradient.
    grad_high : float
        Upper bound of the sample gradient.
    x_total : int
        Total counts for X features per sample.
    y_total : int
        Total counts for Y features per sample.
    coef_mean : float
        Mean of regression coefficient distribution.
    coef_std : float
        Standard deviation of regression coefficient distribution.
    x_noise_std : float
        Standard deviation of noise in X composition.
    x_mean : float
        Mean of X embedding distribution.
    x_std : float
        Standard deviation of X embedding distribution.
    y_mean : float
        Mean of Y embedding distribution.
    y_std : float
        Standard deviation of Y embedding distribution.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

    Returns
    -------
    x_counts : pd.DataFrame of shape (n_samples, n_features_x)
        Count table of X abundances.
    y_counts : pd.DataFrame of shape (n_samples, n_features_y)
        Count table of Y abundances.
    design : np.ndarray of shape (n_samples, 2)
        Design matrix with intercept and gradient.
    coefs : np.ndarray of shape (2, n_features_x)
        Regression coefficients.
    x_main : np.ndarray of shape (n_features_x, n_dimensions)
        True X embedding matrix.
    x_bias : np.ndarray of shape (n_features_x, 1)
        True X bias vector.
    y_main : np.ndarray of shape (n_dimensions, n_features_y - 1)
        True Y embedding matrix.
    y_bias : np.ndarray of shape (1, n_features_y - 1)
        True Y bias vector.

    """
    rng = get_rng(seed)

    # Regression coefficients for gradient model
    coefs = rng.normal(coef_mean, coef_std, size=(2, n_features_x))

    # Design matrix: intercept + gradient
    design = np.vstack(
        (np.ones(n_samples), np.linspace(grad_low, grad_high, n_samples))
    ).T

    # Generate X compositions from ILR with noise
    X = ilr_inv(
        rng.multivariate_normal(
            mean=np.zeros(n_features_x - 1),
            cov=np.diag([x_noise_std] * (n_features_x - 1)),
            size=n_samples,
        )
    )

    # Generate latent embeddings
    x_main = rng.normal(x_mean, x_std, size=(n_features_x, n_components))
    y_main = rng.normal(y_mean, y_std, size=(n_components, n_features_y - 1))

    x_bias = rng.normal(x_mean, x_std, size=(n_features_x, 1))
    y_bias = rng.normal(y_mean, y_std, size=(1, n_features_y - 1))

    # Augmented matrices for computing probabilities
    X_ = np.hstack((np.ones((n_features_x, 1)), x_bias, x_main))
    Y_ = np.vstack((y_bias, np.ones((1, n_features_y - 1)), y_main))

    # Compute conditional probabilities P(Y|X)
    phi = np.hstack((np.zeros((n_features_x, 1)), X_ @ Y_))
    probs = softmax(phi, validate=False)

    # Generate count data
    x_counts = np.zeros((n_samples, n_features_x))
    y_counts = np.zeros((n_samples, n_features_y))

    n1 = x_total
    n2 = y_total // x_total

    for n in range(n_samples):
        # Draw X counts
        x = rng.multinomial(n1, X[n, :])
        # For each x, draw Y conditional on that x
        for i in range(n_features_x):
            y = rng.multinomial(x[i] * n2, probs[i, :])
            y_counts[n, :] += y
        x_counts[n, :] += x

    # Create DataFrames with meaningful IDs
    feature_ids_x = [f"x_feature_{d}" for d in range(x_counts.shape[1])]
    feature_ids_y = [f"y_feature_{d}" for d in range(y_counts.shape[1])]
    sample_ids = [f"sample_{d}" for d in range(y_counts.shape[0])]

    x_counts = pd.DataFrame(x_counts, index=sample_ids, columns=feature_ids_x)
    y_counts = pd.DataFrame(y_counts, index=sample_ids, columns=feature_ids_y)
    return x_counts, y_counts, design, coefs, x_main, x_bias, y_main, y_bias
