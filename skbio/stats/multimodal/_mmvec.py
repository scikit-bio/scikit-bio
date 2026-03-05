# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""MMvec - Microbe-Metabolite Vectors.

This module implements MMvec for learning joint embeddings of microbes
and metabolites from their co-occurrence patterns.
"""

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.sparse import coo_matrix, issparse

from skbio._base import SkbioObject
from skbio.stats.composition import clr_inv as softmax


def _multinomial_loglik_and_grad(logits, y):
    """Compute multinomial log-likelihood and gradient w.r.t. logits.

    Parameters
    ----------
    logits : np.ndarray of shape (batch_size, n_categories)
        Log-odds for each category.
    y : np.ndarray of shape (batch_size, n_categories)
        Observed counts for each category.

    Returns
    -------
    loglik : float
        Sum of log-likelihoods across batch.
    grad : np.ndarray of shape (batch_size, n_categories)
        Gradient of log-likelihood w.r.t. logits.

    Notes
    -----
    The log-likelihood (ignoring multinomial coefficient) is:
        sum_j y_j * log(pi_j) = sum_j y_j * (eta_j - logsumexp(eta))

    The gradient is:
        d/d eta_j = y_j - N * pi_j

    where N = sum(y) and pi = softmax(eta).
    """
    N = y.sum(axis=1, keepdims=True)
    log_norm = logsumexp(logits, axis=1, keepdims=True)
    loglik = np.sum(y * logits) - np.sum(N * log_norm)

    # Gradient: y - N * softmax(logits)
    probs = np.exp(logits - log_norm)
    grad = y - N * probs

    return loglik, grad


class _MMvecModel:
    """Internal model class for MMvec optimization."""

    def __init__(
        self,
        n_microbes,
        n_metabolites,
        n_components=3,
        u_prior_mean=0.0,
        u_prior_scale=1.0,
        v_prior_mean=0.0,
        v_prior_scale=1.0,
        rng=None,
    ):
        """Initialize MMvec model parameters.

        Parameters
        ----------
        n_microbes : int
            Number of microbes (d1).
        n_metabolites : int
            Number of metabolites (d2).
        n_components : int
            Latent dimensionality (p).
        u_prior_mean : float
            Mean of Gaussian prior on U.
        u_prior_scale : float
            Scale of Gaussian prior on U.
        v_prior_mean : float
            Mean of Gaussian prior on V.
        v_prior_scale : float
            Scale of Gaussian prior on V.
        rng : numpy.random.Generator
            Random number generator.
        """
        self.n_microbes = n_microbes
        self.n_metabolites = n_metabolites
        self.n_components = n_components
        self.u_prior_mean = u_prior_mean
        self.u_prior_scale = u_prior_scale
        self.v_prior_mean = v_prior_mean
        self.v_prior_scale = v_prior_scale

        if rng is None:
            rng = np.random.default_rng()

        # Initialize parameters with random normal
        self.U = rng.standard_normal((n_microbes, n_components))
        self.b_U = rng.standard_normal((n_microbes, 1))
        self.V = rng.standard_normal((n_components, n_metabolites - 1))
        self.b_V = rng.standard_normal((1, n_metabolites - 1))

    def _build_augmented_matrices(self):
        """Build augmented U and V matrices for forward pass.

        Returns
        -------
        U_aug : np.ndarray of shape (n_microbes, n_components + 2)
            [1 | b_U | U]
        V_aug : np.ndarray of shape (n_components + 2, n_metabolites - 1)
            [b_V; 1; V]
        """
        d1 = self.n_microbes
        p = self.n_components

        U_aug = np.hstack([np.ones((d1, 1)), self.b_U, self.U])
        V_aug = np.vstack(
            [self.b_V, np.ones((1, self.n_metabolites - 1)), self.V]
        )
        return U_aug, V_aug

    def forward(self, microbe_indices):
        """Compute logits for given microbe indices.

        Parameters
        ----------
        microbe_indices : np.ndarray of shape (batch_size,)
            Indices of microbes in the batch.

        Returns
        -------
        logits : np.ndarray of shape (batch_size, n_metabolites)
            Logits with first column (reference) set to 0.
        """
        U_aug, V_aug = self._build_augmented_matrices()
        u_batch = U_aug[microbe_indices, :]  # (B, p+2)
        logits_nonref = u_batch @ V_aug  # (B, d2-1)
        logits = np.hstack(
            [np.zeros((len(microbe_indices), 1)), logits_nonref]
        )
        return logits

    def loss_and_gradients(
        self,
        X,
        Y,
        batch_size=50,
        batch_normalization="legacy",
        rng=None,
    ):
        """Compute loss and gradients for a mini-batch.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_microbes)
            Microbe counts.
        Y : array-like of shape (n_samples, n_metabolites)
            Metabolite counts.
        batch_size : int
            Mini-batch size.
        batch_normalization : {'legacy', 'unbiased'}
            Batch normalization mode.
        rng : numpy.random.Generator, optional
            Random number generator for batch sampling.

        Returns
        -------
        loss : float
            Negative log posterior.
        grads : dict
            Gradients for U, b_U, V, b_V.
        """
        # Convert to sparse COO format if needed
        if issparse(X):
            X_coo = X.tocoo()
        else:
            X_coo = coo_matrix(X)

        n_samples = X.shape[0]
        total_count = X_coo.data.sum()

        # Compute normalization factor
        if batch_normalization == "legacy":
            norm = n_samples / batch_size
        else:  # unbiased
            norm = total_count / batch_size

        # Sample batch weighted by abundance
        weights = X_coo.data.astype(np.float64)
        weights /= weights.sum()

        if rng is None:
            rng = np.random.default_rng()
        batch_idx = rng.choice(
            len(X_coo.data), size=batch_size, replace=True, p=weights
        )

        sample_ids = X_coo.row[batch_idx]
        microbe_ids = X_coo.col[batch_idx]

        # Forward pass
        Y_arr = np.asarray(Y)
        Y_batch = Y_arr[sample_ids, :]  # (B, d2)
        logits = self.forward(microbe_ids)  # (B, d2)

        # Compute likelihood and gradient w.r.t. logits
        loglik, delta_full = _multinomial_loglik_and_grad(logits, Y_batch)
        delta = delta_full[:, 1:]  # Exclude reference category (B, d2-1)

        # Build augmented matrices
        U_aug, V_aug = self._build_augmented_matrices()
        u_batch = U_aug[microbe_ids, :]  # (B, p+2)

        # Gradient w.r.t. V parameters
        # dV = -norm * (u_batch[:, 2:].T @ delta) + V / sigma^2
        dV = -norm * (u_batch[:, 2:].T @ delta)
        dV += (self.V - self.v_prior_mean) / (self.v_prior_scale**2)

        # db_V = -norm * delta.sum(axis=0, keepdims=True) + b_V / sigma^2
        db_V = -norm * delta.sum(axis=0, keepdims=True)
        db_V += (self.b_V - self.v_prior_mean) / (self.v_prior_scale**2)

        # Gradient w.r.t. U parameters (scatter-add)
        du_batch = delta @ V_aug.T  # (B, p+2)

        dU = np.zeros_like(self.U)
        db_U = np.zeros_like(self.b_U)

        # Scatter-add gradients
        for b in range(batch_size):
            m = microbe_ids[b]
            dU[m, :] -= norm * du_batch[b, 2:]
            db_U[m, 0] -= norm * du_batch[b, 1]

        # Add prior gradients
        dU += (self.U - self.u_prior_mean) / (self.u_prior_scale**2)
        db_U += (self.b_U - self.u_prior_mean) / (self.u_prior_scale**2)

        # Compute total loss (negative log posterior)
        prior_loss = 0.0
        prior_loss += 0.5 * np.sum((self.U - self.u_prior_mean) ** 2) / (
            self.u_prior_scale**2
        )
        prior_loss += 0.5 * np.sum((self.b_U - self.u_prior_mean) ** 2) / (
            self.u_prior_scale**2
        )
        prior_loss += 0.5 * np.sum((self.V - self.v_prior_mean) ** 2) / (
            self.v_prior_scale**2
        )
        prior_loss += 0.5 * np.sum((self.b_V - self.v_prior_mean) ** 2) / (
            self.v_prior_scale**2
        )

        loss = -norm * loglik + prior_loss

        grads = {"U": dU, "b_U": db_U, "V": dV, "b_V": db_V}

        return loss, grads

    def ranks(self):
        """Compute log conditional probabilities (ranks).

        Returns
        -------
        ranks : np.ndarray of shape (n_microbes, n_metabolites)
            Row-centered log conditional probabilities.
        """
        U_aug, V_aug = self._build_augmented_matrices()
        logits = np.hstack(
            [np.zeros((self.n_microbes, 1)), U_aug @ V_aug]
        )
        # Center each row
        ranks = logits - logits.mean(axis=1, keepdims=True)
        return ranks

    def pack_params(self):
        """Flatten all parameters into a single vector.

        Returns
        -------
        theta : np.ndarray
            Flattened parameters [U, b_U, V, b_V].
        """
        return np.concatenate([
            self.U.ravel(),
            self.b_U.ravel(),
            self.V.ravel(),
            self.b_V.ravel(),
        ])

    def unpack_params(self, theta):
        """Restore parameter matrices from flat vector.

        Parameters
        ----------
        theta : np.ndarray
            Flattened parameters.
        """
        d1 = self.n_microbes
        d2 = self.n_metabolites
        p = self.n_components

        idx = 0
        self.U = theta[idx:idx + d1 * p].reshape(d1, p)
        idx += d1 * p
        self.b_U = theta[idx:idx + d1].reshape(d1, 1)
        idx += d1
        self.V = theta[idx:idx + p * (d2 - 1)].reshape(p, d2 - 1)
        idx += p * (d2 - 1)
        self.b_V = theta[idx:idx + (d2 - 1)].reshape(1, d2 - 1)

    def full_batch_loss_and_gradient(self, X_coo, Y):
        """Compute full-batch loss and gradient for L-BFGS.

        Parameters
        ----------
        X_coo : coo_matrix
            Microbe counts in COO format.
        Y : np.ndarray of shape (n_samples, n_metabolites)
            Metabolite counts.

        Returns
        -------
        loss : float
            Negative log-posterior.
        grad : np.ndarray
            Flattened gradient vector.
        """
        # Extract sparse indices and weights
        rows = X_coo.row  # sample indices (nnz,)
        cols = X_coo.col  # microbe indices (nnz,)
        weights = X_coo.data  # microbe counts (nnz,)

        # Build augmented matrices once
        U_aug, V_aug = self._build_augmented_matrices()

        # Precompute all logits: (d1, d2) with reference column = 0
        logits_core = U_aug @ V_aug  # (d1, d2-1)
        logits_all = np.hstack([np.zeros((self.n_microbes, 1)), logits_core])

        # Stable softmax for all microbes
        log_norm = logsumexp(logits_all, axis=1, keepdims=True)  # (d1, 1)
        probs_all = np.exp(logits_all - log_norm)  # (d1, d2)

        # Sample totals
        N = Y.sum(axis=1)  # (n_samples,)

        # === Vectorized loss computation ===
        # Gather logits and log_norm for each (sample, microbe) pair
        logits_batch = logits_all[cols, :]  # (nnz, d2)
        log_norm_batch = log_norm[cols, 0]  # (nnz,)
        Y_batch = Y[rows, :]  # (nnz, d2)
        N_batch = N[rows]  # (nnz,)

        # Log-likelihood: sum_j y_j * eta_j - N * logsumexp(eta)
        log_lik = (Y_batch * logits_batch).sum(axis=1) - N_batch * log_norm_batch
        loss = -np.dot(weights, log_lik)  # weighted sum

        # === Vectorized gradient computation ===
        # Delta: w * (y[1:] - N * pi[1:]) for each (sample, microbe) pair
        probs_batch = probs_all[cols, :]  # (nnz, d2)
        delta = weights[:, None] * (
            Y_batch[:, 1:] - N_batch[:, None] * probs_batch[:, 1:]
        )  # (nnz, d2-1)

        # dV: sum over all entries of outer(U[m], delta)
        # This is U[cols].T @ delta = (p, nnz) @ (nnz, d2-1) = (p, d2-1)
        U_batch = self.U[cols, :]  # (nnz, p)
        dV = -U_batch.T @ delta  # (p, d2-1)

        # db_V: sum of delta over all entries
        db_V = -delta.sum(axis=0, keepdims=True)  # (1, d2-1)

        # dU: scatter-add delta @ V.T to rows indexed by cols
        delta_V = delta @ self.V.T  # (nnz, p)
        dU = np.zeros_like(self.U)
        np.add.at(dU, cols, -delta_V)

        # db_U: scatter-add delta.sum(axis=1) to rows indexed by cols
        delta_sum = delta.sum(axis=1)  # (nnz,)
        db_U = np.zeros_like(self.b_U)
        np.add.at(db_U, (cols, 0), -delta_sum)

        # === Add prior terms (L2 regularization) ===
        u_diff = self.U - self.u_prior_mean
        b_u_diff = self.b_U - self.u_prior_mean
        v_diff = self.V - self.v_prior_mean
        b_v_diff = self.b_V - self.v_prior_mean

        loss += 0.5 * np.sum(u_diff**2) / self.u_prior_scale**2
        loss += 0.5 * np.sum(b_u_diff**2) / self.u_prior_scale**2
        loss += 0.5 * np.sum(v_diff**2) / self.v_prior_scale**2
        loss += 0.5 * np.sum(b_v_diff**2) / self.v_prior_scale**2

        dU += u_diff / self.u_prior_scale**2
        db_U += b_u_diff / self.u_prior_scale**2
        dV += v_diff / self.v_prior_scale**2
        db_V += b_v_diff / self.v_prior_scale**2

        # Pack gradients
        grad = np.concatenate([
            dU.ravel(),
            db_U.ravel(),
            dV.ravel(),
            db_V.ravel(),
        ])

        return loss, grad


class MMvecResults(SkbioObject):
    r"""Results from MMvec analysis.

    This class contains the learned embeddings and co-occurrence patterns
    from fitting an MMvec model. The key outputs enable both interpretation
    (which microbes co-occur with which metabolites) and prediction
    (expected metabolites given a microbial community).

    Attributes
    ----------
    microbe_embeddings : pd.DataFrame
        Microbe coordinates in latent space.
        Shape: (n_microbes, n_components + 1) where +1 is the bias term.

        Each row is a vector representation of a microbe. Microbes with
        similar embedding vectors tend to co-occur with similar sets of
        metabolites. The Euclidean distance or cosine similarity between
        microbe embeddings can be used to identify functionally related
        microbes. The final column ("bias") captures the baseline tendency
        of each microbe to associate with metabolites overall.

    metabolite_embeddings : pd.DataFrame
        Metabolite coordinates in latent space.
        Shape: (n_metabolites, n_components + 1) where +1 is the bias term.

        Each row is a vector representation of a metabolite. Metabolites
        with similar embedding vectors tend to co-occur with similar sets
        of microbes. The first row corresponds to the reference metabolite
        (all zeros) used for identifiability. The distance between
        metabolite embeddings indicates similarity in their microbial
        associations. The final column ("bias") captures the baseline
        abundance of each metabolite.

    ranks : pd.DataFrame
        Log conditional probability matrix (co-occurrence scores).
        Shape: (n_microbes, n_metabolites). Row-centered.

        Entry (i, j) represents the log-odds of observing metabolite j
        given microbe i, relative to the row mean. Higher values indicate
        stronger positive associations. This matrix is row-centered (each
        row sums to zero) for identifiability. To obtain actual conditional
        probabilities, use the :meth:`probabilities` method.

        The ranks matrix is the primary output for identifying
        microbe-metabolite associations. Sorting each row reveals which
        metabolites are most strongly associated with each microbe.

    convergence : pd.DataFrame
        Training diagnostics with columns:

        - ``iteration``: Iteration number (1-indexed).
        - ``loss``: Negative log-posterior (lower is better).

        Use this to diagnose training issues. The loss should generally
        decrease and stabilize. If the loss is still decreasing at the
        final iteration, consider increasing ``max_iter``. If the loss
        oscillates (Adam optimizer), try reducing ``learning_rate``.

    Notes
    -----
    **Detecting Overfitting with Q²**

    Overfitting occurs when the model memorizes training data rather than
    learning generalizable patterns. To detect overfitting:

    1. Split your data into training and test sets before fitting.
    2. Fit the model on training data only.
    3. Use :meth:`score` to compute Q² on held-out test data.

    Interpretation of Q² values:

    - **Q² close to 1**: Excellent predictive performance.
    - **Q² close to 0**: Model predicts no better than the mean.
    - **Q² negative**: Model performs worse than predicting the mean,
      indicating overfitting or model misspecification.

    If Q² is much lower than expected, try:

    - Reducing ``n_components`` (fewer latent dimensions).
    - Increasing regularization via smaller ``u_prior_scale`` and
      ``v_prior_scale`` values.
    - Collecting more training samples.

    **Embedding Interpretation**

    The embeddings place microbes and metabolites in the same latent space.
    The inner product between a microbe embedding and metabolite embedding
    (plus bias terms) gives the log-odds of their co-occurrence:

    .. math::

        \log \frac{P(m_j | \mu_i)}{P(m_{\text{ref}} | \mu_i)} =
        U_i \cdot V_j + b_{U_i} + b_{V_j}

    This means:

    - Microbes pointing in similar directions associate with similar
      metabolites.
    - Metabolites pointing in similar directions are produced/consumed
      by similar microbes.
    - The angle between a microbe and metabolite vector indicates their
      association strength.

    See Also
    --------
    mmvec : Fit an MMvec model.
    probabilities : Convert ranks to conditional probabilities.
    predict : Predict metabolite distributions for new samples.
    score : Evaluate predictive performance with Q².

    """

    def __init__(
        self,
        microbe_embeddings,
        metabolite_embeddings,
        ranks,
        convergence,
    ):
        self.microbe_embeddings = microbe_embeddings
        self.metabolite_embeddings = metabolite_embeddings
        self.ranks = ranks
        self.convergence = convergence

    def __str__(self):
        """Return string representation of MMvecResults."""
        n_microbes, n_metabolites = self.ranks.shape
        n_components = self.microbe_embeddings.shape[1] - 1  # exclude bias
        n_iterations = len(self.convergence)
        return (
            f"MMvecResults\n"
            f"  Microbes: {n_microbes}\n"
            f"  Metabolites: {n_metabolites}\n"
            f"  Components: {n_components}\n"
            f"  Iterations: {n_iterations}"
        )

    def probabilities(self):
        """Convert ranks to probability matrix via softmax.

        Returns
        -------
        probs : pd.DataFrame
            Conditional probabilities P(metabolite | microbe).
            Each row sums to 1.
        """
        probs = softmax(self.ranks.values)
        return pd.DataFrame(
            probs, index=self.ranks.index, columns=self.ranks.columns
        )

    def predict(self, microbes):
        """Predict metabolite distributions given microbe abundances.

        Computes the expected metabolite distribution for each sample by
        marginalizing over microbe compositions:

        P(metabolite) = sum_i P(microbe_i) * P(metabolite | microbe_i)

        Parameters
        ----------
        microbes : pd.DataFrame or array-like of shape (n_samples, n_microbes)
            Microbe abundance counts. Columns must match the microbes used
            during training.

        Returns
        -------
        predictions : pd.DataFrame
            Predicted metabolite proportions for each sample.
            Shape: (n_samples, n_metabolites). Each row sums to 1.

        Examples
        --------
        >>> from skbio.stats.multimodal import mmvec
        >>> import numpy as np
        >>> import pandas as pd
        >>> np.random.seed(42)
        >>> microbes = pd.DataFrame(
        ...     np.random.randint(1, 50, size=(20, 5)),
        ...     columns=[f'OTU_{i}' for i in range(5)]
        ... )
        >>> metabolites = pd.DataFrame(
        ...     np.random.randint(1, 50, size=(20, 8)),
        ...     columns=[f'met_{i}' for i in range(8)]
        ... )
        >>> result = mmvec(microbes, metabolites, n_components=2, max_iter=10)
        >>> # Predict on new samples
        >>> new_microbes = pd.DataFrame(
        ...     np.random.randint(1, 50, size=(5, 5)),
        ...     columns=[f'OTU_{i}' for i in range(5)]
        ... )
        >>> predictions = result.predict(new_microbes)
        >>> predictions.shape
        (5, 8)
        >>> np.allclose(predictions.sum(axis=1), 1.0)
        True

        """
        # Convert to array if needed
        if hasattr(microbes, "values"):
            X = microbes.values.astype(np.float64)
            sample_ids = list(microbes.index)
        else:
            X = np.asarray(microbes, dtype=np.float64)
            sample_ids = [f"sample_{i}" for i in range(X.shape[0])]

        # Normalize to proportions
        row_sums = X.sum(axis=1, keepdims=True)
        if np.any(row_sums == 0):
            raise ValueError(
                "microbes contains samples with all-zero counts. "
                "Remove these samples before calling predict."
            )
        microbe_props = X / row_sums

        # Get conditional probabilities P(metabolite | microbe)
        cond_probs = self.probabilities().values  # (n_microbes, n_metabolites)

        # Marginal: sum_microbe P(metabolite | microbe) * P(microbe)
        predicted = microbe_props @ cond_probs

        return pd.DataFrame(
            predicted,
            index=sample_ids,
            columns=self.ranks.columns,
        )

    def score(self, microbes, metabolites):
        r"""Compute Q² (coefficient of prediction) on held-out data.

        Q² measures predictive performance on test data, analogous to R²
        but for cross-validation. Values range from -inf to 1, where 1
        indicates perfect prediction and 0 indicates prediction no better
        than the mean.

        .. math::

            Q^2 = 1 - \frac{SS_{res}}{SS_{tot}}
                = 1 - \frac{\sum(y - \hat{y})^2}{\sum(y - \bar{y})^2}

        Parameters
        ----------
        microbes : pd.DataFrame or array-like of shape (n_samples, n_microbes)
            Test microbe abundance counts.
        metabolites : pd.DataFrame or array-like of shape (n_samples, n_metabolites)
            Test metabolite abundance counts.

        Returns
        -------
        q2 : float
            Q² score. Higher is better, with 1.0 being perfect prediction.

        See Also
        --------
        predict : Predict metabolite distributions.
        probabilities : Get conditional probability matrix.

        Examples
        --------
        >>> from skbio.stats.multimodal import mmvec
        >>> import numpy as np
        >>> import pandas as pd
        >>> np.random.seed(42)
        >>> # Training data
        >>> microbes = pd.DataFrame(
        ...     np.random.randint(1, 50, size=(30, 5)),
        ...     columns=[f'OTU_{i}' for i in range(5)]
        ... )
        >>> metabolites = pd.DataFrame(
        ...     np.random.randint(1, 50, size=(30, 8)),
        ...     columns=[f'met_{i}' for i in range(8)]
        ... )
        >>> result = mmvec(microbes, metabolites, n_components=2, max_iter=50)
        >>> # Evaluate on test data
        >>> test_microbes = pd.DataFrame(
        ...     np.random.randint(1, 50, size=(10, 5)),
        ...     columns=[f'OTU_{i}' for i in range(5)]
        ... )
        >>> test_metabolites = pd.DataFrame(
        ...     np.random.randint(1, 50, size=(10, 8)),
        ...     columns=[f'met_{i}' for i in range(8)]
        ... )
        >>> q2 = result.score(test_microbes, test_metabolites)
        >>> isinstance(q2, float)
        True

        """
        # Get predictions
        predicted = self.predict(microbes).values

        # Convert actual metabolites to proportions
        if hasattr(metabolites, "values"):
            Y = metabolites.values.astype(np.float64)
        else:
            Y = np.asarray(metabolites, dtype=np.float64)

        row_sums = Y.sum(axis=1, keepdims=True)
        if np.any(row_sums == 0):
            raise ValueError(
                "metabolites contains samples with all-zero counts. "
                "Remove these samples before calling score."
            )
        actual = Y / row_sums

        # Q² = 1 - SS_res / SS_tot
        ss_res = np.sum((actual - predicted) ** 2)
        ss_tot = np.sum((actual - actual.mean()) ** 2)

        if ss_tot == 0:
            # All actual values are the same - return 0 to avoid division by zero
            return 0.0

        return 1 - ss_res / ss_tot


def _adam_update(param, grad, m, v, t, lr, beta_1, beta_2, eps=1e-8):
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
    t : int
        Time step (1-indexed).
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
    m_hat = m / (1 - beta_1**t)
    v_hat = v / (1 - beta_2**t)

    param = param - lr * m_hat / (np.sqrt(v_hat) + eps)

    return param, m, v


def _clip_gradients(grads, clipnorm):
    """Apply global norm gradient clipping.

    Parameters
    ----------
    grads : dict
        Dictionary of gradients.
    clipnorm : float
        Maximum gradient norm.

    Returns
    -------
    grads : dict
        Clipped gradients.
    """
    global_norm = np.sqrt(
        sum(np.sum(g**2) for g in grads.values())
    )

    if global_norm > clipnorm:
        scale = clipnorm / global_norm
        grads = {k: v * scale for k, v in grads.items()}

    return grads


def _validate_inputs(microbes, metabolites, u_prior_scale, v_prior_scale):
    """Validate and convert input data for MMvec.

    Parameters
    ----------
    microbes : pd.DataFrame or array-like
        Microbe abundance counts.
    metabolites : pd.DataFrame or array-like
        Metabolite abundance counts.
    u_prior_scale : float
        Scale of Gaussian prior on U.
    v_prior_scale : float
        Scale of Gaussian prior on V.

    Returns
    -------
    X : np.ndarray
        Microbe counts as float64 array.
    Y : np.ndarray
        Metabolite counts as float64 array.
    microbe_ids : list
        Microbe feature IDs.
    metabolite_ids : list
        Metabolite feature IDs.

    Raises
    ------
    ValueError
        If inputs are invalid.
    """
    # Convert to arrays
    if hasattr(microbes, "values"):
        X = microbes.values.astype(np.float64)
        microbe_ids = list(microbes.columns)
    else:
        X = np.asarray(microbes, dtype=np.float64)
        microbe_ids = [f"microbe_{i}" for i in range(X.shape[1])]

    if hasattr(metabolites, "values"):
        Y = metabolites.values.astype(np.float64)
        metabolite_ids = list(metabolites.columns)
    else:
        Y = np.asarray(metabolites, dtype=np.float64)
        metabolite_ids = [f"metabolite_{i}" for i in range(Y.shape[1])]

    n_samples_X = X.shape[0]
    n_samples_Y = Y.shape[0]

    # Input validation
    if n_samples_X != n_samples_Y:
        raise ValueError(
            f"microbes and metabolites must have the same number of samples. "
            f"Got {n_samples_X} and {n_samples_Y}."
        )

    if u_prior_scale <= 0:
        raise ValueError(
            f"u_prior_scale must be positive, got {u_prior_scale}."
        )

    if v_prior_scale <= 0:
        raise ValueError(
            f"v_prior_scale must be positive, got {v_prior_scale}."
        )

    # Check for all-zero columns in microbes
    microbe_sums = X.sum(axis=0)
    zero_microbes = np.where(microbe_sums == 0)[0]
    if len(zero_microbes) > 0:
        zero_ids = [microbe_ids[i] for i in zero_microbes[:5]]
        msg = f"microbes contains all-zero columns: {zero_ids}"
        if len(zero_microbes) > 5:
            msg += f" and {len(zero_microbes) - 5} more"
        raise ValueError(msg + ". Remove these before calling mmvec.")

    # Check for all-zero columns in metabolites
    metabolite_sums = Y.sum(axis=0)
    zero_metabolites = np.where(metabolite_sums == 0)[0]
    if len(zero_metabolites) > 0:
        zero_ids = [metabolite_ids[i] for i in zero_metabolites[:5]]
        msg = f"metabolites contains all-zero columns: {zero_ids}"
        if len(zero_metabolites) > 5:
            msg += f" and {len(zero_metabolites) - 5} more"
        raise ValueError(msg + ". Remove these before calling mmvec.")

    # Check for all-zero rows (samples with no counts)
    microbe_row_sums = X.sum(axis=1)
    zero_samples = np.where(microbe_row_sums == 0)[0]
    if len(zero_samples) > 0:
        raise ValueError(
            f"microbes contains {len(zero_samples)} samples with all-zero "
            f"counts. Remove these samples before calling mmvec."
        )

    return X, Y, microbe_ids, metabolite_ids


def _train_lbfgs(model, X_coo, Y, max_iter, verbose):
    """Train MMvec model using L-BFGS-B optimization.

    Parameters
    ----------
    model : _MMvecModel
        Initialized model to train.
    X_coo : coo_matrix
        Microbe counts in sparse COO format.
    Y : np.ndarray
        Metabolite counts.
    max_iter : int
        Maximum number of L-BFGS iterations.
    verbose : bool
        Print training progress.

    Returns
    -------
    convergence_data : list of dict
        Training metrics per iteration.
    """
    convergence_data = []
    iteration_count = [0]  # Use list to allow modification in closure

    def objective(theta):
        model.unpack_params(theta)
        loss, grad = model.full_batch_loss_and_gradient(X_coo, Y)
        iteration_count[0] += 1

        convergence_data.append({
            "iteration": iteration_count[0],
            "loss": loss,
        })

        if verbose and iteration_count[0] % 10 == 0:
            print(f"Iteration {iteration_count[0]}, Loss: {loss:.4f}")

        return loss, grad

    # Initial parameters
    theta0 = model.pack_params()

    # Run L-BFGS-B
    result = minimize(
        objective,
        theta0,
        method="L-BFGS-B",
        jac=True,
        options={
            "maxiter": max_iter,
            "ftol": 1e-9,
            "gtol": 1e-5,
            "disp": False,
        },
    )

    # Unpack final parameters
    model.unpack_params(result.x)

    if verbose:
        print(f"L-BFGS-B converged: {result.success}, "
              f"iterations: {result.nit}, message: {result.message}")

    return convergence_data


def _train_adam(
    model,
    X,
    Y,
    X_coo,
    rng,
    max_iter,
    learning_rate,
    batch_size,
    beta_1,
    beta_2,
    clipnorm,
    batch_normalization,
    verbose,
):
    """Train MMvec model using Adam optimization.

    Parameters
    ----------
    model : _MMvecModel
        Initialized model to train.
    X : np.ndarray
        Microbe counts.
    Y : np.ndarray
        Metabolite counts.
    X_coo : coo_matrix
        Microbe counts in sparse COO format.
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
    batch_normalization : str
        Batch normalization mode ('unbiased' or 'legacy').
    verbose : bool
        Print training progress.

    Returns
    -------
    convergence_data : list of dict
        Training metrics per iteration.
    """
    convergence_data = []

    # Initialize Adam moments
    moments = {
        "U": (np.zeros_like(model.U), np.zeros_like(model.U)),
        "b_U": (np.zeros_like(model.b_U), np.zeros_like(model.b_U)),
        "V": (np.zeros_like(model.V), np.zeros_like(model.V)),
        "b_V": (np.zeros_like(model.b_V), np.zeros_like(model.b_V)),
    }

    # Compute number of iterations per epoch
    nnz = len(X_coo.data)
    iterations_per_epoch = max(1, nnz // batch_size)

    t = 0
    for epoch in range(max_iter):
        for _ in range(iterations_per_epoch):
            t += 1

            # Compute loss and gradients
            loss, grads = model.loss_and_gradients(
                X,
                Y,
                batch_size=batch_size,
                batch_normalization=batch_normalization,
                rng=rng,
            )

            # Gradient clipping
            grads = _clip_gradients(grads, clipnorm)

            # Adam updates
            for param_name in ["U", "b_U", "V", "b_V"]:
                param = getattr(model, param_name)
                m, v = moments[param_name]
                param, m, v = _adam_update(
                    param,
                    grads[param_name],
                    m,
                    v,
                    t,
                    learning_rate,
                    beta_1,
                    beta_2,
                )
                setattr(model, param_name, param)
                moments[param_name] = (m, v)

            convergence_data.append({"iteration": t, "loss": loss})

        if verbose:
            print(f"Epoch {epoch + 1}/{max_iter}, Loss: {loss:.4f}")

    return convergence_data


def _build_results(model, microbe_ids, metabolite_ids, convergence_data):
    """Build MMvecResults from trained model.

    Parameters
    ----------
    model : _MMvecModel
        Trained model.
    microbe_ids : list
        Microbe feature IDs.
    metabolite_ids : list
        Metabolite feature IDs.
    convergence_data : list of dict
        Training metrics.

    Returns
    -------
    MMvecResults
        Results object with embeddings and convergence data.
    """
    n_components = model.n_components

    # Compute ranks
    ranks = model.ranks()

    # Microbe embeddings: U with bias
    microbe_emb = np.hstack([model.U, model.b_U])
    pc_cols = [f"PC{i}" for i in range(n_components)] + ["bias"]
    microbe_embeddings = pd.DataFrame(
        microbe_emb, index=microbe_ids, columns=pc_cols
    )

    # Metabolite embeddings: V^T with bias
    # First row is reference (no embedding), rest are actual embeddings
    metabolite_emb = np.vstack([
        np.zeros((1, n_components + 1)),  # Reference category
        np.hstack([model.V.T, model.b_V.T]),
    ])
    metabolite_embeddings = pd.DataFrame(
        metabolite_emb, index=metabolite_ids, columns=pc_cols
    )

    ranks_df = pd.DataFrame(
        ranks, index=microbe_ids, columns=metabolite_ids
    )

    convergence = pd.DataFrame(convergence_data)

    return MMvecResults(
        microbe_embeddings=microbe_embeddings,
        metabolite_embeddings=metabolite_embeddings,
        ranks=ranks_df,
        convergence=convergence,
    )


def mmvec(
    microbes,
    metabolites,
    n_components=3,
    optimizer="lbfgs",
    max_iter=1000,
    learning_rate=1e-3,
    batch_size=50,
    u_prior_mean=0.0,
    u_prior_scale=1.0,
    v_prior_mean=0.0,
    v_prior_scale=1.0,
    beta_1=0.9,
    beta_2=0.95,
    clipnorm=10.0,
    batch_normalization="unbiased",
    random_state=None,
    verbose=False,
):
    r"""Multiomics Microbe-Metabolite Vectors (MMvec).

    Learns joint embeddings of two feature sets from their co-occurrence
    patterns using a multinomial likelihood model.

    While the parameter names use "microbes" and "metabolites" following the
    original publication, this method is **generic** and can be applied to any
    two omics modalities representable as compositional (count-based) data.
    For example: microbes and host transcripts, proteins and metabolites, or
    any pair of feature tables sharing the same samples.

    Parameters
    ----------
    microbes : pd.DataFrame or array-like of shape (n_samples, n_microbes)
        Abundance counts for the first modality (e.g., microbes, proteins).
        This modality is treated as the "conditioning" variable.
    metabolites : pd.DataFrame or array-like of shape (n_samples, n_metabolites)
        Abundance counts for the second modality (e.g., metabolites, transcripts).
        This modality is treated as the "conditioned" variable.
    n_components : int, optional
        Number of latent dimensions for embeddings. Default is 3.
    optimizer : {'lbfgs', 'adam'}, optional
        Optimization algorithm to use. Default is 'lbfgs'.

        - 'lbfgs': L-BFGS-B quasi-Newton method. Recommended for most cases.
          Typically converges in 50-200 iterations. Deterministic.
        - 'adam': Stochastic gradient descent with Adam. Use for very large
          datasets or when stochastic behavior is desired.
    max_iter : int, optional
        Maximum number of iterations. Default is 1000. For 'lbfgs', this is
        the max number of L-BFGS iterations. For 'adam', this is the number
        of epochs.
    learning_rate : float, optional
        Adam optimizer learning rate. Ignored for 'lbfgs'. Default is 1e-3.
    batch_size : int, optional
        Mini-batch size for Adam optimizer. Ignored for 'lbfgs'. Default is 50.
    u_prior_mean : float, optional
        Mean of Gaussian prior on first modality (microbes) embeddings.
        Default is 0.0.
    u_prior_scale : float, optional
        Scale (std) of Gaussian prior on first modality embeddings.
        Default is 1.0. Smaller values increase regularization.
    v_prior_mean : float, optional
        Mean of Gaussian prior on second modality (metabolites) embeddings.
        Default is 0.0.
    v_prior_scale : float, optional
        Scale (std) of Gaussian prior on second modality embeddings.
        Default is 1.0. Smaller values increase regularization.
    beta_1 : float, optional
        Adam exponential decay rate for first moment. Ignored for 'lbfgs'.
        Default is 0.9.
    beta_2 : float, optional
        Adam exponential decay rate for second moment. Ignored for 'lbfgs'.
        Default is 0.95.
    clipnorm : float, optional
        Gradient clipping threshold for Adam (global L2 norm). Ignored for
        'lbfgs'. Default is 10.0.
    batch_normalization : {'unbiased', 'legacy'}, optional
        Method for scaling mini-batch likelihood in Adam. Ignored for 'lbfgs'.
        Default is 'unbiased'.

        - 'unbiased': Uses norm = sum(microbe_counts) / batch_size.
        - 'legacy': Uses norm = n_samples / batch_size.
    random_state : int or numpy.random.Generator, optional
        Seed for random number generation or a Generator instance.
        If an int, creates a Generator with that seed. Default is None.
    verbose : bool, optional
        Print training progress. Default is False.

    Returns
    -------
    MMvecResults
        Object containing:

        - microbe_embeddings: DataFrame (n_microbes, n_components + 1)
        - metabolite_embeddings: DataFrame (n_metabolites, n_components + 1)
        - ranks: DataFrame (n_microbes, n_metabolites)
        - convergence: DataFrame with loss per iteration

    Notes
    -----
    The model learns:

    .. math::

        P(\text{metabolite}_j | \text{microbe}_i) =
        \text{softmax}(U_i \cdot V_j + b_{U_i} + b_{V_j})

    To evaluate model performance on held-out data, use the ``score`` method
    of the returned ``MMvecResults`` object.

    References
    ----------
    .. [1] Morton, J.T., et al. "Learning representations of
           microbe-metabolite interactions." Nature Methods, 2019.

    Examples
    --------
    >>> from skbio.stats.multimodal import mmvec
    >>> import numpy as np
    >>> import pandas as pd
    >>> # Create synthetic data
    >>> np.random.seed(42)
    >>> microbes = pd.DataFrame(
    ...     np.random.randint(0, 100, size=(50, 10)),
    ...     columns=[f'OTU_{i}' for i in range(10)]
    ... )
    >>> metabolites = pd.DataFrame(
    ...     np.random.randint(0, 100, size=(50, 15)),
    ...     columns=[f'metabolite_{i}' for i in range(15)]
    ... )
    >>> result = mmvec(microbes, metabolites, n_components=2, max_iter=10)
    >>> result.ranks.shape
    (10, 15)

    """
    # Create RNG - accept either an int seed or an existing Generator
    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)

    # Validate and convert inputs
    X, Y, microbe_ids, metabolite_ids = _validate_inputs(
        microbes, metabolites, u_prior_scale, v_prior_scale
    )

    n_microbes = X.shape[1]
    n_metabolites = Y.shape[1]

    # Validate optimizer
    optimizer = optimizer.lower()
    if optimizer not in ("lbfgs", "adam"):
        raise ValueError(
            f"optimizer must be 'lbfgs' or 'adam', got '{optimizer}'."
        )

    # Initialize model
    model = _MMvecModel(
        n_microbes=n_microbes,
        n_metabolites=n_metabolites,
        n_components=n_components,
        u_prior_mean=u_prior_mean,
        u_prior_scale=u_prior_scale,
        v_prior_mean=v_prior_mean,
        v_prior_scale=v_prior_scale,
        rng=rng,
    )

    # Convert to sparse COO format
    X_coo = coo_matrix(X)

    # Train model
    if optimizer == "lbfgs":
        convergence_data = _train_lbfgs(model, X_coo, Y, max_iter, verbose)
    else:
        convergence_data = _train_adam(
            model=model,
            X=X,
            Y=Y,
            X_coo=X_coo,
            rng=rng,
            max_iter=max_iter,
            learning_rate=learning_rate,
            batch_size=batch_size,
            beta_1=beta_1,
            beta_2=beta_2,
            clipnorm=clipnorm,
            batch_normalization=batch_normalization,
            verbose=verbose,
        )

    return _build_results(model, microbe_ids, metabolite_ids, convergence_data)
