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
        return_cv=False,
        X_test=None,
        Y_test=None,
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
        return_cv : bool
            Whether to return cross-validation error.
        X_test : array-like, optional
            Test microbe counts for CV.
        Y_test : array-like, optional
            Test metabolite counts for CV.
        rng : numpy.random.Generator, optional
            Random number generator for batch sampling.

        Returns
        -------
        loss : float
            Negative log posterior.
        grads : dict
            Gradients for U, b_U, V, b_V.
        cv_error : float, optional
            Cross-validation MAE if return_cv=True.
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

        if return_cv and X_test is not None and Y_test is not None:
            cv_error = self._compute_cv_error(X_test, Y_test)
            return loss, grads, cv_error

        return loss, grads

    def _compute_cv_error(self, X_test, Y_test):
        """Compute cross-validation mean absolute error.

        Parameters
        ----------
        X_test : array-like
            Test microbe counts.
        Y_test : array-like
            Test metabolite counts.

        Returns
        -------
        cv_error : float
            Mean absolute error on test set.
        """
        if issparse(X_test):
            X_coo = X_test.tocoo()
        else:
            X_coo = coo_matrix(X_test)

        Y_arr = np.asarray(Y_test)
        sample_ids = X_coo.row
        microbe_ids = X_coo.col

        # Forward pass
        logits = self.forward(microbe_ids)
        probs = softmax(logits)

        # Expected counts
        Y_batch = Y_arr[sample_ids, :]
        N = Y_batch.sum(axis=1, keepdims=True)
        expected = N * probs

        # Mean absolute error
        mae = np.mean(np.abs(expected - Y_batch))
        return mae

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


class MMvecResults(SkbioObject):
    """Results from MMvec analysis.

    Attributes
    ----------
    microbe_embeddings : pd.DataFrame
        Microbe coordinates in latent space.
        Shape: (n_microbes, n_components + 1) where +1 is bias.
    metabolite_embeddings : pd.DataFrame
        Metabolite coordinates in latent space.
        Shape: (n_metabolites, n_components + 1).
    ranks : pd.DataFrame
        Log conditional probabilities P(metabolite | microbe).
        Shape: (n_microbes, n_metabolites). Row-centered.
    convergence : pd.DataFrame
        Training metrics per iteration: iteration, loss, cv_error.
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


def mmvec(
    microbes,
    metabolites,
    n_components=3,
    epochs=100,
    batch_size=50,
    learning_rate=1e-3,
    u_prior_mean=0.0,
    u_prior_scale=1.0,
    v_prior_mean=0.0,
    v_prior_scale=1.0,
    beta_1=0.9,
    beta_2=0.95,
    clipnorm=10.0,
    batch_normalization="unbiased",
    test_microbes=None,
    test_metabolites=None,
    random_state=None,
    verbose=False,
):
    r"""Multiomics Microbe-Metabolite Vectors (MMvec).

    Learns joint embeddings of microbes and metabolites from their
    co-occurrence patterns using a multinomial likelihood model.

    Parameters
    ----------
    microbes : pd.DataFrame or array-like of shape (n_samples, n_microbes)
        Microbe abundance counts.
    metabolites : pd.DataFrame or array-like of shape (n_samples, n_metabolites)
        Metabolite abundance/intensity values.
    n_components : int, optional
        Number of latent dimensions for embeddings. Default is 3.
    epochs : int, optional
        Number of training epochs. Default is 100.
    batch_size : int, optional
        Mini-batch size for stochastic optimization. Default is 50.
    learning_rate : float, optional
        Adam optimizer learning rate. Default is 1e-3.
    u_prior_mean : float, optional
        Mean of Gaussian prior on microbe embeddings. Default is 0.0.
    u_prior_scale : float, optional
        Scale (std) of Gaussian prior on microbe embeddings. Default is 1.0.
    v_prior_mean : float, optional
        Mean of Gaussian prior on metabolite embeddings. Default is 0.0.
    v_prior_scale : float, optional
        Scale (std) of Gaussian prior on metabolite embeddings. Default is 1.0.
    beta_1 : float, optional
        Adam exponential decay rate for first moment. Default is 0.9.
    beta_2 : float, optional
        Adam exponential decay rate for second moment. Default is 0.95.
    clipnorm : float, optional
        Gradient clipping threshold (global L2 norm). Default is 10.0.
    batch_normalization : {'unbiased', 'legacy'}, optional
        Method for scaling mini-batch likelihood. Default is 'unbiased'.

        - 'unbiased': Uses norm = sum(microbe_counts) / batch_size.
        - 'legacy': Uses norm = n_samples / batch_size.
    test_microbes : pd.DataFrame, optional
        Test microbe counts for cross-validation.
    test_metabolites : pd.DataFrame, optional
        Test metabolite counts for cross-validation.
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
        - convergence: DataFrame with loss/cv_error per iteration

    Notes
    -----
    The model learns:

    .. math::

        P(\text{metabolite}_j | \text{microbe}_i) =
        \text{softmax}(U_i \cdot V_j + b_{U_i} + b_{V_j})

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
    >>> result = mmvec(microbes, metabolites, n_components=2, epochs=10)
    >>> result.ranks.shape
    (10, 15)

    """
    # Create a single RNG to use throughout the function.
    # Accept either an int seed or an existing Generator.
    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)

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

    n_samples_X, n_microbes = X.shape
    n_samples_Y, n_metabolites = Y.shape

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

    # Prepare test data if provided
    X_test = None
    Y_test = None
    if test_microbes is not None and test_metabolites is not None:
        if hasattr(test_microbes, "values"):
            X_test = test_microbes.values.astype(np.float64)
        else:
            X_test = np.asarray(test_microbes, dtype=np.float64)

        if hasattr(test_metabolites, "values"):
            Y_test = test_metabolites.values.astype(np.float64)
        else:
            Y_test = np.asarray(test_metabolites, dtype=np.float64)

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

    # Initialize Adam moments
    moments = {
        "U": (np.zeros_like(model.U), np.zeros_like(model.U)),
        "b_U": (np.zeros_like(model.b_U), np.zeros_like(model.b_U)),
        "V": (np.zeros_like(model.V), np.zeros_like(model.V)),
        "b_V": (np.zeros_like(model.b_V), np.zeros_like(model.b_V)),
    }

    # Compute number of iterations per epoch
    X_coo = coo_matrix(X)
    nnz = len(X_coo.data)
    iterations_per_epoch = max(1, nnz // batch_size)

    # Training loop
    convergence_data = []
    t = 0

    for epoch in range(epochs):
        for _ in range(iterations_per_epoch):
            t += 1

            # Compute loss and gradients
            if X_test is not None:
                loss, grads, cv_error = model.loss_and_gradients(
                    X,
                    Y,
                    batch_size=batch_size,
                    batch_normalization=batch_normalization,
                    return_cv=True,
                    X_test=X_test,
                    Y_test=Y_test,
                    rng=rng,
                )
            else:
                loss, grads = model.loss_and_gradients(
                    X,
                    Y,
                    batch_size=batch_size,
                    batch_normalization=batch_normalization,
                    rng=rng,
                )
                cv_error = None

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

            convergence_data.append(
                {"iteration": t, "loss": loss, "cv_error": cv_error}
            )

        if verbose:
            msg = f"Epoch {epoch + 1}/{epochs}, Loss: {loss:.4f}"
            if cv_error is not None:
                msg += f", CV Error: {cv_error:.4f}"
            print(msg)

    # Build results
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
