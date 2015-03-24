import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg


def ssvd(A, k=10, p=10, qiter=0, compute_v=True,  **kwargs):
    """Runs the SSVD method

    This method is based on the algorithm described in "Finding Structure with
    Randomness: Probabilistic Algorithms for Constructing Approximate Matrix
    Decompositions [1]. The method implemented is derived from the R
    implementation found in [2].

    Parameters
    ----------
    A : np.ndarray
        The result of PCoABase._A
    k : unsigned int, optional
        The number of eigenvectors and values to find. A lower k will result in
        lower quality resulting eignvectors and values.
    p : unsigned int, optional
        Oversampling parameter, this is added to k to boost accuracy.
    qiter : unsigned int, optional
        The number of iterations to perform.
    compute_v : bool, optional
        Whether or not to compute u v in addition to s. True by default.

    Returns
    -------
    eigvals: np.ndarray
        The resulting k eigenvalues.
    U_ssvd: np.ndarray
        The first set of resulting k eigenvectors.
    V_ssvd: np.ndarray
        The second set of resulting k eigenvectors.
        Returned only if compute_v is true

    References
    ----------
    .. [1] http://epubs.siam.org/doi/abs/10.1137/090771806
    .. [2] https://mahout.apache.org/users/dim-reduction/ssvd.html

    """

    A = np.atleast_2d(A)
    if A.ndim > 2:
        raise ValueError("Input matrix can only have two dimensions or less")

    m, n = A.shape

    p = min(min(m, n) - k, p)      # an mxn matrix M has at most p = min(m,n) unique
                                   # singular values
    r = k + p                      # rank plus oversampling parameter p

    omega = np.random.standard_normal(size = (n, r))   # generate random matrix omega
    y = np.dot(A, omega)           # compute a sample matrix Y: apply A to random
                                   # vectors to identify part of its range corresponding
                                   # to largest singular values
    Q, R = sp.linalg.qr(y)         # find an ON matrix st. Y = QQ'Y
    b = np.dot(Q.T, A)             # multiply A by Q whose columns form
                                   # an orthonormal basis for the range of Y

    #often, no iteration required to small error in eqn. 1.5
    for i in range(qiter):
        y = np.dot(A, b.T)
        Q, R = sp.linalg.qr(y)
        b = np.dot(Q.T, A)

    bbt = np.dot(b, b.T)
    eigvals, U_ssvd = scipy.sparse.linalg.eigsh(bbt, k)
    if compute_v: # compute svd of much smaller matrix b
        V_ssvd = np.dot(b.T, np.dot(U_ssvd, np.diag(1 / np.sqrt(eigvals))))
        U_ssvd = np.dot(Q, U_ssvd)
        return np.sqrt(eigvals), U_ssvd, V_ssvd.T
    else:
        U_ssvd = np.dot(Q, U_ssvd)
        return np.sqrt(eigvals), U_ssvd
