# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

import numpy as np
cimport numpy as cnp
from skbio.alignment import TabularMSA, AlignPath
from skbio.sequence import DNA
from numpy cimport int32_t, uint8_t, int8_t, float32_t, float64_t
cnp.import_array()

# Scoring constants
cdef uint8_t GAP = 1
cdef uint8_t GAP_SCORE = 255
cdef int32_t NEG_INF = -2147483647

# Structure definitions
cdef struct Index:
    Py_ssize_t i 
    Py_ssize_t j

cdef struct TracebackRes:
    float32_t score
    Py_ssize_t idx
    Index start
    Index end

def pairwise_align(seq1, seq2, subMatrix, gap_open, gap_extend, scope):
    """
    Public facing align function, runs Cythonized algorithm.

    Parameters:
    ----------
    seq1 : string or any scikit-bio DNA compatible class
        The first sequence, represented as a string, DNA, or any other DNA-compatible class.
    seq2 : string or any scikit-bio DNA compatible class
        The second sequence, represented as a string, DNA, or any other DNA-compatible class.
    subMatrix : skbio.sequence.SubstitutionMatrix
        Substitution matrix, taken from the scikit-bio sequence module.
    gap_open : float
        The penalty for opening a gap.
    gap_extend : float
        The penalty for extending an existing gap.
    scope : string
        Either "local" or "global", indicates alignment type

    Returns:
    -------
    TabularMSA
        TabularMSA that represents final alignment.
    """

    # Hash seq1 and seq2 through subMatrix
    asDNA = [DNA(seq1), DNA(seq2)]
    seq1_idx = subMatrix._char_hash[asDNA[0]._bytes]
    seq2_idx = subMatrix._char_hash[asDNA[1]._bytes]
    
    # Calculate max alignment length
    max_len = seq1_idx.shape[0] + seq2_idx.shape[0]

    # Initialize final alignment buffers
    aligned_seq1 = np.zeros(max_len, dtype=np.uint8)
    aligned_seq2 = np.zeros(max_len, dtype=np.uint8)

    # Create tuple of arguments and pass into algorithm
    args = (seq1_idx, seq2_idx, aligned_seq1, aligned_seq2, subMatrix._data.astype(float), gap_open, gap_extend)
    if(scope.lower() == "local"):
        args += (local_align_fill, local_align_trace)
    else:
        args += (global_align_fill, global_align_trace)
    res = arg_pass(*args)

    # Format result into a usable format and conver tto alignment path
    result = (aligned_seq1[res['idx']+1:aligned_seq1.size], aligned_seq2[res['idx']+1:aligned_seq2.size])
    if result[0].size == 0 and result[1].size == 1:
        raise Exception("No local alignment found!")
    path = AlignPath.from_bits(np.vstack(result))

    # Cut out local alignment path from originals
    if(scope.lower() == "local"):
        asDNA = [DNA(seq1[res['start']['i']:res['end']['i']]), DNA(seq2[res['start']['j']:res['end']['j']])]
    
    # Convert to TabularMSA and return
    return TabularMSA.from_path_seqs(path, asDNA), res['score']

def arg_pass(const uint8_t[::1] seq1, const uint8_t[::1] seq2, cnp.ndarray aligned_seq1, cnp.ndarray aligned_seq2, const float64_t[:, ::1] subMatrix, const int8_t gap_open, const int8_t gap_extend, fill_func, trace_func):
    '''
    Helper function to connect Python code with Cython code
    '''

    return _pairwise_align(seq1, seq2, aligned_seq1, aligned_seq2, subMatrix, gap_open, gap_extend, fill_func, trace_func)

cdef TracebackRes _pairwise_align(const uint8_t[::1] seq1, const uint8_t[::1] seq2, cnp.ndarray aligned_seq1, cnp.ndarray aligned_seq2, const float64_t[:, ::1] subMatrix, const int8_t gap_open, const int8_t gap_extend, fill_func, trace_func) noexcept:
    """
    Calculate alignment between seq1 and seq2 using given scope and penalties.

    Parameters:
    ----------
    seq1 : const cnp.uint8_t[::1]
        The first sequence hashed through subMatrix.
    seq2 : const cnp.uint8_t[::1]
        The second sequence hashed through subMatrix.
    aligned_seq1 : cnp.ndarray
        Array for storing the aligned version of `seq1`.
    aligned_seq2 : cnp.ndarray
        Array for storing the aligned version of `seq2`.
    subMatrix : const cnp.float64_t[:, :]
        The substitution matrix used for scoring alignments.
    gap_open : int8_t
        The penalty for opening a gap.
    gap_extend : int8_t
        The penalty for extending an existing gap.
    fill_func : function
        Function to fill P, D, and Q matrices.
    trace_func : function
        Function to trace result through filled matrices.

    Returns:
    -------
    TracebackRes
        A structure containing the score as well as the traceback end index.
    """
    
    # Initialize variables
    cdef Py_ssize_t m = seq1.shape[0]
    cdef Py_ssize_t n = seq2.shape[0]
    cdef Py_ssize_t max_len = m + n
    cdef Py_ssize_t idx
    cdef TracebackRes res

    # Initialize matrices
    cdef D = np.empty((m + 1, n + 1), dtype=np.float32)
    cdef P = np.empty((m + 1, n + 1), dtype=np.float32)
    cdef Q = np.empty((m + 1, n + 1), dtype=np.float32)

    # Get memoryview representations for Cython
    cdef float32_t[:, ::1] D_view = D
    cdef float32_t[:, ::1] P_view = P
    cdef float32_t[:, ::1] Q_view = Q
    cdef uint8_t[::1] aligned_seq1_view = aligned_seq1
    cdef uint8_t[::1] aligned_seq2_view = aligned_seq2

    # Run and return fill and trace functions on variables
    fill_func(D_view, P_view, Q_view, subMatrix, seq1, seq2, gap_open, gap_extend)
    res = trace_func(D_view, P_view, Q_view, max_len, aligned_seq1_view, aligned_seq2_view, seq1, seq2, subMatrix, gap_open, gap_extend)
    return res

cdef void global_align_fill(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, const cnp.float64_t[:, :] subMatrix, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const int8_t gap_open, const int8_t gap_extend) noexcept nogil:
    """
    Perform traceback for global sequence alignment using precomputed score matrices.

    Parameters:
    ----------
    D_view : float32_t[:, :]
        Empty match penalty matrix.
    P_view : float32_t[:, :]
        Empty insert penalty matrix.
    Q_view : float32_t[:, :]
        Empty deletion penalty matrix.
    subMatrix : const cnp.float64_t[:, :]
        The substitution matrix used for scoring alignments.
    seq1 : const cnp.uint8_t[::1]
        The first sequence hashed through subMatrix.
    seq2 : const cnp.uint8_t[::1]
        The second sequence hashed through subMatrix.
    gap_open : int8_t
        The penalty for opening a gap.
    gap_extend : int8_t
        The penalty for extending an existing gap.

    Returns:
    -------
    Nothing
    """

    # Initialize variables
    cdef Py_ssize_t i, j
    cdef Py_ssize_t m = seq1.shape[0]
    cdef Py_ssize_t n = seq2.shape[0]

    # Initialize 0, 0 for each matrix
    D_view[0, 0] = 0
    P_view[0, 0] = NEG_INF
    Q_view[0, 0] = NEG_INF

    # Initialize first row for each matrix
    for i in range(1, m + 1):
        D_view[i, 0] = gap_open + gap_extend * i
        Q_view[i, 0] = NEG_INF

    # Initialize first column for each matrix
    for j in range(1, n + 1):
        D_view[0, j] = gap_open + gap_extend * j
        P_view[0, j] = NEG_INF
    
    # Main fill loop
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Fill insertion scores
            P_view[i, j] = max(
                D_view[i - 1, j] + gap_open + gap_extend,
                P_view[i - 1, j] + gap_extend
            )
            # Fill deletion scores
            Q_view[i, j] = max(
                D_view[i, j - 1] + gap_open + gap_extend,
                Q_view[i, j - 1] + gap_extend
            )
            # Fill match/mismatch scores
            D_view[i, j] = max(
                D_view[i - 1, j - 1] + subMatrix[seq1[i-1], seq2[j-1]],
                P_view[i, j],
                Q_view[i, j]
            )

cdef void local_align_fill(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, const cnp.float64_t[:, :] subMatrix, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const int8_t gap_open, const int8_t gap_extend) noexcept nogil:
    """
    Perform traceback for global sequence alignment using precomputed score matrices.

    Parameters:
    ----------
    D_view : float32_t[:, :]
        Empty match penalty matrix.
    P_view : float32_t[:, :]
        Empty insert penalty matrix.
    Q_view : float32_t[:, :]
        Empty deletion penalty matrix.
    subMatrix : const cnp.float64_t[:, :]
        The substitution matrix used for scoring alignments.
    seq1 : const cnp.uint8_t[::1]
        The first sequence hashed through subMatrix.
    seq2 : const cnp.uint8_t[::1]
        The second sequence hashed through subMatrix.
    gap_open : int8_t
        The penalty for opening a gap.
    gap_extend : int8_t
        The penalty for extending an existing gap.

    Returns:
    -------
    Nothing
    """
    
    # Initialize variables
    cdef Py_ssize_t i, j, idx
    cdef Py_ssize_t m = seq1.shape[0]
    cdef Py_ssize_t n = seq2.shape[0]

    # Initialize 0, 0 for each matrix
    D_view[0, 0] = 0
    P_view[0, 0] = NEG_INF
    Q_view[0, 0] = NEG_INF

    # Initialize first row for each matrix
    for i in range(1, m + 1):
        D_view[i, 0] = 0
        Q_view[i, 0] = NEG_INF

    # Initialize first column for each matrix
    for j in range(1, n + 1):
        D_view[0, j] = 0
        P_view[0, j] = NEG_INF

    # Main fill loop
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Fill insertion scores
            P_view[i, j] = max(
                D_view[i - 1, j] + gap_open + gap_extend,
                P_view[i - 1, j] + gap_extend
            )
            # Fill deletion scores
            Q_view[i, j] = max(
                D_view[i, j - 1] + gap_open + gap_extend,
                Q_view[i, j - 1] + gap_extend
            )
            # Fill match/mismatch scores, can't go below 0
            D_view[i, j] = max(
                D_view[i - 1, j - 1] + subMatrix[seq1[i-1], seq2[j-1]],
                P_view[i, j],
                Q_view[i, j],
                0
            )

cdef TracebackRes global_align_trace(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, Py_ssize_t max_len, uint8_t[::1] aligned_seq1_view, uint8_t[::1] aligned_seq2_view, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const cnp.float64_t[:, :] subMatrix, const int8_t gap_open, const int8_t gap_extend) noexcept nogil:
    """
    Perform traceback for global sequence alignment using precomputed score matrices.

    Parameters:
    ----------
    D_view : float32_t[:, :]
        Filled match penalty matrix.
    P_view : float32_t[:, :]
        Filled insertion penalty matrix.
    Q_view : float32_t[:, :]
        Filled deletion penalty matrix.
    max_len : Py_ssize_t
        The maximum allowed length for the traceback.
    aligned_seq1_view : uint8_t[::1]
        Preallocated buffer for storing the aligned version of `seq1`.
    aligned_seq2_view : uint8_t[::1]
        Preallocated buffer for storing the aligned version of `seq2`.
    seq1 : const cnp.uint8_t[::1]
        The first sequence hashed through subMatrix.
    seq2 : const cnp.uint8_t[::1]
        The second sequence hashed through subMatrix.
    subMatrix : const cnp.float64_t[:, :]
        The substitution matrix used for scoring alignments.
    gap_open : int8_t
        The penalty for opening a gap.
    gap_extend : int8_t
        The penalty for extending an existing gap.

    Returns:
    -------
    TracebackRes
        A structure containing the score as well as the traceback end index.
    """

    # Initialize variables
    cdef float32_t score
    cdef uint8_t current_matrix
    cdef TracebackRes res
    cdef Py_ssize_t idx = max_len - 1
    cdef Py_ssize_t m = seq1.shape[0]
    cdef Py_ssize_t n = seq2.shape[0]

    # Check for matrix with highest score and set as starting point
    if D_view[m, n] >= P_view[m, n] and D_view[m, n] >= Q_view[m, n]:
        score = D_view[m, n]
        current_matrix = 0
    elif P_view[m, n] >= Q_view[m, n]:
        score = P_view[m, n]
        current_matrix = 1
    else:
        score = Q_view[m, n]
        current_matrix = 2

    # Main traceback while loop
    while m > 0 and n > 0:
        if current_matrix == 0:
            # Go diagonal if match/mismatch
            if D_view[m, n] == D_view[m-1, n-1] + subMatrix[seq1[m-1], seq2[n-1]]:
                m -= 1
                n -= 1
            # Swap to insert or delete matrix otherwise
            elif D_view[m, n] == P_view[m, n]:
                current_matrix = 1
                continue
            else:
                current_matrix = 2
                continue
        elif current_matrix == 1:
            # If gap is opened, go back to match/mismatch
            aligned_seq2_view[idx] = GAP
            if P_view[m, n] == D_view[m-1, n] + gap_open + gap_extend:
                current_matrix = 0
            m -= 1
        else:
            # If gap is opened, go back to match/mismatch
            aligned_seq1_view[idx] = GAP
            if Q_view[m, n] == D_view[m, n-1] + gap_open + gap_extend:
                current_matrix = 0 
            n -= 1
        # Move back in aligned sequences
        idx -= 1
    
    # Set any remaining spaces as gaps
    while m > 0:
        aligned_seq2_view[idx] = GAP
        m -= 1
        idx -= 1
    while n > 0:
        aligned_seq1_view[idx] = GAP
        n -= 1
        idx -= 1
    
    # Store results in return structure
    res.score = score
    res.idx = idx
    return res

cdef TracebackRes local_align_trace(float32_t[:, :] D_view, float32_t[:, :] P_view, float32_t[:, :] Q_view, Py_ssize_t max_len, uint8_t[::1] aligned_seq1_view, uint8_t[::1] aligned_seq2_view, const cnp.uint8_t[::1] seq1, const cnp.uint8_t[::1] seq2, const cnp.float64_t[:, :] subMatrix, const int8_t gap_open, const int8_t gap_extend) noexcept:
    """
    Perform traceback for local sequence alignment using precomputed score matrices.

    Parameters:
    ----------
    D_view : float32_t[:, :]
        Filled match penalty matrix.
    P_view : float32_t[:, :]
        Filled insertion penalty matrix.
    Q_view : float32_t[:, :]
        Filled deletion penalty matrix.
    loc : Index
        Index object, unused in local trace but required for parity.
    max_len : Py_ssize_t
        The maximum allowed length for the traceback.
    aligned_seq1_view : uint8_t[::1]
        Preallocated buffer for storing the aligned version of `seq1`.
    aligned_seq2_view : uint8_t[::1]
        Preallocated buffer for storing the aligned version of `seq2`.
    seq1 : const cnp.uint8_t[::1]
        The first sequence hashed through subMatrix.
    seq2 : const cnp.uint8_t[::1]
        The second sequence hashed through subMatrix.
    subMatrix : const cnp.float64_t[:, :]
        The substitution matrix used for scoring alignments.
    gap_open : int8_t
        The penalty for opening a gap.
    gap_extend : int8_t
        The penalty for extending an existing gap.

    Returns:
    -------
    TracebackRes
        A structure containing the score as well as the start and end indices of the traceback.
    """

    # Initialize variables
    cdef Index P_loc, D_loc, Q_loc, end
    cdef float32_t score
    cdef uint8_t current_matrix
    cdef TracebackRes res
    cdef Py_ssize_t idx = max_len - 1

    # Get highest score from each matrix
    P_loc.i, P_loc.j = np.unravel_index(np.argmax(P_view), np.shape(P_view))
    D_loc.i, D_loc.j = np.unravel_index(np.argmax(D_view), np.shape(D_view))
    Q_loc.i, Q_loc.j = np.unravel_index(np.argmax(Q_view), np.shape(Q_view))

    # Check for matrix with highest overall score and set as starting point
    if D_view[D_loc.i, D_loc.j] >= P_view[P_loc.i, P_loc.j] and D_view[D_loc.i, D_loc.j] >= Q_view[Q_loc.i, Q_loc.j]:
        current_matrix = 0
        score = D_view[D_loc.i, D_loc.j]
        loc = D_loc
    elif P_view[P_loc.i, P_loc.j] >= Q_view[Q_loc.i, Q_loc.j]:
        current_matrix = 1
        score = P_view[P_loc.i, P_loc.j]
        loc = P_loc
    else:
        current_matrix = 2
        score = Q_view[Q_loc.i, Q_loc.j]
        loc = Q_loc
    
    # Save endpoint for return
    end = loc

    # Main traceback while loop
    while loc.i > 0 and loc.j > 0:
        if current_matrix == 0:
            # Go diagonal if match/mismatch
            if D_view[loc.i, loc.j] == D_view[loc.i-1, loc.j-1] + subMatrix[seq1[loc.i-1], seq2[loc.j-1]]:
                loc.i -= 1
                loc.j -= 1
            elif D_view[loc.i, loc.j] == P_view[loc.i, loc.j]:
                current_matrix = 1
                continue
            else:
                current_matrix = 2
                continue
            # Break once the end of the traceback is reached
            if D_view[loc.i, loc.j] == 0:
                idx -= 1
                break
        elif current_matrix == 1:
            # If gap is opened, go back to match/mismatch
            aligned_seq2_view[idx] = GAP
            if P_view[loc.i, loc.j] == D_view[loc.i-1, loc.j] + gap_open + gap_extend:
                current_matrix = 0
            loc.i -= 1
        else:
            # If gap is opened, go back to match/mismatch
            aligned_seq1_view[idx] = GAP
            if Q_view[loc.i, loc.j] == D_view[loc.i, loc.j-1] + gap_open + gap_extend:
                current_matrix = 0 
            loc.j -= 1
        # Move back in aligned sequences
        idx -= 1

    # Store results in return structure
    res.score = score
    res.idx = idx
    res.start = loc
    res.end = end
    return res

# Alignment score wrapper function for export
def align_score(seq1, seq2, subMatrix, gap_open, gap_extend):
    """
    Computes the alignment score between two sequences using the given substitution matrix 
    and gap penalties.

    Parameters:
    ----------
    seq1 : string or any scikit-bio DNA compatible class
        The first sequence, represented as a string, DNA, or any other DNA-compatible class.
    seq2 : string or any scikit-bio DNA compatible class
        The second sequence, represented as a string, DNA, or any other DNA-compatible class.
    subMatrix : skbio.sequence.SubstitutionMatrix
        Substitution matrix, taken from the scikit-bio sequence module.
    gap_open : float
        The penalty for opening a gap.
    gap_extend : float
        The penalty for extending an existing gap.

    Returns:
    -------
    float
        The computed alignment score.
    """
    # Has sequences through substitution matrix
    seq1_idx = subMatrix._char_hash[DNA(seq1)._bytes]
    seq2_idx = subMatrix._char_hash[DNA(seq2)._bytes]
    # Run though helper function to get ouput
    return _align_score(seq1_idx, seq2_idx, subMatrix._data.astype(np.float32), gap_open, gap_extend)

# Alignment score main function
cdef float32_t _align_score(const uint8_t[::1] seq1, const uint8_t[::1] seq2, const float32_t[:, :] subMatrix, const float32_t gap_open, const float32_t gap_extend) noexcept nogil:
    # Initialize  variables
    cdef Py_ssize_t i
    cdef uint8_t state = 0
    cdef float32_t score = 0

    # Iterate through sequences
    for i in range(seq1.shape[0]):
        # Gap in seq1
        if seq1[i] == GAP_SCORE:
            # Either no current gap or gap in seq2
            score += gap_extend
            if state == 0 or state == 2:
                score += gap_open
                state = 1
        # Gap in seq2
        elif seq2[i] == GAP_SCORE:
            # Either no current gap or gap in seq1
            score += gap_extend
            if state == 0 or state == 1:
                score += gap_open
                state = 2
        # No gaps
        else:
            score += subMatrix[seq1[i], seq2[i]] 
            state = 0
    # Return final result
    return score