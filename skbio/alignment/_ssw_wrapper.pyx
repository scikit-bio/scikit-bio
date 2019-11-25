# -----------------------------------------------------------------------------
#  Copyright (c) 2013--, scikit-bio development team.
#
#  Distributed under the terms of the Modified BSD License.
#
#  The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from cpython cimport bool
import numpy as np
cimport numpy as cnp
from skbio.sequence import Protein, Sequence

cdef extern from "_lib/ssw.h":

    ctypedef struct s_align:
        cnp.uint16_t score1
        cnp.uint16_t score2
        cnp.int32_t ref_begin1
        cnp.int32_t ref_end1
        cnp.int32_t read_begin1
        cnp.int32_t read_end1
        cnp.int32_t ref_end2
        cnp.uint32_t* cigar
        cnp.int32_t cigarLen

    ctypedef struct s_profile:
        pass

    cdef s_profile* ssw_init(const cnp.int8_t* read,
                             const cnp.int32_t readLen,
                             const cnp.int8_t* mat,
                             const cnp.int32_t n,
                             const cnp.int8_t score_size)

    cdef void init_destroy(s_profile* p)

    cdef s_align* ssw_align(const s_profile* prof,
                            const cnp.int8_t* ref,
                            cnp.int32_t refLen,
                            const cnp.uint8_t weight_gapO,
                            const cnp.uint8_t weight_gapE,
                            const cnp.uint8_t flag,
                            const cnp.uint16_t filters,
                            const cnp.int32_t filterd,
                            const cnp.int32_t maskLen)

    cdef void align_destroy(s_align* a)

np_aa_table = np.array([
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
    23,  0, 20,  4,  3,  6, 13,  7,  8,  9, 23, 11, 10, 12,  2, 23,
    14,  5,  1, 15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
    23,  0, 20,  4,  3,  6, 13,  7,  8,  9, 23, 11, 10, 12,  2, 23,
    14,  5,  1, 15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23])

np_nt_table = np.array([
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  3,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,
    4,  4,  4,  4,  3,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4])

mid_table = np.array(['M', 'I', 'D'])


cdef class AlignmentStructure:
    """Wraps the result of an alignment c struct so it is accessible to Python

    Notes
    -----
    `cigar` may be empty depending on parameters used.

    `target_begin` and `query_begin` may be -1 depending on parameters used.

    Developer note: `read_sequence` is an alias for `query_sequence` used by
    ssw.c as is `reference_sequence` for `target_sequence`
    """
    cdef s_align *p
    cdef str read_sequence
    cdef str reference_sequence
    cdef int index_starts_at
    cdef str _cigar_string

    def __cinit__(self, read_sequence, reference_sequence, index_starts_at):
        # We use `read_sequence` and `reference_sequence` here as they are
        # treated sematically as a private output of ssw.c like the `s_align`
        # struct
        self.read_sequence = read_sequence
        self.reference_sequence = reference_sequence
        self.index_starts_at = index_starts_at

    cdef __constructor(self, s_align* pointer):
        self.p = pointer

    def __dealloc__(self):
        if self.p is not NULL:
            align_destroy(self.p)

    def __getitem__(self, key):
        return getattr(self, key)

    def __repr__(self):
        data = ['optimal_alignment_score', 'suboptimal_alignment_score',
                'query_begin', 'query_end', 'target_begin',
                'target_end_optimal', 'target_end_suboptimal', 'cigar',
                'query_sequence', 'target_sequence']
        return "{\n%s\n}" % ',\n'.join([
            "    {!r}: {!r}".format(k, self[k]) for k in data])

    def __str__(self):
        score = "Score: %d" % self.optimal_alignment_score
        if self.query_sequence and self.cigar:
            target = self.aligned_target_sequence
            query = self.aligned_query_sequence
            align_len = len(query)
            if align_len > 13:
                target = target[:10] + "..."
                query = query[:10] + "..."

            length = "Length: %d" % align_len
            return "\n".join([query, target, score, length])
        return score

    @property
    def optimal_alignment_score(self):
        """Optimal alignment score

        Returns
        -------
        int
            The optimal alignment score

        """
        return self.p.score1

    @property
    def suboptimal_alignment_score(self):
        """Suboptimal alignment score

        Returns
        -------
        int
            The suboptimal alignment score

        """
        return self.p.score2

    @property
    def target_begin(self):
        """Character index where the target's alignment begins

        Returns
        -------
        int
            The character index of the target sequence's alignment's beginning

        Notes
        -----
        The result is a 0-based index by default

        """
        return self.p.ref_begin1 + self.index_starts_at if (self.p.ref_begin1
                                                            >= 0) else -1

    @property
    def target_end_optimal(self):
        """Character index where the target's optimal alignment ends

        Returns
        -------
        int
            The character index of the target sequence's optimal alignment's
             end

        Notes
        -----
        The result is a 0-based index by default

        """
        return self.p.ref_end1 + self.index_starts_at

    @property
    def target_end_suboptimal(self):
        """Character index where the target's suboptimal alignment ends

        Returns
        -------
        int
            The character index of the target sequence's suboptimal alignment's
             end

        Notes
        -----
        The result is a 0-based index by default

        """
        return self.p.ref_end2 + self.index_starts_at

    @property
    def query_begin(self):
        """Returns the character index at which the query sequence begins

        Returns
        -------
        int
            The character index of the query sequence beginning

        Notes
        -----
        The result is a 0-based index by default

        """
        return self.p.read_begin1 + self.index_starts_at if (self.p.read_begin1
                                                             >= 0) else -1

    @property
    def query_end(self):
        """Character index at where query sequence ends

        Returns
        -------
        int
            The character index of the query sequence ending

        Notes
        -----
        The result is a 0-based index by default

        """
        return self.p.read_end1 + self.index_starts_at

    @property
    def cigar(self):
        """Cigar formatted string for the optimal alignment

        Returns
        -------
        str
            The cigar string of the optimal alignment

        Notes
        -----
        The cigar string format is described in [1]_ and [2]_.

        If there is no cigar or optimal alignment, this will return an empty
        string

        References
        ----------
        .. [1] http://genome.sph.umich.edu/wiki/SAM
        .. [2] http://samtools.github.io/hts-specs/SAMv1.pdf

        """
        # Memoization! (1/2)
        if self._cigar_string is not None:
            return self._cigar_string
        cigar_list = []
        for i in range(self.p.cigarLen):
            # stored the same as that in BAM format,
            # high 28 bits: length, low 4 bits: M/I/D (0/1/2)

            # Length, remove first 4 bits
            cigar_list.append(str(self.p.cigar[i] >> 4))
            # M/I/D, lookup first 4 bits in the mid_table
            cigar_list.append(mid_table[self.p.cigar[i] & 0xf])
        # Memoization! (2/2)
        self._cigar_string = "".join(cigar_list)
        return self._cigar_string

    @property
    def query_sequence(self):
        """Query sequence

        Returns
        -------
        str
            The query sequence

        """
        return self.read_sequence

    @property
    def target_sequence(self):
        """Target sequence

        Returns
        -------
        str
            The target sequence

        """
        return self.reference_sequence

    @property
    def aligned_query_sequence(self):
        """Returns the query sequence aligned by the cigar

        Returns
        -------
        str
            Aligned query sequence

        Notes
        -----
        This will return `None` if `suppress_sequences` was True when this
        object was created

        """
        if self.query_sequence:
            return self._get_aligned_sequence(self.query_sequence,
                                              self._tuples_from_cigar(),
                                              self.query_begin, self.query_end,
                                              "D")
        return None

    @property
    def aligned_target_sequence(self):
        """Returns the target sequence aligned by the cigar

        Returns
        -------
        str
            Aligned target sequence

        Notes
        -----
        This will return `None` if `suppress_sequences` was True when this
        object was created

        """
        if self.target_sequence:
            return self._get_aligned_sequence(self.target_sequence,
                                              self._tuples_from_cigar(),
                                              self.target_begin,
                                              self.target_end_optimal,
                                              "I")
        return None

    def set_zero_based(self, is_zero_based):
        """Set the aligment indices to start at 0 if True else 1 if False

        """
        if is_zero_based:
            self.index_starts_at = 0
        else:
            self.index_starts_at = 1

    def is_zero_based(self):
        """Returns True if alignment inidices start at 0 else False

        Returns
        -------
        bool
            Whether the alignment inidices start at 0

        """
        return self.index_starts_at == 0

    def _get_aligned_sequence(self, sequence, tuple_cigar, begin, end,
                              gap_type):
        # Save the original index scheme and then set it to 0 (1/2)
        orig_z_base = self.is_zero_based()
        self.set_zero_based(True)
        aligned_sequence = []
        seq = sequence[begin:end + 1]
        index = 0
        for length, mid in tuple_cigar:
            if mid == gap_type:
                aligned_sequence += ['-' * length]
            else:
                aligned_sequence += [seq[index:index + length]]
                index += length
        # Our sequence end is sometimes beyond the cigar:
        aligned_sequence += [seq[index:end - begin + 1]]
        # Revert our index scheme to the original (2/2)
        self.set_zero_based(orig_z_base)
        return "".join(aligned_sequence)

    def _tuples_from_cigar(self):
        tuples = []
        length_stack = []
        for character in self.cigar:
            if character.isdigit():
                length_stack.append(character)
            else:
                tuples.append((int("".join(length_stack)), character))
                length_stack = []
        return tuples

cdef class StripedSmithWaterman:
    """Performs a striped (banded) Smith Waterman Alignment.

    First a StripedSmithWaterman object must be instantiated with a query
    sequence. The resulting object is then callable with a target sequence and
    may be reused on a large collection of target sequences.

    Parameters
    ----------
    query_sequence : string
        The query sequence, this may be upper or lowercase from the set of
        {A, C, G, T, N} (nucleotide) or from the set of
        {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, B, Z, X, *
        } (protein)
    gap_open_penalty : int, optional
        The penalty applied to creating a gap in the alignment. This CANNOT
        be 0.
        Default is 5.
    gap_extend_penalty : int, optional
        The penalty applied to extending a gap in the alignment. This CANNOT
        be 0.
        Default is 2.
    score_size : int, optional
        If your estimated best alignment score is < 255 this should be 0.
        If your estimated best alignment score is >= 255, this should be 1.
        If you don't know, this should be 2.
        Default is 2.
    mask_length : int, optional
        The distance between the optimal and suboptimal alignment ending
        position >= mask_length. We suggest to use len(query_sequence)/2, if
        you don't have special concerns.
        Detailed description of mask_length: After locating the optimal
        alignment ending position, the suboptimal alignment score can be
        heuristically found by checking the second largest score in the array
        that contains the maximal score of each column of the SW matrix. In
        order to avoid picking the scores that belong to the alignments
        sharing the partial best alignment, SSW C library masks the reference
        loci nearby (mask length = mask_length) the best alignment ending
        position and locates the second largest score from the unmasked
        elements.
        Default is 15.
    mask_auto : bool, optional
        This will automatically set the used mask length to be
        max(int(len(`query_sequence`)/2), `mask_length`).
        Default is True.
    score_only : bool, optional
        This will prevent the best alignment beginning positions (BABP) and the
        cigar from being returned as a result. This overrides any setting on
        `score_filter`, `distance_filter`, and `override_skip_babp`. It has the
        highest precedence.
        Default is False.
    score_filter : int, optional
        If set, this will prevent the cigar and best alignment beginning
        positions (BABP) from being returned if the optimal alignment score is
        less than `score_filter` saving some time computationally. This filter
        may be overridden by `score_only` (prevents BABP and cigar, regardless
        of other arguments), `distance_filter` (may prevent cigar, but will
        cause BABP to be calculated), and `override_skip_babp` (will ensure
        BABP) returned.
        Default is None.
    distance_filter : int, optional
        If set, this will prevent the cigar from being returned if the length
        of the `query_sequence` or the `target_sequence` is less than
        `distance_filter` saving some time computationally. The results of
        this filter may be overridden by `score_only` (prevents BABP and cigar,
        regardless of other arguments), and `score_filter` (may prevent cigar).
        `override_skip_babp` has no effect with this filter applied, as BABP
        must be calculated to perform the filter.
        Default is None.
    override_skip_babp : bool, optional
        When True, the best alignment beginning positions (BABP) will always be
        returned unless `score_only` is set to True.
        Default is False.
    protein : bool, optional
        When True, the `query_sequence` and `target_sequence` will be read as
        protein sequence. When False, the `query_sequence` and
        `target_sequence` will be read as nucleotide sequence. If True, a
        `substitution_matrix` must be supplied.
        Default is False.
    match_score : int, optional
        When using a nucleotide sequence, the match_score is the score added
        when a match occurs. This is ignored if `substitution_matrix` is
        provided.
        Default is 2.
    mismatch_score : int, optional
        When using a nucleotide sequence, the mismatch is the score subtracted
        when a mismatch occurs. This should be a negative integer.
        This is ignored if `substitution_matrix` is provided.
        Default is -3.
    substitution_matrix : 2D dict, optional
        Provides the score for each possible substitution of sequence
        characters. This may be used for protein or nucleotide sequences. The
        entire set of possible combinations for the relevant sequence type MUST
        be enumerated in the dict of dicts. This will override `match_score`
        and `mismatch_score`. Required when `protein` is True.
        Default is None.
    suppress_sequences : bool, optional
        If True, the query and target sequences will not be returned for
        convenience.
        Default is False.
    zero_index : bool, optional
        If True, all inidices will start at 0. If False, all inidices will
        start at 1.
        Default is True.

    Notes
    -----
    This is a wrapper for the SSW package [1]_.

    `mask_length` has to be >= 15, otherwise the suboptimal alignment
    information will NOT be returned.

    `match_score` is a positive integer and `mismatch_score` is a negative
    integer.

    `match_score` and `mismatch_score` are only meaningful in the context of
    nucleotide sequences.

    A substitution matrix must be provided when working with protein sequences.

    References
    ----------
    .. [1] Zhao, Mengyao, Wan-Ping Lee, Erik P. Garrison, & Gabor T.
       Marth. "SSW Library: An SIMD Smith-Waterman C/C++ Library for
       Applications". PLOS ONE (2013). Web. 11 July 2014.
       http://www.plosone.org/article/info:doi/10.1371/journal.pone.0082138

    """
    cdef s_profile *profile
    cdef cnp.uint8_t gap_open_penalty
    cdef cnp.uint8_t gap_extend_penalty
    cdef cnp.uint8_t bit_flag
    cdef cnp.uint16_t score_filter
    cdef cnp.int32_t distance_filter
    cdef cnp.int32_t mask_length
    cdef str read_sequence
    cdef int index_starts_at
    cdef bool is_protein
    cdef bool suppress_sequences
    cdef cnp.ndarray __KEEP_IT_IN_SCOPE_read
    cdef cnp.ndarray __KEEP_IT_IN_SCOPE_matrix

    def __cinit__(self, query_sequence,
                  gap_open_penalty=5,  # BLASTN Default
                  gap_extend_penalty=2,  # BLASTN Default
                  score_size=2,  # BLASTN Default
                  mask_length=15,  # Minimum length for a suboptimal alignment
                  mask_auto=True,
                  score_only=False,
                  score_filter=None,
                  distance_filter=None,
                  override_skip_babp=False,
                  protein=False,
                  match_score=2,  # BLASTN Default
                  mismatch_score=-3,  # BLASTN Default
                  substitution_matrix=None,
                  suppress_sequences=False,
                  zero_index=True):
        # initalize our values
        self.read_sequence = query_sequence
        if gap_open_penalty <= 0:
            raise ValueError("`gap_open_penalty` must be > 0")
        self.gap_open_penalty = gap_open_penalty
        if gap_extend_penalty <= 0:
            raise ValueError("`gap_extend_penalty` must be > 0")
        self.gap_extend_penalty = gap_extend_penalty
        self.distance_filter = 0 if distance_filter is None else \
            distance_filter
        self.score_filter = 0 if score_filter is None else score_filter
        self.suppress_sequences = suppress_sequences
        self.is_protein = protein
        self.bit_flag = self._get_bit_flag(override_skip_babp, score_only)
        # http://www.cs.utexas.edu/users/EWD/transcriptions/EWD08xx/EWD831.html
        # Dijkstra knows what's up:
        self.index_starts_at = 0 if zero_index else 1
        # set up our matrix
        cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] matrix
        if substitution_matrix is None:
            if protein:
                raise Exception("Must provide a substitution matrix for"
                                " protein sequences")
            matrix = self._build_match_matrix(match_score, mismatch_score)
        else:
            matrix = self._convert_dict2d_to_matrix(substitution_matrix)
        # Set up our mask_length
        # Mask is recommended to be max(query_sequence/2, 15)
        if mask_auto:
            self.mask_length = len(query_sequence) / 2
            if self.mask_length < mask_length:
                self.mask_length = mask_length
        else:
            self.mask_length = mask_length

        cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] read_seq
        read_seq = self._seq_converter(query_sequence)

        cdef cnp.int32_t read_length
        read_length = len(query_sequence)

        cdef cnp.int8_t s_size
        s_size = score_size

        cdef cnp.int32_t m_width
        m_width = 24 if self.is_protein else 5

        cdef s_profile* p
        self.profile = ssw_init(<cnp.int8_t*> read_seq.data,
                                read_length,
                                <cnp.int8_t*> matrix.data,
                                m_width,
                                s_size)

        # A hack to keep the python GC from eating our data
        self.__KEEP_IT_IN_SCOPE_read = read_seq
        self.__KEEP_IT_IN_SCOPE_matrix = matrix

    def __call__(self, target_sequence):
        """Align `target_sequence` to `query_sequence`

        Parameters
        ----------
        target_sequence : str

        Returns
        -------
        skbio.alignment.AlignmentStructure
            The resulting alignment.

        """
        reference_sequence = target_sequence
        cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] reference
        reference = self._seq_converter(reference_sequence)

        cdef cnp.int32_t ref_length
        ref_length = len(reference_sequence)

        cdef s_align *align
        align = ssw_align(self.profile, <cnp.int8_t*> reference.data,
                          ref_length, self.gap_open_penalty,
                          self.gap_extend_penalty, self.bit_flag,
                          self.score_filter, self.distance_filter,
                          self.mask_length)

        # Cython won't let me do this correctly, so duplicate code ahoy:
        if self.suppress_sequences:
            alignment = AlignmentStructure("", "", self.index_starts_at)
        else:
            alignment = AlignmentStructure(self.read_sequence,
                                           reference_sequence,
                                           self.index_starts_at)
        alignment.__constructor(align)  # Hack to get a pointer through
        return alignment

    def __dealloc__(self):
        if self.profile is not NULL:
            init_destroy(self.profile)

    def _get_bit_flag(self, override_skip_babp, score_only):
        bit_flag = 0
        if score_only:
            return bit_flag
        if override_skip_babp:
            bit_flag = bit_flag | 0x8
        if self.distance_filter != 0:
            bit_flag = bit_flag | 0x4
        if self.score_filter != 0:
            bit_flag = bit_flag | 0x2
        if bit_flag == 0 or bit_flag == 8:
            bit_flag = bit_flag | 0x1
        return bit_flag

    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] _seq_converter(
            self,
            sequence):
        cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] seq
        seq = np.empty(len(sequence), dtype=np.int8)
        if self.is_protein:
            for i, char in enumerate(sequence):
                seq[i] = np_aa_table[ord(char)]
        else:
            for i, char in enumerate(sequence):
                seq[i] = np_nt_table[ord(char)]
        return seq

    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] \
            _build_match_matrix(self, match_score, mismatch_score):
        sequence_order = "ACGTN"
        dict2d = {}
        for row in sequence_order:
            dict2d[row] = {}
            for column in sequence_order:
                if column == 'N' or row == 'N':
                    dict2d[row][column] = 0
                else:
                    dict2d[row][column] = match_score if row == column \
                        else mismatch_score
        return self._convert_dict2d_to_matrix(dict2d)

    cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] \
            _convert_dict2d_to_matrix(self, dict2d):
        if self.is_protein:
            sequence_order = "ARNDCQEGHILKMFPSTWYVBZX*"
        else:
            sequence_order = "ACGTN"
        cdef int i = 0
        length = len(sequence_order)
        cdef cnp.ndarray[cnp.int8_t, ndim = 1, mode = "c"] py_list_matrix = \
            np.empty(length*length, dtype=np.int8)
        for row in sequence_order:
            for column in sequence_order:
                py_list_matrix[i] = dict2d[row][column]
                i += 1
        return py_list_matrix
