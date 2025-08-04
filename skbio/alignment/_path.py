# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
from typing import Optional, Union, Tuple, List, TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray

from skbio._base import SkbioObject
from skbio.util._decorator import classonlymethod

if TYPE_CHECKING:  # pragma: no cover
    from skbio.sequence import Sequence


# CIGAR codes indexed by states in PairAlignPath
_cigar_codes = np.array(["M", "I", "D", "P"])

# Mapping of CIGAR codes to states in PairAlignPath
_cigar_mapping = {
    77: 0,  # M
    73: 1,  # I
    68: 2,  # D
    80: 3,  # P
    61: 0,  # =
    88: 0,  # X
    78: 2,  # N
    83: 1,  # S
    72: 3,  # H
}


class AlignPath(SkbioObject):
    r"""Store an alignment path between sequences.

    It is a compact data structure that stores the operations of sequence alignment:
    inserting gaps into designated positions between subsequences. It does not store
    the sequence data.

    Parameters
    ----------
    lengths : array_like of int of shape (n_segments,)
        Length of each segment in the alignment.
    states : array_like of uint8 of shape (n_segments,) or (n_packs, n_segments)
        Packed bits representing character (0) or gap (1) status per sequence per
        segment in the alignment.
    ranges : array_like of int of shape (n_sequences, 2), optional
        Start and stop positions of each sequence in the alignment.
    starts : array_like of int of shape (n_sequences,), optional
        Start position of each sequence in the alignment.
    stops : array_like of int of shape (n_sequences,), optional
        Stop position of each sequence in the alignment.

        .. note::
           The ranges can be defined by supplying ``starts``, ``stops`` or ``ranges``.
           With either of the first two, the program will locate the other side based
           on ``lengths`` and ``states``. For the last (``ranges``), the program will
           NOT validate the correctness of the values but simply take them.

    See Also
    --------
    PairAlignPath
    skbio.alignment.TabularMSA
    skbio.sequence.Sequence

    Notes
    -----
    The underlying data structure of the ``AlignPath`` class efficiently represents a
    sequence alignment as two equal-length vectors: ``lengths`` and ``states``. The
    lengths vector contains the lengths of individual segments of the alignment with
    consistent gap status. The states vector contains the packed bits of gap (1) and
    and character (0) status for each position in the alignment.

    This data structure is calculated by performing
    :wiki:`run-length encoding <Run-length_encoding>` (RLE) on the alignment,
    considering each segment with consistent gap status to be a unit in the encoding.
    This resembles the CIGAR string (see :meth:`PairAlignPath.to_cigar`), and is
    generalized to an arbitrary number of sequences.

    An ``AlignPath`` object is detached from the original or aligned sequences and is
    highly memory efficient. The more similar the sequences are (i.e., the fewer gaps),
    the more compact this data structure is. In the worst case, this object consumes
    1/8 memory space of the aligned sequences.

    This class permits fully vectorized operations and enables efficient conversion
    between various formats such as aligned sequences, indices (Biotite), and
    coordinates (Biopython).

    In addition to alignment operations, an ``AlignPath`` object also stores the ranges
    (start and stop positions) of the aligned region in the original sequences. This
    facilitates extraction of aligned sequences. The positions are 0-based and
    half-open, consistent with Python indexing, and compatible with the
    :wiki:`BED format <BED_(file_format)>`.

    Examples
    --------
    Create an ``AlignPath`` object from a ``TabularMSA`` object with three DNA
    sequences and 20 positions.

    >>> from skbio import DNA, TabularMSA
    >>> from skbio.alignment import AlignPath
    >>> msa = TabularMSA([
    ...    DNA('CGGTCGTAACGCGTA---CA'),
    ...    DNA('CAG--GTAAG-CATACCTCA'),
    ...    DNA('CGGTCGTCAC-TGTACACTA'),
    ... ])
    >>> path = AlignPath.from_tabular(msa)
    >>> path
    <AlignPath, sequences: 3, positions: 20, segments: 7>

    >>> path.lengths
    array([3, 2, 5, 1, 4, 3, 2])

    >>> path.states
    array([[0, 2, 0, 6, 0, 1, 0]], dtype=uint8)

    In the above example, the first three positions of the alignment contain no gaps,
    so the first value in the lengths array is 3, and that in the states array is 0.
    The fourth segment, which has length 1, would have gap status (0, 1, 1), which
    then becomes 6 after bit packing.

    An ``AlignPath`` object is rarely created from scratch. But one still could, like:

    >>> path = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
    ...                  states=[0, 2, 0, 6, 0, 1, 0],
    ...                  starts=[5, 1, 0])

    The parameter ``starts`` defines the start positions of the aligned region of each
    sequence. The program will automatically calculate the stop positions.

    >>> path.starts
    array([5, 1, 0])

    >>> path.stops
    array([22, 18, 19])

    >>> path.ranges
    array([[ 5, 22],
           [ 1, 18],
           [ 0, 19]])

    With the ranges, one can extract aligned subsequences from the original sequences.

    >>> seqs = [
    ...     DNA("NNNNNCGGTCGTAACGCGTACANNNNNNN"),
    ...     DNA("NCAGGTAAGCATACCTCA"),
    ...     DNA("CGGTCGTCACTGTACACTANN"),
    ... ]
    >>> for seq, (start, stop) in zip(seqs, path.ranges):
    ...     print(seq[start:stop])
    CGGTCGTAACGCGTACA
    CAGGTAAGCATACCTCA
    CGGTCGTCACTGTACACTA

    Alternatively, one can extract the aligned sequences with gap characters:

    >>> print(*path.to_aligned(seqs), sep='\n')
    CGGTCGTAACGCGTA---CA
    CAG--GTAAG-CATACCTCA
    CGGTCGTCAC-TGTACACTA

    """

    def __init__(
        self,
        lengths: ArrayLike,
        states: ArrayLike,
        *,
        ranges: Optional[ArrayLike] = None,
        starts: Optional[ArrayLike] = None,
        stops: Optional[ArrayLike] = None,
    ) -> None:
        self._lengths = np.asarray(lengths, dtype=np.intp)
        if self._lengths.ndim > 1:
            raise TypeError("`lengths` must be a 1-D array.")

        self._states = np.atleast_2d(np.asarray(states, dtype=np.uint8))
        if self._states.ndim > 2:
            raise TypeError("`states` must be a 1-D or 2-D array.")

        if self._lengths.shape[0] != self._states.shape[1]:
            raise ValueError(
                f"Numbers of segments in `lengths` ({self._lengths.shape[0]}) "
                f"and `states` ({self._states.shape[1]}) do not match."
            )

        if ranges is not None:
            self._ranges = np.asarray(ranges, dtype=np.intp)
            if self._ranges.ndim != 2 or self._ranges.shape[1] != 2:
                raise TypeError("`ranges` must be a 2-column array.")
        else:
            if starts is not None:
                starts_ = np.asarray(starts, dtype=np.intp)
                if starts_.ndim != 1:
                    raise TypeError("`starts` must be a 1-D array.")
                self._ranges = np.column_stack(
                    (starts_, starts_ + self._to_sizes(starts_.shape[0]))
                )
            elif stops is not None:
                stops_ = np.asarray(stops, dtype=np.intp)
                if stops_.ndim != 1:
                    raise TypeError("`stops` must be a 1-D array.")
                self._ranges = np.column_stack(
                    (stops_ - self._to_sizes(stops_.shape[0]), stops_)
                )
            else:
                raise ValueError("`ranges`, `starts` or `stops` must be provided.")

        # Number of sequences needs to be explicitly provided, because the packed bits
        # does not contain this information. (It is merely in multiples of 8.)
        nseq = self._ranges.shape[0]
        if np.ceil(nseq / 8) != self._states.shape[0]:
            max_seqs = self._states.shape[0] * 8
            raise ValueError(
                f"Number of sequences in ranges ({nseq}) and capacity of "
                f"states ({max_seqs - 7} to {max_seqs}) do not match."
            )

        # Shape is n_sequences (rows) x n_positions (columns), which is consistent with
        # TabularMSA
        npos = int(self._lengths.sum())
        self._shape = (nseq, npos)

    def __str__(self) -> str:
        r"""Return string representation of this alignment path."""
        return self.__repr__()

    def __repr__(self) -> str:
        r"""Return summary of the alignment path."""
        return (
            f"<{self.__class__.__name__}, "
            f"sequences: {self._shape[0]}, "
            f"positions: {self._shape[1]}, "
            f"segments: {self._lengths.size}>"
        )

    @property
    def lengths(self) -> NDArray[np.intp]:
        """Array of lengths of segments in alignment path."""
        return self._lengths

    @property
    def states(self) -> NDArray[np.uint8]:
        """Array of gap status of segments in alignment path."""
        return self._states

    @property
    def ranges(self) -> NDArray[np.intp]:
        """Array of (start, stop) positions of sequences in the alignment."""
        return self._ranges

    @property
    def starts(self) -> NDArray[np.intp]:
        """Array of start positions of sequences in the alignment."""
        return self._ranges[:, 0]

    @property
    def stops(self) -> NDArray[np.intp]:
        """Array of stop positions of sequences in the alignment."""
        return self._ranges[:, 1]

    @property
    def shape(self) -> Tuple[int, int]:
        """Number of sequences (rows) and positions (columns)."""
        return self._shape

    def _to_bits(self, count=None):
        r"""Unpack the alignment path into an array of bits by segment."""
        if count is None:
            count = self._shape[0]
        return np.unpackbits(self._states, axis=0, count=count, bitorder="little")

    def _to_sizes(self, count=None):
        r"""Calculate the size of aligned region within each sequence."""
        return (self._lengths * (1 - self._to_bits(count=count))).sum(axis=1)

    def to_bits(self, expand: bool = True) -> NDArray[np.uint8]:
        r"""Unpack the alignment path into an array of bits.

        .. versionchanged:: 0.7.0
            The default behavior now returns positions.

        Parameters
        ----------
        expand : bool, optional
            If True (default), each column in the returned array represents a position
            in the original alignment. If False, each column represents a segment in
            the alignment path.

            .. versionadded:: 0.7.0

        Returns
        -------
        ndarray of (0, 1) of shape (n_sequences, n_positions or n_segments)
            Array of zeros (character) and ones (gap) which represent the alignment.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> from skbio.alignment import AlignPath
        >>> seqs = [
        ...    DNA('CGTCGTGC'),
        ...    DNA('CA--GT-C'),
        ...    DNA('CGTCGT-T')
        ... ]
        >>> msa = TabularMSA(seqs)
        >>> path = AlignPath.from_tabular(msa)

        Return a bit array representing positions:

        >>> path.to_bits()
        array([[0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 1, 1, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 1, 0]], dtype=uint8)

        Return a bit array representing segments:

        >>> path.to_bits(expand=False)
        array([[0, 0, 0, 0, 0],
               [0, 1, 0, 1, 0],
               [0, 0, 0, 1, 0]], dtype=uint8)

        """
        bits = self._to_bits()
        return np.repeat(bits, self._lengths, axis=1) if expand else bits

    @classonlymethod
    def from_bits(
        cls, bits: ArrayLike, starts: Optional[ArrayLike] = None
    ) -> "AlignPath":
        r"""Create an alignment path from a bit array (0 - character, 1 - gap).

        Parameters
        ----------
        bits : array_like of (0, 1) of shape (n_sequences, n_positions)
            Array of zeros (character) and ones (gap) which represent the alignment.
        starts : array_like of int of shape (n_sequences,), optional
            Start position of each sequence in the alignment. If omitted, will set as
            zeros.

        Returns
        -------
        AlignPath
            The alignment path created from the given bit array.

        Examples
        --------
        >>> import numpy as np
        >>> from skbio.alignment import AlignPath
        >>> bit_arr = np.array([[0, 0, 0, 0, 0, 0, 0, 0],
        ...                     [0, 0, 1, 1, 0, 0, 1, 0],
        ...                     [0, 0, 0, 0, 0, 0, 1, 0]])
        >>> path = AlignPath.from_bits(bit_arr)
        >>> path
        <AlignPath, sequences: 3, positions: 8, segments: 5>

        """
        # pack bits into integers
        ints = np.packbits(bits, axis=0, bitorder="little")

        # get indices where segments start
        idx = np.append(0, np.where((ints[:, :-1] != ints[:, 1:]).any(axis=0))[0] + 1)

        # get lengths of segments
        lens = np.append(idx[1:] - idx[:-1], ints.shape[1] - idx[-1])

        # keep indices of segment starts
        ints = ints[:, idx]

        # set start positions as zeros if not specified
        if starts is None:
            starts = np.zeros(bits.shape[0], dtype=int)

        # return per-segment lengths and states
        return cls(lens, ints, starts=starts)

    def to_aligned(self, seqs, gap_char="-", flanking=None):
        r"""Extract aligned regions from original sequences.

        .. versionadded:: 0.7.0

        Parameters
        ----------
        seqs : iterable of Sequence or str
            Original sequences.
        gap_char : str, optional
            Character to be placed in each gap position. Default is "-". Set as "" to
            suppress gaps in the output.
        flanking : int or (int, int), optional
            Length of flanking regions in the original sequences to be included in the
            output. Can be two numbers (leading and trailing, respectively) or one
            number (same for leading and trailing). If the specified flanking region
            is longer than a sequence actually has, the remaining space will be filled
            with white spaces (" ").

        Returns
        -------
        list of str
            Aligned regions of the sequences.

        Raises
        ------
        ValueError
            If there are more sequences than in the path.
        ValueError
            If any sequence is shorter than in the path.

        See Also
        --------
        from_aligned
        skbio.alignment.TabularMSA.from_path_seqs

        Notes
        -----
        This method provides a convenient way to process and display alignments,
        without invoking the explicit ``TabularMSA`` class. Both ``Sequence`` objects
        and plain strings are valid input sequences.

        However, it only outputs strings without retaining the ``Sequence`` objects and
        their metadata. For the later purpose, please use ``TabularMSA``'s
        :meth:`~skbio.alignment.TabularMSA.from_path_seqs` method instead.

        Examples
        --------
        >>> from skbio.sequence import DNA
        >>> from skbio.alignment import AlignPath
        >>> path = AlignPath(
        ...     lengths=[2, 2, 2, 1, 1],
        ...     states=[0, 2, 0, 6, 0],
        ...     starts=[0, 3, 0],
        ... )
        >>> seqs = [
        ...    DNA('CGTCGTGC'),
        ...    DNA('ATTCAGTCGG'),
        ...    DNA('CGTCGTTAA')
        ... ]
        >>> path.to_aligned(seqs)
        ['CGTCGTGC',
         'CA--GT-C',
         'CGTCGT-T']

        """
        starts = self._ranges[:, 0]
        lens = self._lengths
        bits = self._to_bits()

        if isinstance(flanking, tuple):
            leading, trailing = flanking
        else:
            leading, trailing = flanking, flanking

        res = []
        for i, seq in enumerate(map(str, seqs)):
            try:
                pos = starts[i]
            except IndexError:
                raise ValueError("There are more sequences than in the path.")
            aln = ""

            # leading flanking region
            if leading:
                offset = pos - leading
                if offset >= 0:
                    aln += seq[offset:pos]
                else:
                    aln += " " * -offset + seq[:pos]

            # alignment region
            for L, gap in zip(lens, bits[i]):
                if gap:
                    aln += gap_char * L
                else:
                    new_pos = pos + L
                    aln += seq[pos:new_pos]
                    pos = new_pos

            remaining = len(seq) - pos
            if remaining < 0:
                raise ValueError(f"Sequence {i} is shorter than in the path.")

            # trailing flanking region
            if trailing:
                offset = remaining - trailing
                if offset >= 0:
                    aln += seq[pos : pos + trailing]
                else:
                    aln += seq[pos:] + " " * -offset

            res.append(aln)
        return res

    @classonlymethod
    def from_aligned(cls, aln, gap_chars="-", starts=None) -> "AlignPath":
        r"""Create an alignment path from aligned sequences.

        Parameters
        ----------
        aln : iterable of :class:`~skbio.sequence.Sequence`, str or sequence
            Aligned sequences. Can be skbio sequences, strings or sequences of any
            scalars.
        gap_chars : str or container, optional
            Characters that should be treated as gaps in aligned sequences. Default
            is "-".
        starts : array_like of int of shape (2,), optional
            Start positions of sequences. If omitted, will be all zeros.

        Returns
        -------
        AlignPath
            The alignment path created from the aligned sequences.

        See Also
        --------
        to_aligned
        from_tabular

        Notes
        -----
        This method is more general but less efficient than ``from_tabular``. It works
        with various sequence formats.

        Examples
        --------
        >>> from skbio.alignment import AlignPath
        >>> aln = [
        ...    'CGTCGTGC',
        ...    'CA--GT-C',
        ...    'CGTCGT-T'
        ... ]
        >>> path = AlignPath.from_aligned(aln)
        >>> path
        <AlignPath, sequences: 3, positions: 8, segments: 5>

        """
        from skbio.sequence import Sequence

        gaps = []
        for seq in aln:
            if isinstance(seq, Sequence):
                seq = str(seq)
            row = [x in gap_chars for x in seq]
            gaps.append(np.array(row, dtype=int))
        try:
            gaps = np.vstack(gaps)
        except ValueError:
            raise ValueError("Sequence lengths do not match.")
        return cls.from_bits(gaps, starts=starts)

    def _to_matrices(self, seqs, gap_code=45):
        r"""Generate matrices representing the alignment.

        Parameters
        ----------
        seqs : list of ndarray
            Original sequences as bytes.
        gap_code : uint8, optional
            Code to fill in gap positions in the character matrix. Default is 45 (-).

        Returns
        -------
        ndarray of shape (n_sequences, n_positions)
            Matrix of character codes.
        ndarray of bool of shape (n_sequences, n_positions)
            Matrix of gap status.
        ndarray of shape (n_sequences, n_segments)
            Matrix of segment status.
        ndarray of int of shape (n_segments,)
            Vector of segment lengths.

        Notes
        -----
        This method is currently private. Its only use case is to prepare data for
        `align_score`. The code is efficient for this purpose because it avoids all
        unnecessary calculations.

        If one intends to index the output with an ``SubstitutionMatrix`` object,
        having a ``gap_code`` < 128 is necessary. Otherwise, setting ``gap_code`` as
        None can save some compute.

        """
        # See also `to_bits` and `stops`. The following code mixes both to enhance
        # performance.
        bits = self._to_bits()
        lens = self._lengths

        # locate aligned region
        starts = self._ranges[:, 0]
        stops = self._ranges[:, 1]

        # create gap array
        gaps = np.repeat(bits, lens, axis=1).astype(bool)

        # allocate byte array
        chars = np.empty(self._shape, dtype=seqs[0].dtype)

        # fill in gaps (optional)
        if gap_code is not None:
            chars[gaps] = gap_code

        # fill in characters
        try:
            chars[~gaps] = np.concatenate(
                [seqs[i][starts[i] : stops[i]] for i in range(self._shape[0])]
            )
        except IndexError:
            raise ValueError("Fewer sequences were provided than in alignment path.")
        except ValueError:
            raise ValueError("Some sequences are shorter than in alignment path.")

        return chars, gaps, bits, lens

    @classonlymethod
    def from_tabular(cls, msa, starts=None) -> "AlignPath":
        r"""Create an alignment path from a `TabularMSA` object.

        Parameters
        ----------
        msa : TabularMSA
            TabularMSA to be converted into AlignPath object.
        starts : array_like of int of shape (2,), optional
            Start positions of sequences. If omitted, will be all zeros.

        Returns
        -------
        AlignPath
            The alignment path created from the TabularMSA object.

        See Also
        --------
        from_aligned
        skbio.TabularMSA.from_path_seqs

        Notes
        -----
        This method is more efficient and more specific than ``from_aligned``.

        """
        # Convert TabularMSA into a 2D array of bytes.
        # TODO: TabularMSA itself should be refactored to have this as the default data
        # structure.
        byte_arr = np.stack([x._bytes for x in msa._seqs])

        # Get gap character code.
        gap_chars = [ord(x) for x in msa.dtype.gap_chars]

        # Identify gap positions, and convert them into a bit array, then create an
        # alignment path based on it.
        return cls.from_bits(np.isin(byte_arr, gap_chars), starts=starts)

    def to_indices(self, gap=-1):
        r"""Generate an array of indices of characters in the original sequences.

        Parameters
        ----------
        gap : int, np.nan, np.inf, "del", or "mask", optional
            Method to encode gaps in the alignment. If numeric, replace gaps with this
            value. If "del", delete columns that have any gap. If "mask", return an
            ``np.ma.MaskedArray``, with gaps masked. Default is -1.

        Returns
        -------
        ndarray of int of shape (n_sequences, n_positions)
            Array of indices of characters in the original sequences.

        See Also
        --------
        from_indices

        Notes
        -----
        The transpose of the output matches the underlying data structure of Biotite's
        ``Alignment`` class [1]_. Therefore, one can convert scikit-bio alignments into
        Biotite alignments, and vice versa.

        References
        ----------
        .. [1] `biotite.sequence.align.Alignment <https://www.biotite-python.org/
           latest/apidoc/biotite.sequence.align.Alignment.html>`_

        Examples
        --------
        >>> from skbio.alignment import AlignPath
        >>> path = AlignPath(lengths=[2, 1, 2, 1],
        ...                  states=[0, 6, 0, 1],
        ...                  starts=[0, 1, 2])
        >>> idx = path.to_indices()
        >>> idx
        array([[ 0,  1,  2,  3,  4, -1],
               [ 1,  2, -1,  3,  4,  5],
               [ 2,  3, -1,  4,  5,  6]])

        One can create a Biotite ``Alignment`` object from the transposed indices and
        the original sequences.

        >>> from biotite.sequence import NucleotideSequence  # doctest: +SKIP
        >>> from biotite.sequence.align import Alignment  # doctest: +SKIP
        >>> seqs = [NucleotideSequence("ACGTGA"),  # doctest: +SKIP
        ...         NucleotideSequence("TACTCA"),  # doctest: +SKIP
        ...         NucleotideSequence("GGACTGA")]  # doctest: +SKIP
        >>> aln = Alignment(seqs, idx.T)  # doctest: +SKIP
        >>> print(aln)  # doctest: +SKIP
        ACGTG-
        AC-TCA
        AC-TGA

        """
        errmsg = "Gap must be an integer, 'del', or 'mask'."
        valid_gaps = {"del", "mask"}
        if isinstance(gap, str):
            if gap not in valid_gaps:
                raise TypeError(errmsg)
        elif not np.issubdtype(type(gap), np.integer):
            raise TypeError(errmsg)

        bits = np.squeeze(self._to_bits())
        # TODO: Consider optimization using np.arange.
        # thought: initiate [-1, -1, -1 ... -1], then add slices of arange into it
        pos = np.repeat(1 - bits, self._lengths, axis=1)
        idx = np.cumsum(pos, axis=1, dtype=int) - 1
        if (starts := self._ranges[:, 0]).any():
            idx += starts.reshape(-1, 1)
        if gap == "del":
            keep = np.repeat(self._states == 0, self._lengths)
            return idx[:, keep]
        elif gap == "mask":
            return np.ma.array(idx, mask=(1 - pos))
        else:
            return np.where(pos, idx, gap)

    @classonlymethod
    def from_indices(cls, indices, gap=-1) -> "AlignPath":
        r"""Create an alignment path from character indices in the original sequences.

        Parameters
        ----------
        indices : array_like of int of shape (n_sequences, n_positions)
            Each element in the array is the index in the corresponding sequence.
        gap : int or "mask", optional
            The value which represents a gap in the alignment. Defaults to -1, but
            can be other values. If "mask", ``indices`` must be an
            ``np.ma.MaskedArray``. Cannot use "del".

        Returns
        -------
        AlignPath
            The alignment path created from the given indices.

        See Also
        --------
        to_indices

        Notes
        -----
        If a sequence in the alignment consists of entirely gap characters, its start
        position will be equal to the gap character.

        The input is equivalent to the transpose of the underlying data structure of
        Biotite's ``Alignment`` class [1]_.

        References
        ----------
        .. [1] `biotite.sequence.align.Alignment <https://www.biotite-python.org/
           latest/apidoc/biotite.sequence.align.Alignment.html>`_

        Examples
        --------
        >>> import numpy as np
        >>> from skbio.alignment import AlignPath
        >>> idx = np.array([[0, -1, -1,  1,  2,  3],
        ...                 [0,  1,  2, -1, -1, -1],
        ...                 [0, -1, -1,  1,  2, -1]])
        >>> path = AlignPath.from_indices(idx)
        >>> path
        <AlignPath, sequences: 3, positions: 6, segments: 4>

        One can convert a Biotite's ``Alignment`` object into a scikit-bio alignment
        path using this method. For example:

        >>> from biotite.sequence import NucleotideSequence  # doctest: +SKIP
        >>> from biotite.sequence.align import SubstitutionMatrix  # doctest: +SKIP
        >>> from biotite.sequence.align import align_optimal  # doctest: +SKIP
        >>> submat = SubstitutionMatrix.std_nucleotide_matrix()  # doctest: +SKIP
        >>> seq1 = NucleotideSequence("GATCGTC")  # doctest: +SKIP
        >>> seq2 = NucleotideSequence("ATCGCTC")  # doctest: +SKIP
        >>> res = align_optimal(seq1, seq2, submat)  # doctest: +SKIP
        >>> print(res[0])  # doctest: +SKIP
        GATCG-TC
        -ATCGCTC

        >>> trace = res[0].trace  # doctest: +SKIP
        >>> trace  # doctest: +SKIP
        array([[ 0, -1],
               [ 1,  0],
               [ 2,  1],
               [ 3,  2],
               [ 4,  3],
               [-1,  4],
               [ 5,  5],
               [ 6,  6]])

        >>> from skbio.alignment import PairAlignPath  # doctest: +SKIP
        >>> path = PairAlignPath.from_indices(trace.T)  # doctest: +SKIP
        >>> path  # doctest: +SKIP
        <PairAlignPath, positions: 8, segments: 4, CIGAR: '1D4M1I2M'>

        """
        if gap == "mask":
            return cls.from_bits(
                np.ma.getmask(indices),
                indices[
                    np.arange(indices.shape[0]),
                    np.argmax(indices != indices.fill_value, axis=1),
                ],
            )
        else:
            if isinstance(indices, np.ma.MaskedArray):
                raise TypeError("For masked arrays, gap must be 'mask'.")
            indices = np.asarray(indices)
            return cls.from_bits(
                indices == gap,
                indices[np.arange(indices.shape[0]), np.argmax(indices != gap, axis=1)],
            )
        # TODO
        # n. optimization

    def to_coordinates(self):
        r"""Generate an array of segment coordinates in the original sequences.

        Returns
        -------
        ndarray of int of shape (n_sequences, n_segments)
            Array where each value defines the start positions (index) of each segment
            for each sequence.

        See Also
        --------
        from_coordinates

        Notes
        -----
        The output is consistent with the underlying data structure of BioPython's
        ``Alignment`` class [1]_. Therefore, one can convert scikit-bio alignments into
        Biopython alignments, and vice versa.

        References
        ----------
        .. [1] `Bio.Align.Alignment <https://biopython.org/docs/latest/api/
           Bio.Align.html#Bio.Align.Alignment>`_

        Examples
        --------
        >>> from skbio.alignment import AlignPath
        >>> path = AlignPath(lengths=[2, 1, 2, 1],
        ...                  states=[0, 6, 0, 1],
        ...                  starts=[0, 1, 2])
        >>> coords = path.to_coordinates()
        >>> coords # doctest: +ELLIPSIS
        array([[0, 2, 3, 5, 5],
               [1, 3, 3, 5, 6],
               [2, 4, 4, 6, 7]]...

        One can create a Biopython ``Alignment`` object from the coordinates and the
        original sequences.

        >>> from Bio.Align import Alignment  # doctest: +SKIP
        >>> seqs = ["ACGTGA", "TACTCA", "GGACTGA"]  # doctest: +SKIP
        >>> aln = Alignment(seqs, coords)  # doctest: +SKIP
        >>> aln.coordinates is coords  # doctest: +SKIP
        True

        >>> aln.counts()  # doctest: +SKIP
        AlignmentCounts(gaps=5, identities=11, mismatches=2)

        >>> print(aln)  # doctest: +SKIP
                          0 ACGTG- 5
                          1 AC-TCA 6
                          2 AC-TGA 7

        """
        lens = self._lengths * (1 - self._to_bits())
        col0 = np.zeros((self._shape[0], 1), dtype=int)
        lens = np.append(col0, lens, axis=1)
        if self.starts.any():
            return lens.cumsum(axis=1) + self.starts.reshape(-1, 1)
        else:
            return lens.cumsum(axis=1)

    @classonlymethod
    def from_coordinates(cls, coords) -> "AlignPath":
        r"""Create an alignment path from an array of segment coordinates.

        Parameters
        ----------
        coords : array_like of int of shape (n_sequences, n_segments)
            Array where each value defines the start positions (index) of each segment
            for each sequence.

        Returns
        -------
        AlignPath
            The alignment path created from the given coordinates.

        See Also
        --------
        to_coordinates

        Notes
        -----
        The input is compatible with the underlying data structure of BioPython's
        ``Alignment`` class [1]_.

        References
        ----------
        .. [1] `Bio.Align.Alignment <https://biopython.org/docs/latest/api/
           Bio.Align.html#Bio.Align.Alignment>`_

        Examples
        --------
        >>> import numpy as np
        >>> from skbio.alignment import AlignPath
        >>> coordinates = np.array([[0, 1, 1, 3, 4],
        ...                         [0, 1, 3, 3, 3],
        ...                         [0, 1, 1, 3, 3]])
        >>> path = AlignPath.from_coordinates(coordinates)
        >>> path
        <AlignPath, sequences: 3, positions: 6, segments: 4>

        One can convert a Biopython's ``Alignment`` object into a scikit-bio alignment
        path using this method.

        >>> from Bio import Align  # doctest: +SKIP
        >>> a = Align.PairwiseAligner()  # doctest: +SKIP
        >>> res = a.align("GATCGTC", "ATCGCTC")  # doctest: +SKIP
        >>> print(res[0])  # doctest: +SKIP
        target            0 GATCG-TC 7
                          0 -||||-|| 8
        query             0 -ATCGCTC 7

        >>> coords = res[0].coordinates  # doctest: +SKIP
        >>> coords  # doctest: +SKIP
        array([[0, 1, 5, 5, 7],
               [0, 0, 4, 5, 7]])

        >>> from skbio.alignment import PairAlignPath  # doctest: +SKIP
        >>> path = PairAlignPath.from_coordinates(coords)  # doctest: +SKIP
        >>> path  # doctest: +SKIP
        <PairAlignPath, positions: 8, segments: 4, CIGAR: '1D4M1I2M'>

        """
        starts = coords[:, 0]
        diff = np.diff(coords)
        bits = diff == 0
        lens = diff[bits.argmin(axis=0), np.arange(diff.shape[1])]
        ints = np.packbits(bits, axis=0, bitorder="little")
        if ints.shape[0] == 1:
            ints = np.squeeze(ints, axis=0)
        return cls(lens, ints, starts=starts)


class PairAlignPath(AlignPath):
    r"""Store a pairwise alignment path between two sequences.

    ``PairAlignPath`` is a subclass of ``AlignPath``, with additional methods specific
    to pairwise alignments, such as the processing of CIGAR strings.

    Parameters
    ----------
    lengths : array_like of int of shape (n_segments,)
        Length of each segment in the alignment.
    states : array_like of uint8 of shape (n_segments,)
        Bits representing character (0) or gap (1) status per sequence per segment in
        the alignment.
    ranges : array_like of int of shape (n_sequences, 2), optional
        Start and stop positions of each sequence in the alignment.
    starts : array_like of int of shape (n_sequences,), optional
        Start position of each sequence in the alignment.
    stops : array_like of int of shape (n_sequences,), optional
        Stop position of each sequence in the alignment.

        .. note::
           If none of ``ranges``, ``starts`` or ``stops`` are provided,
           ``starts=[0, 0]`` will be used.

    See Also
    --------
    AlignPath
    skbio.sequence.Sequence
    skbio.alignment.TabularMSA

    Notes
    -----
    ``PairAlignPath`` uses a compact data structure to store alignment operations.
    Specifically, it encodes gap status in the two sequences in ``states``, a 2-D array
    with just one row of packed bits. The elements may be:

    * 0: Gap in neither sequence.
    * 1: Gap in sequence 1.
    * 2: Gap in sequence 2.
    * 3: Gap in both sequences.

    Meanwhile, it stores the length of segment per gap status in a 1-D array
    ``lengths``. For example, the following alignment:

    .. code-block:: none

       GAGCCAT-AC
       GC--CATAAC

    Can be represented by:

    .. code-block:: none

       lengths: 2 2 3 1 2
        states: 0 2 0 1 0

    This data structure resembles the :wiki:`CIGAR string <Sequence_alignment#CIGAR_
    Format>`, as defined in the SAM format specification [1]_. One can convert a
    pairwise alignment path to/from a CIGAR string using the :meth:`to_cigar` /
    :meth:`from_cigar` methods.

    The translation from CIGAR codes to ``states`` elements is as follows:

    +-----+----+------+----------------------------------+
    |Code |BAM |State |Description                       |
    +=====+====+======+==================================+
    |``M``|0   |0     |Alignment match                   |
    +-----+----+------+----------------------------------+
    |``I``|1   |1     |Insertion to the reference        |
    +-----+----+------+----------------------------------+
    |``D``|2   |2     |Deletion from the reference       |
    +-----+----+------+----------------------------------+
    |``N``|3   |2     |Skipped region from the reference |
    +-----+----+------+----------------------------------+
    |``S``|4   |1     |Soft clipping                     |
    +-----+----+------+----------------------------------+
    |``H``|5   |3     |Hard clipping                     |
    +-----+----+------+----------------------------------+
    |``P``|6   |3     |Padding                           |
    +-----+----+------+----------------------------------+
    |``=``|7   |0     |Sequence match                    |
    +-----+----+------+----------------------------------+
    |``X``|8   |0     |Sequence mismatch                 |
    +-----+----+------+----------------------------------+

    .. note::

       Sequences 1 and 2 are referred to as "query" and "reference" in the SAM format.

    See also the superclass :class:`AlignPath`, a generalization of this data structure
    to an arbitrary number of sequences.

    References
    ----------
    .. [1] https://samtools.github.io/hts-specs/SAMv1.pdf

    Examples
    --------
    >>> from skbio.alignment import pair_align
    >>> seqs = 'GATCGTC', 'ATCGCTC'
    >>> path = pair_align(*seqs).paths[0]
    >>> path
    <PairAlignPath, positions: 8, segments: 4, CIGAR: '1D4M1I2M'>

    >>> path.to_cigar()
    '1D4M1I2M'

    >>> path.lengths
    array([1, 4, 1, 2])

    >>> path.states
    array([[2, 0, 1, 0]], dtype=uint8)

    >>> path.to_aligned(seqs)
    ['GATCG-TC', '-ATCGCTC']

    """

    def __init__(
        self,
        lengths: ArrayLike,
        states: ArrayLike,
        *,
        ranges: Optional[ArrayLike] = None,
        starts: Optional[ArrayLike] = None,
        stops: Optional[ArrayLike] = None,
    ):
        if ranges is None and starts is None and stops is None:
            starts = [0, 0]
        super().__init__(lengths, states, ranges=ranges, starts=starts, stops=stops)
        if (self._states[0] > 3).any():
            raise ValueError(
                "For pairwise alignment, `states` must only contain zeros, ones, "
                "twos, or threes."
            )
        if self._shape[0] != 2:
            raise ValueError(
                "A pairwise alignment must represent exactly two sequences, but "
                f"{self._shape[0]} were given."
            )

    def __str__(self) -> str:
        r"""Return string representation of this alignment path."""
        return self.__repr__()

    def __repr__(self) -> str:
        r"""Return summary of the alignment path."""
        repr_ = (
            f"{self.__class__.__name__}, "
            f"positions: {self._shape[1]}, "
            f"segments: {self._lengths.size}"
        )
        width = 66 - len(repr_)  # maximum line width: 79
        cigar = self.to_cigar()
        if len(cigar) > width:
            cigar = cigar[: max(width - 3, 0)] + "..."
        return f"<{repr_}, CIGAR: '{cigar}'>"

    @classonlymethod
    def from_bits(cls, bits, starts=None) -> "PairAlignPath":
        r"""Create a pairwise alignment path from a bit array.

        Refer to :meth:`AlignPath.from_bits` for usage.

        Parameters
        ----------
        bits : array_like of (0, 1) of shape (2, n_positions)
            Bit array representing the alignment.
        starts : array_like of int of shape (2,), optional
            Start positions of sequences.

        Returns
        -------
        PairAlignPath
            The pairwise alignment path created from the given bit array.

        See Also
        --------
        AlignPath.from_bits

        """
        bits = np.asarray(bits)
        if bits.shape[0] != 2:
            raise ValueError(
                "A pairwise alignment must represent exactly two sequences, but "
                f"{bits.shape[0]} were given."
            )
        return super().from_bits(bits, starts)

    def to_cigar(self, seqs: Optional[List[Union["Sequence", str]]] = None):
        r"""Generate a CIGAR string representing the pairwise alignment path.

        Parameters
        ----------
        seqs : list of skbio.Sequence or string
            A pair of sequences to generate CIGAR string. If provided, will
            distinguish match (``=``) and mismatch (``X``). Otherwise, will uniformly
            note them as (mis)match (``M``). The first sequence in the list is the
            query sequence, the second is the reference sequence.

        Returns
        -------
        str
            CIGAR string representing the pairwise alignment path.

        Examples
        --------
        >>> from skbio.alignment import PairAlignPath
        >>> path = PairAlignPath(lengths=[2, 5, 3, 1],
        ...                      states=[0, 3, 2, 1],
        ...                      starts=[0, 0])
        >>> path.to_cigar()
        '2M5P3D1I'

        """
        from skbio.sequence import Sequence

        cigar = []
        states = np.squeeze(self._states)

        # TODO: Make this compatible with `SequenceLike`.
        if seqs is not None:
            # test if seqs is strings or Sequence object or something else
            if isinstance(seqs[0], str) and isinstance(seqs[1], str):
                seq1 = np.frombuffer(seqs[0].encode("ascii"), dtype=np.uint8)
                seq2 = np.frombuffer(seqs[1].encode("ascii"), dtype=np.uint8)
            elif isinstance(seqs[0], Sequence) and isinstance(seqs[1], Sequence):
                seq1 = seqs[0]._bytes
                seq2 = seqs[1]._bytes
            else:
                raise TypeError("`seqs` must be strings or Sequence objects.")

            idx1, idx2 = self._ranges[:, 0]

            for length, state in zip(self._lengths, states):
                if state == 0:
                    match_arr = seq1[idx1 : idx1 + length] == seq2[idx2 : idx2 + length]
                    char_arr = np.where(match_arr, "=", "X")
                    n = len(char_arr)
                    count = 1
                    curr_char = char_arr[0]
                    for i in range(1, n):
                        if char_arr[i] == curr_char:
                            count += 1
                        else:
                            cigar.append(str(count) + curr_char)
                            curr_char = char_arr[i]
                            count = 1
                    cigar.append(str(count) + curr_char)
                    idx1 += length
                    idx2 += length
                elif state == 1:
                    cigar.append(str(length) + "I")
                    idx2 += length
                elif state == 2:
                    cigar.append(str(length) + "D")
                    idx1 += length
                else:  # state == 3:
                    cigar.append(str(length) + "P")
            return "".join(cigar)
        else:
            return "".join(
                f"{L}{C}" for L, C in zip(self._lengths, _cigar_codes[states])
            )

    @classonlymethod
    def from_cigar(
        cls, cigar: Union[str, bytes], starts: Optional[ArrayLike] = None
    ) -> "PairAlignPath":
        r"""Create a pairwise alignment path from a CIGAR string.

        Parameters
        ----------
        cigar : str or bytes
            CIGAR format string used to build the PairAlignPath.
        starts : array_like of (int, int), optional
            Start position of each sequence in the alignment. If omitted, will set as
            zeros.

        Returns
        -------
        PairAlignPath
            The pairwise alignment path created from the given CIGAR string.

        Raises
        ------
        ValueError
            CIGAR string is empty.
        ValueError
            CIGAR string contains invalid characters.

        Examples
        --------
        >>> from skbio.alignment import PairAlignPath
        >>> cigar = "2M5P3D1I"
        >>> path = PairAlignPath.from_cigar(cigar)
        >>> path
        <PairAlignPath, positions: 11, segments: 4, CIGAR: '2M5P3D1I'>

        """
        if not cigar:
            raise ValueError("CIGAR string must not be empty.")
        if isinstance(cigar, str):
            try:
                cigar = cigar.encode("ascii")
            except UnicodeEncodeError:
                raise ValueError("CIGAR string contains invalid character(s).")
        lens, gaps = _parse_cigar(cigar)
        return cls(lens, gaps, starts=(starts or [0, 0]))


def _parse_cigar(cigar: bytes) -> Tuple[NDArray, NDArray]:
    r"""Convert a CIGAR string into run-length vectors.

    Parameters
    ----------
    cigar : bytes
        CIGAR string.

    Returns
    -------
    ndarray of intp
        Lengths.
    ndarray of uint8
        Gap states.

    """
    lens, gaps = [], []
    L, no_ones, last = 0, 1, None
    for char in cigar:
        if char in _cigar_mapping:  # letter (operation)
            curr = _cigar_mapping[char]
            if curr == last:
                lens[-1] += L + no_ones
            else:
                lens.append(L + no_ones)
                gaps.append(curr)
                last = curr
            L, no_ones = 0, 1
        elif 48 <= char <= 57:  # number (length)
            no_ones = 0
            L = L * 10 + char - 48
        else:
            raise ValueError("CIGAR string contains invalid character(s).")
    return (
        np.array(lens, dtype=np.intp),
        np.array(gaps, dtype=np.uint8),
    )


def _merge_same_1d(lens, gaps):
    r"""Merge consecutive same gap states and sum corresponding lengths.

    Parameters
    ----------
    lens : ndarray of int of shape (n_segments,)
        Segment lengths.
    gaps : ndarray of uint8 of shape (n_segments,)
        Segment states.

    Returns
    -------
    ndarray of int of shape (n_merged_segments,)
        Merged segment lengths.
    ndarray of uint8 of shape (n_merged_segments,)
        Merged segment states.

    """
    idx = np.flatnonzero(np.r_[True, gaps[1:] != gaps[:-1]])
    return np.add.reduceat(lens, idx, dtype=lens.dtype), gaps[idx]


def _merge_same_2d(lens, gaps):
    r"""Merge consecutive same gap states and sum corresponding lengths.

    Similar to ``_merge_same_1d`` but applies to a 2-D gaps array.

    """
    idx = np.flatnonzero(np.r_[True, (gaps[:, 1:] != gaps[:, :-1]).any(axis=0)])
    return np.add.reduceat(lens, idx, dtype=lens.dtype), gaps[:, idx]


def _run_length_encode(s):
    r"""Perform run length encoding on a string.

    Parameters
    ----------
    s : str
        String on which to perform run length encoding.
    """
    input_arr = np.array(list(s))
    idx = np.append(0, np.where(input_arr[:-1] != input_arr[1:])[0] + 1)
    count = np.diff(np.concatenate((idx, [len(s)])))
    unique = input_arr[idx]
    return "".join(str(c) + u for c, u in zip(count, unique))
