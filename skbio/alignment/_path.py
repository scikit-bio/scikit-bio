# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections

import numpy as np

from skbio._base import SkbioObject
from skbio.util._decorator import classonlymethod
from skbio.sequence import Sequence


_Shape = collections.namedtuple("Shape", ["sequence", "position"])

# CIGAR codes indexed by states in PairAlignPath
_cigar_codes = np.array(["M", "I", "D", "P"])

# Mapping of CIGAR codes to states in PairAlignPath
_cigar_mapping = {
    "M": 0,
    "I": 1,
    "D": 2,
    "P": 3,
    "=": 0,
    "X": 0,
    "N": 2,
    "S": 1,
    "H": 3,
}


class AlignPath(SkbioObject):
    r"""Store an alignment path between sequences.

    It defines the operations of aligning sequences by inserting gaps into designated
    positions between characters.

    The underlying data structure of the ``AlignPath`` class efficiently represents a
    sequence alignment as two equal-length vectors: lengths and gap status. The lengths
    vector contains the lengths of individual segments of the alignment with consistent
    gap status. The gap status vector contains the encoded bits of the gap (1) and
    character (0) status for each position in the alignment.

    This data structure is detached from the original sequences and is highly memory
    efficient. It permits fully vectorized operations and enables efficient conversion
    between various formats such as CIGAR, tabular, indices (Biotite), and coordinates
    (Biopython).

    Parameters
    ----------
    lengths : array_like of int of shape (n_segments,)
        Length of each segment in the alignment.
    states : array_like of uint8 of shape (n_segments,) or (n_packs, n_segments)
        Packed bits representing character (0) or gap (1) status per sequence per
        segment in the alignment.
    starts : array_like of int of shape (n_sequences,)
        Start position (0-based) of each sequence in the alignment.

    See Also
    --------
    PairAlignPath
    skbio.sequence.Sequence
    skbio.alignment.TabularMSA

    Notes
    -----
    The underlying logic of the ``AlignPath`` data structure is rooted in two concepts:
    run length encoding and bit arrays.

    The lengths array is calculated by performing run length encoding on the alignment,
    considering each segment with consistent gap status to be an individual unit in the
    encoding. In the above example, the first three positions of the alignment contain
    no gaps, so the first value in the lengths array is 3, and so on.

    The states array is calculated by turning the alignment segments into a bit array
    where gaps become 1's, and characters become zeros. Then, the 0's and 1's are
    converted into bytes. In the above example, the fourth segment, which has length 1,
    would become [0, 1, 1], which then becomes 6.

    Examples
    --------
    Create an ``AlignPath`` object from a ``TabularMSA`` object with three DNA
    sequences and 20 positions.

    >>> from skbio import DNA, TabularMSA
    >>> from skbio.alignment import AlignPath
    >>> seqs = [
    ...    DNA('CGGTCGTAACGCGTA---CA'),
    ...    DNA('CAG--GTAAG-CATACCTCA'),
    ...    DNA('CGGTCGTCAC-TGTACACTA')
    ... ]
    >>> msa = TabularMSA(seqs)
    >>> msa
    TabularMSA[DNA]
    ----------------------
    Stats:
        sequence count: 3
        position count: 20
    ----------------------
    CGGTCGTAACGCGTA---CA
    CAG--GTAAG-CATACCTCA
    CGGTCGTCAC-TGTACACTA
    >>> path = AlignPath.from_tabular(msa)
    >>> path
    AlignPath
    Shape(sequence=3, position=20)
    lengths: [3 2 5 1 4 3 2]
    states: [0 2 0 6 0 1 0]

    """

    def __init__(self, lengths, states, starts):
        self._lengths = np.asarray(lengths, dtype=np.int64)
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

        # Number of sequences needs to be explicitly provided, because the packed bits
        # does not contain this information. (It is merely in multiples of 8.)
        self._starts = np.asarray(starts, dtype=np.int64)
        if self._starts.ndim > 1:
            raise TypeError("`starts` must be a 1-D array.")

        n_sequences = len(self._starts)
        if np.ceil(n_sequences / 8) != self._states.shape[0]:
            max_seqs = self._states.shape[0] * 8
            raise ValueError(
                f"Number of sequences in `starts` ({n_sequences}) and capacity of "
                f"`states` ({max_seqs - 7} to {max_seqs}) do not match."
            )

        # Shape is n_sequences (rows) x n_positions (columns), which is consistent with
        # TabularMSA
        n_positions = int(self._lengths.sum())
        self._shape = _Shape(sequence=n_sequences, position=n_positions)

    def __str__(self):
        r"""Return string representation of this alignment path."""
        # Not sure if this makes sense for this class, but it is needed for all
        # SkbioObjects.
        return self.__repr__()

    def __repr__(self):
        r"""Return summary of the alignment path."""
        return (
            f"{self.__class__.__name__}\n{self._shape}\nlengths: "
            f"{self._lengths}\nstates: {np.squeeze(self._states)}"
        )

    @property
    def lengths(self):
        """Array of lengths of segments in alignment path."""
        return self._lengths

    @property
    def states(self):
        """Array of gap status of segments in alignment path."""
        return self._states

    @property
    def starts(self):
        """Array of start positions of sequences in the alignment."""
        return self._starts

    @property
    def shape(self):
        """Number of sequences (rows) and positions (columns).

        Notes
        -----
        This property is not writeable.

        """
        return self._shape

    def to_bits(self, expand=True):
        r"""Unpack the alignment path into an array of bits.

        .. versionchanged:: 0.6.4
            The default behavior now returns positions.

        Parameters
        ----------
        expand : bool, optional
            If True (default), each column in the returned array represents a position
            in the original alignment. If False, each column represents a segment in
            the alignment path.

            .. versionadded:: 0.6.4

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
        bits = np.unpackbits(
            self._states, axis=0, count=self._shape[0], bitorder="little"
        )
        if expand:
            bits = np.repeat(bits, self._lengths, axis=1)
        return bits

    @classonlymethod
    def from_bits(cls, bits, starts=None):
        r"""Create an alignment path from a bit array (0 - character, 1 - gap).

        Parameters
        ----------
        bits : array_like of (0, 1) of shape (n_sequences, n_positions)
            Array of zeros (character) and ones (gap) which represent the alignment.
        starts : array_like of int of shape (n_sequences,), optional
            Start position (0-based) of each sequence in the alignment. If omitted,
            will set as zeros.

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
        AlignPath
        Shape(sequence=3, position=8)
        lengths: [2 2 2 1 1]
        states: [0 2 0 6 0]

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
        return cls(lens, ints, starts)

    @classonlymethod
    def from_tabular(cls, msa):
        r"""Create an alignment path from a `TabularMSA` object.

        Parameters
        ----------
        msa : TabularMSA
            TabularMSA to be converted into AlignPath object.

        Returns
        -------
        AlignPath
            The alignment path created from the TabularMSA object.

        Notes
        -----
        The returned alignment path will span across the entire tabular MSA. Its start
        positions will be uniformly zeros.

        See Also
        --------
        skbio.TabularMSA.from_path_seqs

        """
        # Convert TabularMSA into a 2D array of bytes.
        # TODO: TabularMSA itself should be refactored to have this as the default data
        # structure.
        byte_arr = np.stack([x._bytes for x in msa._seqs])

        # Get gap character code.
        gap_chars = [ord(x) for x in msa.dtype.gap_chars]

        # Identify gap positions, and convert them into a bit array, then create an
        # alignment path based on it.
        return cls.from_bits(np.isin(byte_arr, gap_chars))

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

        Examples
        --------
        >>> from skbio.alignment import AlignPath
        >>> path = AlignPath(lengths=[1, 2, 2, 1],
        ...                  states=[0, 5, 2, 6],
        ...                  starts=[0, 0, 0])
        >>> path.to_indices()
        array([[ 0, -1, -1,  1,  2,  3],
               [ 0,  1,  2, -1, -1, -1],
               [ 0, -1, -1,  1,  2, -1]])

        """
        errmsg = "Gap must be an integer, np.nan, np.inf, 'del', or 'mask'."
        valid_gaps = {"del", "mask"}
        if isinstance(gap, str):
            if gap not in valid_gaps:
                raise TypeError(errmsg)
        elif isinstance(gap, float):
            if not (np.isnan(gap) or np.isinf(gap)):
                raise TypeError(errmsg)
        elif not np.issubdtype(type(gap), np.integer):
            raise TypeError(errmsg)

        bits = np.squeeze(self.to_bits(expand=False))
        # TODO: Consider optimization using np.arange.
        # thought: initiate [-1, -1, -1 ... -1], then add slices of arange into it
        pos = np.repeat(1 - bits, self._lengths, axis=1)
        idx = np.cumsum(pos, axis=1, dtype=int) - 1
        if self._starts.any():
            idx += self._starts.reshape(-1, 1)
        if gap == "del":
            keep = np.repeat(self._states == 0, self._lengths)
            return idx[:, keep]
        elif gap == "mask":
            return np.ma.array(idx, mask=(1 - pos))
        else:
            return np.where(pos, idx, gap)

    @classonlymethod
    def from_indices(cls, indices, gap=-1):
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

        Notes
        -----
        If a sequence in the alignment consists of entirely gap characters, its start
        position will be equal to the gap character.

        Examples
        --------
        >>> import numpy as np
        >>> from skbio.alignment import AlignPath
        >>> indices = np.array([[0, -1, -1,  1,  2,  3],
        ...                     [0,  1,  2, -1, -1, -1],
        ...                     [0, -1, -1,  1,  2, -1]])
        >>> path = AlignPath.from_indices(indices)
        >>> path
        AlignPath
        Shape(sequence=3, position=6)
        lengths: [1 2 2 1]
        states: [0 5 2 6]

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

        Examples
        --------
        >>> from skbio.alignment import AlignPath
        >>> path = AlignPath(lengths=[1, 2, 2, 1],
        ...                  states=[0, 5, 2, 6],
        ...                  starts=[0, 0, 0])
        >>> path.to_coordinates() # doctest: +ELLIPSIS
        array([[0, 1, 1, 3, 4],
               [0, 1, 3, 3, 3],
               [0, 1, 1, 3, 3]]...

        """
        lens = self._lengths * (1 - self.to_bits(expand=False))
        col0 = np.zeros((self._shape[0], 1), dtype=int)
        lens = np.append(col0, lens, axis=1)
        if self.starts.any():
            return lens.cumsum(axis=1) + self.starts.reshape(-1, 1)
        else:
            return lens.cumsum(axis=1)

    @classonlymethod
    def from_coordinates(cls, coords):
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

        Examples
        --------
        >>> import numpy as np
        >>> from skbio.alignment import AlignPath
        >>> coordinates = np.array([[0, 1, 1, 3, 4],
        ...                         [0, 1, 3, 3, 3],
        ...                         [0, 1, 1, 3, 3]])
        >>> path = AlignPath.from_coordinates(coordinates)
        >>> path
        AlignPath
        Shape(sequence=3, position=6)
        lengths: [1 2 2 1]
        states: [0 5 2 6]

        """
        starts = coords[:, 0]
        diff = np.diff(coords)
        bits = diff == 0
        lens = diff[bits.argmin(axis=0), np.arange(diff.shape[1])]
        ints = np.packbits(bits, axis=0, bitorder="little")
        if ints.shape[0] == 1:
            ints = np.squeeze(ints, axis=0)
        return cls(lens, ints, starts)


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
    starts : array_like of (int, int), optional
        Start position (0-based) of each sequence in the alignment. Default is (0, 0).

    See Also
    --------
    AlignPath
    skbio.sequence.Sequence
    skbio.alignment.TabularMSA

    """

    def __init__(self, lengths, states, starts=(0, 0)):
        super().__init__(lengths, states, starts)
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

    def __str__(self):
        r"""Return string representation of this alignment path."""
        return self.__repr__()

    def __repr__(self):
        r"""Return summary of the alignment path."""
        return (
            f"<{self.__class__.__name__}, position count: {self._shape[1]}, "
            f"CIGAR: '{self.to_cigar()}'>"
        )

    @classonlymethod
    def from_bits(cls, bits, starts=None):
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

    def to_cigar(self, seqs=None):
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
        cigar = []

        states = np.squeeze(self._states)

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

            idx1, idx2 = self._starts

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
                elif state == 3:
                    cigar.append(str(length) + "P")
            return "".join(cigar)
        else:
            return "".join(
                f"{L}{C}" for L, C in zip(self._lengths, _cigar_codes[states])
            )

    @classonlymethod
    def from_cigar(cls, cigar, starts=None):
        r"""Create a pairwise alignment path from a CIGAR string.

        Parameters
        ----------
        cigar : str
            CIGAR format string used to build the PairAlignPath.
        starts : array_like of (int, int), optional
            Start position (0-based) of each sequence in the alignment. If omitted,
            will set as zeros.

        Returns
        -------
        PairAlignPath
            The pairwise alignment path created from the given CIGAR string.

        Examples
        --------
        >>> from skbio.alignment import PairAlignPath
        >>> cigar = "2M5P3D1I"
        >>> path = PairAlignPath.from_cigar(cigar)
        >>> path
        <PairAlignPath, position count: 11, CIGAR: '2M5P3D1I'>

        """
        # Make sure cigar is not empty.
        if not cigar:
            raise ValueError("CIGAR string must not be empty.")

        lengths = []
        gaps = []
        current_length = 0
        no_ones = True
        for char in cigar:
            if char.isdigit():
                no_ones = False
                current_length = current_length * 10 + int(char)
            elif char in _cigar_mapping:
                if no_ones:
                    lengths.append(current_length + 1)
                else:
                    lengths.append(current_length)
                gaps.append(_cigar_mapping[char])
                current_length = 0
                no_ones = True
            else:
                raise ValueError("CIGAR string contains invalid character(s).")
        lengths, gaps = _fix_arrays(lengths=np.array(lengths), gaps=np.array(gaps))

        return cls(lengths, gaps, [0, 0] if starts is None else starts)


def _fix_arrays(lengths, gaps):
    r"""Merge consecutive same values from gaps array and sum corresponding values
    in lengths array.

    Parameters
    ----------
    lengths : array_like of int of shape (n_segments,)
        Length of each segment in the alignment.
    gaps : array_like of uint8 of shape (n_segments,) or (n_packs, n_segments)
        Packed bits representing character (0) or gap (1) status per sequence per
        segment in the alignment.
    """
    idx = np.diff(gaps, prepend=np.array([True])) != 0
    gaps_out = np.asarray(gaps[idx])
    groups = np.cumsum(idx)
    lengths_out = np.asarray(np.bincount(groups, weights=lengths).astype(int)[1:])

    return lengths_out, gaps_out


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
