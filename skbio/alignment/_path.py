# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio._base import SkbioObject
from skbio.util._decorator import classonlymethod
from skbio.sequence import Sequence


class AlignPath(SkbioObject):
    r"""Create an alignment path from segment lengths and states.

    Parameters
    ----------
    lengths : array_like of int of shape (n_segments,)
        Length of each segment in the alignment.
    states : array_like of unit8 of shape (n_segments,) or (n_packs, n_segments)
        Packed bits representing character (0) or gap (1) status per sequence per
        segment in the alignment.
    starts : array_like of int of shape (n_sequences,), optional
        Start position (0-based) of each sequence in the alignment.

    See Also
    --------
    skbio.sequence.Sequence
    skbio.alignment.TabularMSA

    """

    def __init__(self, lengths, states, starts):
        # Number of sequences needs to be explicitly provided, because the packed bits
        # does not contain this information. (It is merely in multiples of 8.)
        self.lengths = np.asarray(lengths, dtype=np.int64)
        self.n_positions = self.lengths.sum()
        self.states = np.atleast_2d(np.asarray(states, dtype=np.uint8))

        # start positions
        self.starts = np.asarray(starts, dtype=np.int64)
        if self.starts.ndim > 1:
            raise ValueError("`starts` must be a 1-D vector.")
        self.n_seqs = self.starts.size
        if np.ceil(self.n_seqs / 8) != self.states.shape[0]:
            raise ValueError("Sizes of `starts` and `states` do not match.")

        # Shape is n_seqs (rows) x n_positions (columns), which is consistent with
        # TabularMSA
        self.shape = (self.n_seqs, self.n_positions)

    def __str__(self):
        r"""Return string representation of this AlignPath."""
        # Not sure if this makes sense for this class, but it is needed for all
        # SkbioObjects.
        return self.__repr__()

    def __repr__(self):
        r"""Return summary of the alignment path."""
        return (
            f"<{self.__class__.__name__}, shape: {self.shape}, lengths: "
            f"{self.lengths}, states: {self.states}>"
        )

    def to_bits(self):
        r"""Unpack states into an array of bits."""
        return np.unpackbits(
            np.atleast_2d(self.states), axis=0, count=self.shape[0], bitorder="little"
        )

    @classonlymethod
    def from_bits(cls, bits):
        r"""Create an alignment path from a bit array (0 - character, 1 - gap).

        Parameters
        ----------
        bits : array_like of 0's and 1's of shape (n_seqs, n_positions)
            Array of zeros (character) and ones (gap) which represent the alignment.

        Returns
        -------
        AlignPath
            The alignment path created from the given bit array.
        """
        # Pack bits into integers.
        ints = np.packbits(bits, axis=0, bitorder="little")

        # If there are 8 or less sequences, squeeze the 2D array into 1D.
        # This is an optional optimization especially for pairwise alignment.
        if ints.shape[0] == 1:
            ints = np.squeeze(ints, axis=0)

            # Get indices where segments start. Equivalent to but faster than:
            # idx = np.where(np.diff(ints, prepend=np.nan))[0]
            idx = np.append(0, np.where(ints[:-1] != ints[1:])[0] + 1)

            # Get lengths of segments. Equivalent to but faster than:
            # lens = np.diff(idx, append=ints.size)
            lens = np.append(idx[1:] - idx[:-1], ints.size - idx[-1])

            # Keep indices of segment starts.
            ints = ints[idx]

        # This is the 2D equivalent of the above code.
        else:
            idx = np.append(
                0, np.where((ints[:, :-1] != ints[:, 1:]).any(axis=0))[0] + 1
            )
            lens = np.append(idx[1:] - idx[:-1], ints.shape[1] - idx[-1])
            ints = ints[:, idx]

        # return per-segment lengths and bits
        return cls(lens, ints, np.zeros(bits.shape[0], dtype=int))

    @classonlymethod
    def from_tabular(cls, msa):
        r"""Create an alignment path from a `TabularMSA` object.

        Parameters
        ----------
        msa : TabularMSA object
            TabularMSA to be converted into AlignPath object.

        Returns
        -------
        AlignPath
            The alignment path created from the TabularMSA object.
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
            If numeric, replace gaps with this value. If "del", delete columns
            that have any gap. If "mask", mask gaps. Default is -1.
        """
        valid_gaps = {"del", "mask"}
        if isinstance(gap, str):
            if gap not in valid_gaps:
                raise ValueError(
                    "Gap must be an integer, np.nan, np.inf, 'del', " "or 'mask'."
                )
        elif isinstance(gap, float):
            if not (np.isnan(gap) or np.isinf(gap)):
                raise ValueError(
                    "Gap must be an integer, np.nan, np.inf, 'del', " "or 'mask'."
                )
        elif not np.issubdtype(type(gap), np.integer):
            raise ValueError(
                "Gap must be an integer, np.nan, np.inf, 'del', or " "'mask'."
            )

        bits = np.squeeze(self.to_bits())
        # TODO: Consider optimization using np.arange.
        # thought: initiate [-1, -1, -1 ... -1], then add slices of arange into it
        pos = np.repeat(1 - bits, self.lengths, axis=1)
        idx = np.cumsum(pos, axis=1, dtype=int) - 1
        if gap == "del":
            keep = np.repeat(self.states == 0, self.lengths)
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
        indices : array_like of int of shape (n_seqs, n_positions)
            Each element in the array is the index in the corresponding sequence.
        gap : int or "mask", optional
            The value which represents a gap in the alignment. Defaults to -1, but
            can be other values. If "mask", `indices` must be a masked array. Cannot
            use "del".

        Returns
        -------
        AlignPath
            The alignment path created from the given indices.
        """
        if gap == "mask":
            return cls.from_bits(np.ma.getmask(indices))
        else:
            return cls.from_bits(indices == gap)

    def to_coordinates(self):
        r"""Generate an array of segment coordinates in the original sequences."""
        lens = self.lengths * (1 - self.to_bits())
        col0 = np.zeros((self.shape[0], 1), dtype=int)
        lens = np.append(col0, lens, axis=1)
        return lens.cumsum(axis=1)

    @classonlymethod
    def from_coordinates(cls, coords):
        r"""Generate the an alignment path from an array of segment coordinates.

        Parameters
        ----------
        coords : array_like of int of shape (n_seqs, n_segments)
            Array where each value defines the start positions (index) of each segment
            for each sequence.

        Returns
        -------
        AlignPath
            The alignment path created from the given coordinates.
        """
        diff = np.diff(coords)
        bits = diff == 0
        lens = diff[bits.argmin(axis=0), np.arange(diff.shape[1])]
        ints = np.packbits(bits, axis=0, bitorder="little")
        if ints.shape[0] == 1:
            ints = np.squeeze(ints, axis=0)
        return cls(lens, ints, np.zeros(diff.shape[0], dtype=int))


class PairAlignPath(AlignPath):
    r"""Create a pairwise alignment path from segment lengths and states.

    Parameters
    ----------
    lengths : array_like of int of shape (n_segments,)
        Length of each segment in the alignment.
    states : array_like of unit8 of shape (n_segments,) or (n_packs, n_segments)
        Packed bits representing character (0) or gap (1) status per sequence per
        segment in the alignment.
    starts : array_like of int of shape (n_sequences,), optional
        Start position (0-based) of each sequence in the alignment.

    See Also
    --------
    skbio.sequence.Sequence
    skbio.alignment.TabularMSA

    """

    def __str__(self):
        r"""Return string representation of this AlignPath."""
        return self.__repr__()

    def __repr__(self):
        r"""Return summary of the alignment path."""
        return (
            f"<{self.__class__.__name__}, shape: {self.shape}, "
            f"CIGAR: '{self.to_cigar()}'>"
        )

    @classonlymethod
    def from_bits(cls, bits):
        r"""Create an alignment path from a bit array.

        Parameters
        ----------
        bits : array_like of 0's and 1's of shape (n_seqs, n_positions)
            Array of zeros (character) and ones (gap) which represent the alignment.

        Returns
        -------
        PairAlignPath
            The alignment path of the provided bit array.
        """
        # Ensure bits is a 2D array-like of ones and zeros.
        if not isinstance(bits, np.ndarray):
            bits = np.atleast_2d(bits)
        if bits.ndim != 2 or bits.shape[1] == 0:
            raise ValueError("Input 'bits' must be a non-empty 2D numpy array.")
        if not (np.logical_or(bits == 0, bits == 1).all()):
            raise ValueError("Input 'bits' must contain only zeros and ones.")

        ints = bits[0] + (bits[1] << 1)
        idx = np.append(0, np.where(ints[:-1] != ints[1:])[0] + 1)
        lens = np.append(idx[1:] - idx[:-1], ints.size - idx[-1])
        return cls(lens, ints[idx], np.zeros(bits.shape[0], dtype=int))

    def to_bits(self):
        r"""Unpack states into an array of bits."""
        if not np.all(np.isin(self.states, [0, 1, 2, 3])):
            raise ValueError(
                "For pairwise alignment, `states` must only contain "
                "zeros, ones, twos, or threes."
            )
        return np.stack([self.states & 1, self.states >> 1])

    def to_cigar(self, seqs=None):
        r"""Generate a CIGAR string representing the pairwise alignment path.

        Parameters
        ----------
        seqs : list of skbio.Sequence or string
            A pair of sequences to generate CIGAR string. If provided, will
            distinguish match (``=``) and mismatch (``X``). Otherwise, will uniformly
            note them as (mis)match (``M``). The first sequence in the list is the
            query sequence, the second is the reference sequence.
        """
        cigar = []

        states = np.squeeze(self.states)

        if seqs is not None:
            # test if seqs is strings or Sequence object or something else
            if isinstance(seqs[0], str) and isinstance(seqs[1], str):
                seq1 = np.frombuffer(seqs[0].encode("ascii"), dtype=np.uint8)
                seq2 = np.frombuffer(seqs[1].encode("ascii"), dtype=np.uint8)
            elif isinstance(seqs[0], Sequence) and isinstance(seqs[1], Sequence):
                seq1 = seqs[0]._bytes
                seq2 = seqs[1]._bytes
            else:
                raise ValueError("`seqs` must be strings or Sequence objects.")

            idx1, idx2 = self.starts

            for length, state in zip(self.lengths, states):
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
                f"{L}{C}" for L, C in zip(self.lengths, _cigar_codes[states])
            )

    @classonlymethod
    def from_cigar(cls, cigar):
        r"""Create a pairwise alignment path from a CIGAR string.

        Parameters
        ----------
        cigar : str
            CIGAR format string used to build the PairAlignPath.

        Returns
        -------
        AlignPath
            The alignment path created from the given CIGAR string.
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
                raise ValueError("Invalid characters in CIGAR string.")
        lengths, gaps = cls._fix_arrays(lengths=np.array(lengths), gaps=np.array(gaps))
        return cls(lengths, gaps, [0, 0])

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


def _run_length_encode(string_in):
    r"""Perform run length encoding on a string.

    Parameters
    ----------
    string_in : str
        String on which to perform run length encoding.
    """
    input_arr = np.array(list(string_in))
    idx = np.append(0, np.where(input_arr[:-1] != input_arr[1:])[0] + 1)
    count = np.diff(np.concatenate((idx, [len(string_in)])))
    unique = input_arr[idx]
    encoded_str = "".join(str(c) + u for c, u in zip(count, unique))

    return encoded_str


_cigar_codes = np.array(["M", "I", "D", "P"])

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
