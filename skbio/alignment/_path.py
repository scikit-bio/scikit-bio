# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio._base import SkbioObject
from skbio.util._decorator import overrides, classonlymethod


class AlignPath(SkbioObject):
    def __init__(self, lengths, states, starts):
        """Create an alignment path from segment lengths and states.

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
        # Number of sequences needs to be explicitly provided, because the packed bits
        # does not contain this information. (It is merely in multiples of 8.)
        self.lengths = np.asarray(lengths, dtype=int)
        n_positions = self.lengths.sum()
        self.states = np.atleast_2d(np.asarray(states, dtype=np.uint8))

        # start positions
        self.starts = np.asarray(starts, dtype=int)
        if self.starts.ndim > 1:
            raise ValueError("`starts` must be a 1-D vector.")
        n_seqs = self.starts.size
        if np.ceil(n_seqs / 8) != self.states.shape[0]:
            raise ValueError("Sizes of `starts` and `states` do not match.")

        # Shape is n_seqs (rows) x n_positions (columns), which is consistent with
        # TabularMSA
        self.shape = (n_seqs, n_positions)

        # TODO: An additional parameter `starts` should record the starting positions
        # of the alignment in each sequence (shape: (n_seqs,)). This is important
        # especially for local alignments (e.g., search for a short read in a genome
        # sequence). If it is provided, `n_seqs` is no longer necessary.
        # self.starts = starts

        # TODO: Needs to think about whether reverse complemented (Boolean array of
        # (n_seqs,)) should be included as a parameter. It is only relevant for
        # nucleotide sequences.

    def __str__(self):
        """Return string representation of this AlignPath."""
        # Not sure if this makes sense for this class, but it is needed for all
        # SkbioObjects.
        return self.__repr__()

    def __repr__(self):
        """Return summary of the alignment path."""
        return (
            f"<{self.__class__.__name__}, shape: {self.shape}, lengths: "
            f"{self.lengths}, states: {self.states}"
        )

    def subset(self):
        """Select subset of sequences from AlignPath.

        Better to have ability to index an AlignPath for particular sequences,
        then compress these sequences using something like _fix_arrays"""
        pass

    def to_bits(self):
        """Unpack states into an array of bits."""
        return np.unpackbits(
            np.atleast_2d(self.states), axis=0, count=self.shape[0], bitorder="little"
        )

    @classonlymethod
    def from_bits(cls, bits):
        """Create an alignment path from a bit array (0 - char, 1 - gap)."""
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
        """Create an alignment path from a `TabularMSA` object."""
        # Convert TabularMSA in to a 2D array of bytes.
        # TODO: TabularMSA itself should be refactored to have this as the default data
        # structure.
        byte_arr = np.stack([x._bytes for x in msa._seqs])

        # Get gap character code.
        gap_chars = [ord(x) for x in msa.dtype.gap_chars]

        # Identify gap positions, and convert them into a bit array, then create an
        # alignment path based on it.
        return cls.from_bits(np.isin(byte_arr, gap_chars))

    def to_indices(self, gap=-1):
        """Generate an array of indices of characters in the original sequences.

        gap = "del": delete columns that have any gap, "mask": mask gaps, others
        (default: -1):
        fill gaps with this value.
        """
        valid_gaps = {-1, "del", "mask"}
        if gap not in valid_gaps:
            raise ValueError("Gap must be -1, 'del', or 'mask'.")

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
        """Create an alignment path from character indices in the original sequences.

        gap can be a value or "mask", but cannot be "del".
        """
        if gap == "mask":
            return cls.from_bits(np.ma.getmask(indices))
        else:
            return cls.from_bits(indices == gap)

    def to_coordinates(self):
        """Generate an array of segment coordinates in the original sequences."""
        lens = self.lengths * (1 - self.to_bits())
        col0 = np.zeros((self.shape[0], 1), dtype=int)
        lens = np.append(col0, lens, axis=1)
        return lens.cumsum(axis=1)

    @classonlymethod
    def from_coordinates(cls, coords):
        diff = np.diff(coords)
        bits = diff == 0
        lens = diff[bits.argmin(axis=0), np.arange(diff.shape[1])]
        ints = np.packbits(bits, axis=0, bitorder="little")
        if ints.shape[0] == 1:
            ints = np.squeeze(ints, axis=0)
        return cls(lens, ints, np.zeros(diff.shape[0], dtype=int))


class PairAlignPath(AlignPath):
    @overrides(AlignPath)
    @classonlymethod
    def from_bits(cls, bits):
        """Create an alignment path from a bit array."""
        # This should be faster than the generic solution I guess.
        # TODO: Pending benchmark and optimization.

        # Ensure bits is a 2D numpy array of ones and zeros.
        if not isinstance(bits, np.ndarray) or bits.ndim != 2 or bits.shape[1] == 0:
            raise ValueError("Input 'bits' must be a non-empty 2D numpy array.")
        if not np.all(np.logical_or(bits == 0, bits == 1)):
            raise ValueError("Input 'bits' must contain only zeros and ones.")

        ints = bits[0] + (bits[1] << 1)
        idx = np.append(0, np.where(ints[:-1] != ints[1:])[0] + 1)
        lens = np.append(idx[1:] - idx[:-1], ints.size - idx[-1])
        return cls(lens, ints[idx], np.zeros(bits.shape[0], dtype=int))

    @overrides(AlignPath)
    def to_bits(self):
        """Unpack states into an array of bits."""
        # This should be faster than the generic solution I guess.
        # TODO: Pending benchmark and optimization.
        if not np.all(np.isin(self.states, [0, 1, 2])):
            raise ValueError(
                "For pairwise alignment, `states` must only contain "
                "zeros, ones, or twos."
            )
        return np.stack([self.states & 1, self.states >> 1])

    def to_cigar(self, seqs=None):
        """Get a CIGAR string representing the pairwise alignment path.

        seqs: If provided, will distinguish match (=) and mismatch (X). Otherwise,
            will uniformly note them as (mis)match (M).
        """
        cigar = ""
        lengths = self.lengths
        gaps = np.squeeze(self.states)
        codes = ["M", "I", "D"]
        if seqs is not None:
            query = seqs[0]
            ref = seqs[1]
            for qchar, rchar in zip(query, ref):
                if qchar == rchar:
                    cigar += "="
                elif qchar == "-":
                    cigar += "I"
                elif rchar == "-":
                    cigar += "D"
                else:
                    cigar += "X"
            return self._run_length_encode(cigar)
        else:
            for i, length in enumerate(lengths):
                cigar += str(length) + codes[gaps[i]]
            return cigar

    @classonlymethod
    def from_cigar(cls, cigar):
        """Create a pairwise alignment path from a CIGAR string."""
        # need to have ability for strings without 1's
        # Make sure cigar is not empty.
        if not cigar:
            raise ValueError("CIGAR string must not be empty.")

        # Determine whether or not to fix the arrays.
        fix_arrs = False
        if "=" in cigar or "X" in cigar:
            fix_arrs = True

        cigar = cigar.replace("=", "M").replace("X", "M")
        lengths = []
        gaps = []
        current_length = 0
        mapping = {"M": 0, "I": 1, "D": 2}
        no_ones = True
        for char in cigar:
            if char.isdigit():
                no_ones = False
                current_length = current_length * 10 + int(char)
            elif char in mapping:
                if no_ones:
                    lengths.append(current_length + 1)
                else:
                    lengths.append(current_length)
                gaps.append(mapping[char])
                current_length = 0
                no_ones = True
            else:
                raise ValueError("Invalid characters in CIGAR string.")
        if fix_arrs:
            lengths, gaps = cls._fix_arrays(
                lengths=np.array(lengths), gaps=np.array(gaps)
            )
        return cls(np.asarray(lengths), np.asarray(gaps), [0, 0])

    def _run_length_encode(self, input_string):
        """Perform run length encoding on a string."""
        input_array = np.array(list(input_string))
        change_indices = np.append(
            0, np.where(input_array[:-1] != input_array[1:])[0] + 1
        )
        count = np.diff(np.concatenate((change_indices, [len(input_string)])))
        unique_chars = input_array[change_indices]
        encoded_string = "".join(str(c) + u for c, u in zip(count, unique_chars))

        return encoded_string

    def _fix_arrays(lengths, gaps):
        """Merge consecutive same values from gaps array and sum corresponding values
        in lengths array."""
        idx = np.diff(gaps, prepend=np.array([True])) != 0
        gaps_out = gaps[idx]
        groups = np.cumsum(idx)
        lengths_out = np.bincount(groups, weights=lengths).astype(int)[1:]

        return lengths_out, gaps_out
