# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.util._decorator import overrides, classonlymethod


class AlignPath:
    def __init__(self, lengths, states, n_seqs, starts=None):
        """Create an alignment path from segment lengths and states."""
        # Number of sequences needs to be explicitly provided, because the packed bits
        # does not contain this information. (It is merely in multiples of 8.)
        self.lengths = np.asarray(lengths, dtype=int)
        self.states = np.asarray(states, dtype=np.uint8)

        # Shape is n_seqs (rows) x n_positions (columns), which is consistent with
        # TabularMSA
        self.shape = (n_seqs, self.lengths.sum())

        # TODO: An additional parameter `starts` should record the starting positions
        # of the alignment in each sequence (shape: (n_seqs,)). This is important
        # especially for local alignments (e.g., search for a short read in a genome
        # sequence). If it is provided, `n_seqs` is no longer necessary.

        # TODO: Needs to think about whether reverse complemented (Boolean array of
        # (n_seqs,)) should be included as a parameter. It is only relevant for
        # nucleotide sequences.

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
        return cls(lens, ints, bits.shape[0])

    @classonlymethod
    def from_tabular(cls, msa):
        """Create an alignment path from a `TabularMSA` object."""
        # Convert TabularMSA in to a 2D array of bytes.
        # TODO: TabularMSA itself should be refactored to have this as the default data
        # structure.
        byte_arr = np.stack([x._bytes for x in msa._seqs])

        # Get gap character code.
        # TODO: current code only handles default gap character. Not sure if other gap
        # characters matter. Should add.
        gap_char = ord(msa.dtype.default_gap_char)

        # Identify gap positions, and convert them into a bit array, then create an
        # alignment path based on it.
        return cls.from_bits(byte_arr == gap_char)

    def to_indices(self, gap=-1):
        """Generate an array of indices of characters in the original sequences.

        gap = "del": delete columns that have any gap, "mask": mask gaps, others
        (default: -1):
        fill gaps with this value.
        """
        bits = self.to_bits()
        # TODO: Consider optimization using np.arange.
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
        return cls(lens, ints, diff.shape[0])


class PairAlignPath(AlignPath):
    @overrides(AlignPath)
    @classonlymethod
    def from_bits(cls, bits):
        """Create an alignment path from a bit array."""
        # This should be faster than the generic solution I guess.
        # TODO: Pending benchmark and optimization.
        ints = bits[0] + (bits[1] << 1)
        idx = np.append(0, np.where(ints[:-1] != ints[1:])[0] + 1)
        lens = np.append(idx[1:] - idx[:-1], ints.size - idx[-1])
        return cls(lens, ints[idx], bits.shape[0])

    @overrides(AlignPath)
    def to_bits(self):
        """Unpack states into an array of bits."""
        # This should be faster than the generic solution I guess.
        # TODO: Pending benchmark and optimization.
        return np.stack([self.states & 1, self.states >> 1])

    def to_cigar(self, seqs=None):
        """Get a CIGAR string representing the pairwise alignment path.

        seqs: If provided, will distinguish match (=) and mismatch (X). Otherwise,
            will uniformly note them as (mis)match (M).
        """
        cigar = ""
        lengths = self.lengths
        gaps = self.states
        codes = ["M", "I", "D"]
        # Leaving this out for now, but need to implement "=" and "X" on a per
        # character basis.
        if seqs is not None:
            pass
        else:
            for i, length in enumerate(lengths):
                cigar += str(length) + codes[gaps[i]]
            return cigar

    @classonlymethod
    def from_cigar(self, cigar):
        """Create a pairwise alignment path from a CIGAR string."""
        if "=" in cigar or "X" in cigar:
            self.fix_arrs = True
        cigar = cigar.replace("=", "M").replace("X", "M")
        lengths = []
        gaps = []
        current_length = 0
        mapping = {"M": 0, "I": 1, "D": 2}
        for char in cigar:
            if char.isdigit():
                current_length = current_length * 10 + int(char)
            elif char in mapping:
                lengths.append(current_length)
                gaps.append(mapping[char])
                current_length = 0
        if self.fix_arrs:
            lengths, gaps = self._fix_arrays(np.array(lengths), np.array(gaps))
        return lengths, gaps

    def _fix_arrays(lengths, gaps):
        """Fix output arrays if subsequent '=', 'X', or 'M' are present in the input
        CIGAR string."""
        is_zero = np.concatenate(([0], np.equal(gaps, 0).view(np.int8), [0]))
        abs_diff = np.abs(np.diff(is_zero))

        ranges = np.where(abs_diff == 1)[0].reshape(-1, 2)

        for start, stop in ranges:
            if stop - start > 1:
                lengths[start] += np.sum(lengths[start + 1 : stop])

        identical_indices = np.where(gaps[:-1] == gaps[1:])[0] + 1
        print(identical_indices)
        gaps = np.delete(gaps, identical_indices)
        lengths = np.delete(lengths, identical_indices)

        return lengths, gaps
