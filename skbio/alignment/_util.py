def align_score(path, seqs, match=2, mismatch=-3, gap_open=-5, gap_extend=-2):
    """Calculate alignment score of a pairwise alignment."""
    # The current code only accepts an alignment path and the original sequences. It is
    # optimized for this data structure.
    # TODO: It should also accept a tabular MSA object. For optimal performance, it
    # should be a separate function. A dispatcher will unite the two functions.
    if path.shape[0] != 2:
        raise ValueError("Only pairwise alignment is accepted.")

    # step 1: calculate gap penalties
    # based on the number and lengths of gap segments
    # agnostic of the sequences
    gap_segs = path.states != 0
    n_gap_opens = gap_segs.sum()
    n_gap_extends = (gap_segs * path.lengths).sum()

    # step 2: calculate match & mismatch scores
    # based on the characters within the non-gap segments
    idx = path.to_indices(gap="del")
    seq1, seq2 = (x._bytes[i] for x, i in zip(seqs, idx))
    n_matches = (seq1 == seq2).sum()
    n_mismatches = seq1.size - n_matches

    # get total score
    return (
        match * n_matches
        + mismatch * n_mismatches
        + gap_open * n_gap_opens
        + gap_extend * n_gap_extends
    )


def align_scores(path, seqs, match, mismatch, gap_open, gap_extend):
    """Calculate pairwise alignment scores of a multiple alignment and return a
    matrix."""
    pass


def align_sp_scores(path, seqs, match, mismatch, gap_open, gap_extend):
    """Calculate sum-of-pairs (SP) score of a multiple alignment."""
    pass


def separate_rle(msa):
    byte_arr = np.stack([x._bytes for x in msa._seqs])
    gap_char = ord(msa.dtype.default_gap_char)
    diff = np.diff(byte_arr == gap_char, prepend=np.nan, axis=1)
    return [np.diff(np.where(x)[0], append=byte_arr.shape[1]) for x in diff]
