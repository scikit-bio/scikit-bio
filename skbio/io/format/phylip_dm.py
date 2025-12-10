"""PHYLIP distance matrix format (:mod:`skbio.io.format.phylip_dm`)
===============================================================

.. currentmodule:: skbio.io.format.phylip_dm

The PHYLIP file format can store a distance matrix.

An example PHYLIP-formatted distance matrix taken from

        7
Bovine    0.0000 1.6866 1.7198 1.6606 1.5243 1.6043 1.5905
Mouse     1.6866 0.0000 1.5232 1.4841 1.4465 1.4389 1.4629
Gibbon    1.7198 1.5232 0.0000 0.7115 0.5958 0.6179 0.5583
Orang     1.6606 1.4841 0.7115 0.0000 0.4631 0.5061 0.4710
Gorilla   1.5243 1.4465 0.5958 0.4631 0.0000 0.3484 0.3083
Chimp     1.6043 1.4389 0.6179 0.5061 0.3484 0.0000 0.2692
Human     1.5905 1.4629 0.5583 0.4710 0.3083 0.2692 0.0000"""

import numpy as np

from skbio.stats.distance import DistanceMatrix
from skbio.io import create_format, PhylipFormatError

phylip_dm = create_format("phylip_dm")


# @phylip_dm.sniffer()
def _phylip_dm_sniffer(fh):
    # Strategy:
    #  Read the header and a single sequence; verify that the sequence length matches
    #  the header information. Don't verify that  the total number of lines matches the
    #  header information, since that would require reading the whole file.
    try:
        header = next(_line_generator(fh))
        n_seqs = _validate_header(header)
        # We check the first three lines after the header to enable sniffing of
        # lower triangle AND square matrices.
        for line_no in range(3):
            line = next(_line_generator(fh))
            _validate_line(line, n_seqs, line_no)
    except (StopIteration, PhylipFormatError):
        return False, {}
    return True, {}


# phylip_dm.reader(DistanceMatrix)
def _phylip_dm_to_distance_matrix(fh, cls=None):
    if cls is None:
        cls = DistanceMatrix
    data = _parse_phylip_dm_raw(fh)
    dists = [x[0] for x in data]
    ids = [x[1] for x in data]
    # If it's in lower triangular form we convert it into condensed form to pass to
    # DistanceMatrix.
    if len(dists[0]) == 0:
        dists = [
            dists[row][col]
            for col in range(len(dists))
            for row in range(col + 1, len(dists))
        ]
    return cls(np.array(dists, dtype=float), ids)


def _line_generator(fh):
    """Remove linebreak characters and yield lines."""
    for line in fh:
        yield line.rstrip("\n")


def _validate_header(header):
    header_vals = header.split()
    try:
        (n_seqs,) = [int(x) for x in header_vals]
        if n_seqs < 1:
            raise PhylipFormatError("The number of sequences must be positive.")
    except ValueError:
        raise PhylipFormatError(
            "Found non-header line when attempting to read the 1st record "
            f"(header line should have a single integer): {header}"
        )
    return n_seqs


def _validate_line(line, n_seqs, n_dists):
    if not line:
        raise PhylipFormatError("Empty lines are not allowed.")
    split_line = line.split()
    id = split_line[0]
    dists = split_line[1:]
    # This check handles lower triangle matrices. We expects 0 distances on the first
    # non-header line, a single distance on the second non-header line, and two
    # distances on the third non-header line.
    if len(dists) != n_dists:
        # If there are more distances than expected for a lower triangle matrix, we
        # expect that it is a square matrix. In this case we check that the number of
        # distances matches the value specified in the header.
        if len(dists) != n_seqs:
            raise PhylipFormatError(
                f"The number of distances {dists} is not {n_seqs} as specified in the "
                "header."
            )

    return (dists, id)


def _parse_phylip_dm_raw(fh):
    """Raw parser for PHYLIP formatted distance matrix files.

    Returns a list of raw (dists, id) values."""
    # File should have a single header on the first line.
    try:
        header = next(_line_generator(fh))
    except StopIteration:
        raise PhylipFormatError("This file is empty.")
    n_seqs = _validate_header(header)

    # All following lines should be ID+distances. No blank lines are allowed.
    n_dists = 0
    data = []
    for line in _line_generator(fh):
        data.append(_validate_line(line, n_seqs, n_dists))
        n_dists += 1
    if len(data) != n_seqs:
        raise PhylipFormatError(
            f"The number of sequences is not {n_seqs} as specified in the header."
        )
    return data
