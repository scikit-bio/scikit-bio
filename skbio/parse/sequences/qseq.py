#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function

import collections

from skbio.util.io import open_file

from .fastq import ascii_to_phred33, ascii_to_phred64


def parse_qseq(infile, phred_offset=33):
    r"""Generator of seq ids, seqs, quals and other records from a qseq file.

    Parameters
    ----------
    infile : open file object or str
        An open qseq file or a path to a qseq file.
    phred_offset : {33, 64}, optional
        What Phred offset to use when converting qual score symbols
        to integers.

    Returns
    -------
    four-item tuple: (str, str, np.array(dtype=int), namedtuple)
        yields the sequence id, sequence, qual array and other record
        information for each entry.  The sequence ID format is:
        <Machine name>_<Run number>:<Lane number>:<Tile number>:<x>:<y>#
        <Index>/<Read number>.  The namedtuple attributes are:
        machine_name, run, lane, tile, x, y, index, read and filtered.

    Examples
    --------
    Assume we have a qseq-formatted file with the following contents::

        CRESSIA       242     1       2204    1453    1918    0       1
            .TTAATAAGAATGTCTGTTGTGGCTTAAAA  B[[[W][Y[Zccccccccc\cccac_____  1
        CRESSIA       242     1       2204    1490    1921    0       2
            ..GTAAAACCCATATATTGAAAACTACAAA  BWUTWcXVXXcccc_cccccccccc_cccc  1

    >>> from future.utils.six import StringIO
    >>> qseq_f = StringIO('CRESSIA\t242\t1\t2204\t1453\t1918\t0\t1\t'
    ...   '.TTAATAAGAATGTCTGTTGTGGCTTAAAA\tB[[[W][Y[Zccccccccc\cccac_____\t1\n'
    ...                   'CRESSIA\t242\t1\t2204\t1490\t1921\t0\t2\t'
    ...   '..GTAAAACCCATATATTGAAAACTACAAA\tBWUTWcXVXXcccc_cccccccccc_cccc\t1\n'
    ... )

    We can parse this as follows:

    >>> from skbio.parse.sequences import parse_qseq
    >>> for seq_id, seq, qual, record in parse_qseq(qseq_f, phred_offset=64):
    ...     print(seq_id)
    ...     print(seq)
    ...     print(qual[:10])
    ...     print(record.run)
    ...     print(record.lane)
    CRESSIA_242:1:2204:1453:1918#0/1
    .TTAATAAGAATGTCTGTTGTGGCTTAAAA
    [ 2 27 27 27 23 29 27 25 27 26]
    242
    1
    CRESSIA_242:1:2204:1490:1921#0/2
    ..GTAAAACCCATATATTGAAAACTACAAA
    [ 2 23 21 20 23 35 24 22 24 24]
    242
    1
    """
    if phred_offset == 33:
        phred_f = ascii_to_phred33
    elif phred_offset == 64:
        phred_f = ascii_to_phred64
    else:
        raise ValueError("Unknown PHRED offset of %s" % phred_offset)

    # namedtuple to store all other record information
    Record = collections.namedtuple(
        'Record',
        ['machine_name',
         'run',
         'lane',
         'tile',
         'x',
         'y',
         'index',
         'read',
         'filtered'])

    with open_file(infile) as lines:
        for rec in lines:
            try:
                rec = str(rec.decode('utf-8'))
            except AttributeError:
                pass
            (machine_name, run, lane, tile, x, y, index, read, seq, qual,
             filtered) = rec.split()
            # sequence ID is formatted using the first eight items.
            seq_id = '%s_%s:%s:%s:%s:%s#%s/%s' % (
                machine_name, run, lane, tile, x, y, index, read)
            # qual string is converted to an array of ints.
            qual = phred_f(qual)
            # other items are returned as a namedtuple
            record = Record(
                machine_name=machine_name,
                run=int(run),
                lane=int(lane),
                tile=int(tile),
                x=int(x),
                y=int(y),
                index=int(index),
                read=int(read),
                filtered=bool(int(filtered)))

            yield seq_id, seq, qual, record
