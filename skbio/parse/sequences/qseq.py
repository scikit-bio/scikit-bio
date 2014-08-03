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

from skbio.core.exception import RecordError
from skbio.util.io import open_file

from .fastq import ascii_to_phred64


def parse_qseq(infile, strict=True):
    r"""Generator of sequence ids, sequences, quals and other records
    from a qseq file.

    Parameters
    ----------
    infile : open file object or str
        An open qseq file or a path to a qseq file.

    strict : bool
        If ``True`` a ``RecordError`` will be raised if there is a record
        without enough items.  If ``False``, partial records will be skipped.

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

    >>> from StringIO import StringIO
    >>> qseq_f = StringIO('CRESSIA\t242\t1\t2204\t1453\t1918\t0\t1\t'
    ...   '.TTAATAAGAATGTCTGTTGTGGCTTAAAA\tB[[[W][Y[Zccccccccc\cccac_____\t1\n'
    ...                   'CRESSIA\t242\t1\t2204\t1490\t1921\t0\t2\t'
    ...   '..GTAAAACCCATATATTGAAAACTACAAA\tBWUTWcXVXXcccc_cccccccccc_cccc\t1\n'
    ... )

    We can parse this as follows:

    >>> from skbio.parse.sequences import parse_qseq
    >>> for seq_id, seq, qual, record in parse_qseq(qseq_f):
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
            rec = rec.strip('\n').split('\t')
            # record must have at least ten items
            if len(rec) < 11:
                if strict:
                    raise RecordError(
                        "Found line without enough items: %s" % rec)
                else:
                    continue
            # sequence ID is formatted using the first eight items.
            seq_id = '%s_%s:%s:%s:%s:%s#%s/%s' % (
                rec[0], rec[1], rec[2], rec[3], rec[4], rec[5], rec[6], rec[7])
            # sequence is the eight item.
            seq = rec[8]
            # qual string is converted to an array of ints.
            qual = ascii_to_phred64(rec[9])
            # other items are returned as a namedtuple
            record = Record(
                machine_name=rec[0],
                run=int(rec[1]),
                lane=int(rec[2]),
                tile=int(rec[3]),
                x=int(rec[4]),
                y=int(rec[5]),
                index=int(rec[6]),
                read=int(rec[7]),
                filtered=bool(int(rec[10])))

            yield seq_id, seq, qual, record
