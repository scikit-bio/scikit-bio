# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import functools
import contextlib

import pandas as pd

_possible_columns = {
    "qseqid": str,
    "qgi": float,
    "qacc": str,
    "qaccver": str,
    "qlen": float,
    "sseqid": str,
    "sallseqid": str,
    "sgi": float,
    "sallgi": float,
    "sacc": str,
    "saccver": str,
    "sallacc": str,
    "slen": float,
    "qstart": float,
    "qend": float,
    "sstart": float,
    "send": float,
    "qseq": str,
    "sseq": str,
    "evalue": float,
    "bitscore": float,
    "score": float,
    "length": float,
    "pident": float,
    "nident": float,
    "mismatch": float,
    "positive": float,
    "gapopen": float,
    "gaps": float,
    "ppos": float,
    "frames": str,
    "qframe": float,
    "sframe": float,
    "btop": float,
    "staxids": str,
    "sscinames": str,
    "scomnames": str,
    "sblastnames": str,
    "sskingdoms": str,
    "stitle": str,
    "salltitles": str,
    "sstrand": str,
    "qcovs": float,
    "qcovhsp": float,
}


def _parse_blast_data(fh, columns, error, error_message, comment=None, skiprows=None):
    read_csv = functools.partial(
        pd.read_csv,
        na_values="N/A",
        sep="\t",
        header=None,
        keep_default_na=False,
        comment=comment,
        skiprows=skiprows,
    )

    # HACK for https://github.com/pandas-dev/pandas/issues/14418
    # this avoids closing the `fh`, whose lifetime isn't the responsibility
    # of this parser
    with _noop_close(fh) as fh:
        lineone = read_csv(fh, nrows=1)

        if len(lineone.columns) != len(columns):
            raise error(error_message % (len(columns), len(lineone.columns)))

        fh.seek(0)

        return read_csv(fh, names=columns, dtype=_possible_columns)


# HACK for https://github.com/pandas-dev/pandas/issues/14418
@contextlib.contextmanager
def _noop_close(fh):
    backup = fh.close
    fh.close = lambda: None
    try:
        yield fh
    finally:
        fh.close = backup
