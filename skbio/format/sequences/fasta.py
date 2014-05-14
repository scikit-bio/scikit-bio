#!/usr/bin/env python
"""Writer for FASTA sequence format"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.core.alignment import Alignment
from skbio.core.sequence import BiologicalSequence


def fasta_from_sequences(seqs, make_seqlabel=None, line_wrap=None):
    """Returns a FASTA string given a list of sequence objects.

    A ``sequence.Label`` attribute takes precedence over ``sequence.Name``.

    Parameters
    ----------
    seqs : list
        seqs can be a list of sequence objects or strings.
    make_seqlabel : function, optional
        callback function that takes the seq object and returns a label
        ``str``. If ``None`` is passed, the following attributes will try to be
        retrieved in this order and the first to exist will be used:
        ``id``, ``Label`` or ``Name``. In any other case an integer
        with the position of the sequence object will be used.
    line_wrap : int, optional
        line_wrap: a integer for maximum line width, if ``None`` is passed the
        full sequence will be used.

    Returns
    -------
    str
        FASTA formatted string composed of the objects passed in via `seqs`.

    See Also
    --------
    skbio.parse.sequences.parse_fasta

    Examples
    --------
    Formatting a list of sequence objects

    >>> from skbio.format.sequences import fasta_from_sequences
    >>> from skbio.core.sequence import DNASequence
    >>> seqs = [DNASequence('ACTCGAGATC', 'seq1'),
    ...         DNASequence('GGCCT', 'seq2')]
    >>> print fasta_from_sequences(seqs)
    >seq1
    ACTCGAGATC
    >seq2
    GGCCT

    """
    fasta_list = []
    for i, seq in enumerate(seqs):
        # Check if it has a label, or one is to be created
        label = str(i)
        if make_seqlabel is not None:
            label = make_seqlabel(seq)
        elif hasattr(seq, 'id') and seq.id:
            label = seq.id
        elif hasattr(seq, 'Label') and seq.Label:
            label = seq.Label
        elif hasattr(seq, 'Name') and seq.Name:
            label = seq.Name

        # wrap sequence lines
        seq_str = str(seq)
        if line_wrap is not None:
            numlines, remainder = divmod(len(seq_str), line_wrap)
            if remainder:
                numlines += 1
            body = [seq_str[j * line_wrap:(j + 1) * line_wrap]
                    for j in range(numlines)]
        else:
            body = [seq_str]

        fasta_list.append('>' + label)
        fasta_list += body

    return '\n'.join(fasta_list)


def fasta_from_alignment(aln, make_seqlabel=None, line_wrap=None, sort=True):
    """Returns a FASTA string given an alignment object

    Parameters
    ----------
    aln : Alignment, dict
        alignment or dictionary where the keys are the sequence ids and
        the values are the sequences themselves.
    make_seqlabel : function, optional
        callback function that takes the seq object and returns a label
        ``str``. If ``None`` is passed, the following attributes will try to be
        retrieved in this order and the first to exist will be used:
        ``id``, ``Label`` or ``Name``. In any other case an integer
        with the position of the sequence object will be used.
    line_wrap : int, optional
        line_wrap: a integer for maximum line width, if ``None`` is passed the
        full sequence will be used.
    sort : bool, optional
        Whether or not the sequences should be sorted by their sequence
        id, default value is ``True``.

    Returns
    -------
    str
        FASTA formatted string composed of the objects passed in via `seqs`.

    See Also
    --------
    skbio.parse.sequences.parse_fasta
    skbio.core.alignment.Alignment

    Examples
    --------
    Formatting a sequence alignment object into a FASTA file.
    >>> from skbio.core.alignment import Alignment
    >>> from skbio.core.sequence import DNA
    >>> from skbio.format.sequences import fasta_from_alignment
    >>> seqs = [DNA("ACC--G-GGTA..", id="seq1"),
    ...     DNA("TCC--G-GGCA..", id="seqs2")]
    >>> a1 = Alignment(seqs)
    >>> print fasta_from_alignment(a1)
    >seq1
    ACC--G-GGTA..
    >seqs2
    TCC--G-GGCA..

    """
    # check if it's an Alignment object or a dictionary
    if isinstance(aln, Alignment):
        order = aln.ids()
    else:
        order = aln.keys()

    if sort:
        order = sorted(order)

    ordered_seqs = []
    for label in order:
        seq = aln[label]
        if isinstance(seq, str):
            seq = BiologicalSequence(seq, label)
        ordered_seqs.append(seq)
    return fasta_from_sequences(ordered_seqs, make_seqlabel=make_seqlabel,
                                line_wrap=line_wrap)
