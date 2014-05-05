#!/usr/bin/env python
"""Writer for FASTA sequence format"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class _fake_seq(str):

    """a holder for string sequences that allows provision of a seq.Label
    attribute, required by fasta formatting funcs."""

    def __new__(cls, Label, Seq):
        new = str.__new__(cls, Seq)
        new.Label = Label
        return new

    def __getslice__(self, *args, **kwargs):
        new_seq = str.__getslice__(self, *args, **kwargs)
        return self.__new__(self.__class__, self.Label, new_seq)


def fasta_from_sequences(seqs, make_seqlabel=None, line_wrap=None):
    """Returns a FASTA string given a list of sequence objects.

    A ``sequence.Label`` attribute takes precedence over ``sequence.Name``.

    Parameters
    ----------
    seqs : list
        seqs can be a list of sequence objects or strings.
    make_seqlabel : function, optional
        callback function that takes the seq object and returns a label ``str``
        . If ``None`` is passed, the following attributes will try to be
        retrieved in this order and the first to exist will be used:
        ``identifier``, ``Label`` or ``Name``. In any other case an integer
        with the position of the sequence object will be used.
    line_wrap : int, optional
        line_wrap: a integer for maximum line width, if ``None`` is passed the
        full sequence will be used.

    Returns
    -------
    str
        FASTA formatted string composed of the objects passed in via `seqs`.

    Examples
    --------
    Formatting a list of sequence objects

    >>> from skbio.format.sequences import fasta_from_sequences
    >>> from skbio.core.sequence import DNASequence
    >>> seqs = [DNASequence('ACTCGAGATC', 'seq1'),
    ...         DNASequence('GGCCT', 'seq2')]
    >>> fasta_from_sequences(seqs)
    >seq1
    ACTCGAGATC
    >seq2
    GGCCT

    See Also
    --------
    skbio.parse.sequences.parse_fasta

    """
    fasta_list = []
    for i, seq in enumerate(seqs):
        # Check if it has a label, or one is to be created
        label = str(i)
        if make_seqlabel is not None:
            label = make_seqlabel(seq)
        elif hasattr(seq, 'identifier') and seq.identifier:
            label = seq.identifier
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
            body = ["%s" % seq_str[j * line_wrap:(j + 1) * line_wrap]
                    for j in range(numlines)]
        else:
            body = ["%s" % seq_str]

        fasta_list.append('>' + label)
        fasta_list += body

    return '\n'.join(fasta_list)


def fasta_from_alignment(aln, make_seqlabel=None, line_wrap=None, sort=True):
    """Returns a FASTA string given an alignment.

        - aln can be an Alignment object or dict.
        - make_seqlabel: callback function that takes the seq object and
          returns a label str
        - line_wrap: a integer for maximum line width
    """
    # get seq output order
    try:
        order = aln.Names[:]
    except AttributeError:
        order = aln.keys()

    if sort:
        order = sorted(order)

    try:
        seq_dict = aln.NamedSeqs
    except AttributeError:
        seq_dict = aln

    ordered_seqs = []
    for label in order:
        seq = seq_dict[label]
        if isinstance(seq, str):
            seq = _fake_seq(label, seq)
        ordered_seqs.append(seq)
    return fasta_from_sequences(ordered_seqs, make_seqlabel=make_seqlabel,
                                line_wrap=line_wrap)
