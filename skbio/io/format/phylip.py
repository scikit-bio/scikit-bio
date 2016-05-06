"""
PHYLIP multiple sequence alignment format (:mod:`skbio.io.format.phylip`)
=========================================================================

.. currentmodule:: skbio.io.format.phylip

The PHYLIP file format stores a multiple sequence alignment. The format was
originally defined and used in Joe Felsenstein's PHYLIP package [1]_, and has
since been supported by several other bioinformatics tools (e.g., RAxML [2]_).
See [3]_ for the original format description, and [4]_ and [5]_ for additional
descriptions.

An example PHYLIP-formatted file taken from [3]_::

          5    42
    Turkey    AAGCTNGGGC ATTTCAGGGT GAGCCCGGGC AATACAGGGT AT
    Salmo gairAAGCCTTGGC AGTGCAGGGT GAGCCGTGGC CGGGCACGGT AT
    H. SapiensACCGGTTGGC CGTTCAGGGT ACAGGTTGGC CGTTCAGGGT AA
    Chimp     AAACCCTTGC CGTTACGCTT AAACCGAGGC CGGGACACTC AT
    Gorilla   AAACCCTTGC CGGTACGCTT AAACCATTGC CGGTACGCTT AA

.. note:: Original copyright notice for the above PHYLIP file:

   *(c) Copyright 1986-2008 by The University of Washington. Written by Joseph
   Felsenstein. Permission is granted to copy this document provided that no
   fee is charged for it and that this copyright notice is not removed.*

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.alignment.TabularMSA`                              |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
PHYLIP format is a plain text format containing exactly two sections: a header
describing the dimensions of the alignment, followed by the multiple sequence
alignment itself.

The format described here is "strict" PHYLIP, as described in [4]_. Strict
PHYLIP requires that each sequence identifier is exactly 10 characters long
(padded with spaces as necessary). Other bioinformatics tools (e.g., RAxML) may
relax this rule to allow for longer sequence identifiers. See the
**Alignment Section** below for more details.

The format described here is "sequential" format. The original PHYLIP format
specification [3]_ describes both sequential and interleaved formats.

.. note:: scikit-bio currently supports reading and writing strict, sequential
   PHYLIP-formatted files. Relaxed and/or interleaved PHYLIP formats are not
   supported.

Header Section
^^^^^^^^^^^^^^
The header consists of a single line describing the dimensions of the
alignment. It **must** be the first line in the file. The header consists of
optional spaces, followed by two positive integers (``n`` and ``m``) separated
by one or more spaces. The first integer (``n``) specifies the number of
sequences (i.e., the number of rows) in the alignment. The second integer
(``m``) specifies the length of the sequences (i.e., the number of columns) in
the alignment. The smallest supported alignment dimensions are 1x1.

.. note:: scikit-bio will write the PHYLIP format header *without* preceding
   spaces, and with only a single space between ``n`` and ``m``.

   PHYLIP format *does not* support blank line(s) between the header and the
   alignment.

Alignment Section
^^^^^^^^^^^^^^^^^
The alignment section immediately follows the header. It consists of ``n``
lines (rows), one for each sequence in the alignment. Each row consists of a
sequence identifier (ID) and characters in the sequence, in fixed width format.

The sequence ID can be up to 10 characters long. IDs less than 10 characters
must have spaces appended to them to reach the 10 character fixed width. Within
an ID, all characters except newlines are supported, including spaces,
underscores, and numbers.

.. note:: When reading a PHYLIP-formatted file into an
   ``skbio.alignment.TabularMSA`` object, sequence identifiers/labels are
   stored as ``TabularMSA`` index labels (``index`` property).

   When writing an ``skbio.alignment.TabularMSA`` object as a PHYLIP-formatted
   file, ``TabularMSA`` index labels will be converted to strings and written
   as sequence identifiers/labels.

   scikit-bio supports the empty string (``''``) as a valid sequence ID. An
   empty ID will be padded with 10 spaces when writing.

Sequence characters immediately follow the sequence ID. They *must* start at
the 11th character in the line, as the first 10 characters are reserved for the
sequence ID. While PHYLIP format does not explicitly restrict the set of
supported characters that may be used to represent a sequence, the original
format description [3]_ specifies the IUPAC nucleic acid lexicon for DNA or RNA
sequences, and the IUPAC protein lexicon for protein sequences. The original
PHYLIP specification uses ``-`` as a gap character, though older versions also
supported ``.``. The sequence characters may contain optional spaces (e.g., to
improve readability), and both upper and lower case characters are supported.

.. note:: scikit-bio will read/write a PHYLIP-formatted file as long as the
   alignment's sequence characters are valid for the type of in-memory sequence
   object being read into or written from. This differs from the PHYLIP
   specification, which states that a PHYLIP-formatted file can only contain
   valid IUPAC characters. See the ``constructor`` format parameter below for
   details.

   Since scikit-bio supports both ``-`` and ``.`` as gap characters (e.g., in
   ``DNA``, ``RNA``, and ``Protein`` sequence objects), both are supported when
   reading/writing a PHYLIP-formatted file.

   When writing a PHYLIP-formatted file, scikit-bio will split up each sequence
   into chunks that are 10 characters long. Each chunk will be separated by a
   single space. The sequence will always appear on a single line (sequential
   format). It will *not* be wrapped across multiple lines. Sequences are
   chunked in this manner for improved readability, and because most example
   PHYLIP files are chunked in a similar way (e.g., see the example file
   above). Note that this chunking is not required when reading
   PHYLIP-formatted files, nor by the PHYLIP format specification itself.

Format Parameters
-----------------
The only supported format parameter is ``constructor``, which specifies the
type of in-memory sequence object to read each aligned sequence into. This must
be a subclass of ``GrammaredSequence`` (e.g., ``DNA``, ``RNA``, ``Protein``)
and is a required format parameter. For example, if you know that the PHYLIP
file you're reading contains DNA sequences, you would pass ``constructor=DNA``
to the reader call.

Examples
--------
Let's create a ``TabularMSA`` with three DNA sequences:

>>> from skbio import TabularMSA, DNA
>>> seqs = [DNA('ACCGTTGTA-GTAGCT', metadata={'id':'seq1'}),
...         DNA('A--GTCGAA-GTACCT', metadata={'id':'sequence-2'}),
...         DNA('AGAGTTGAAGGTATCT', metadata={'id':'3'})]
>>> msa = TabularMSA(seqs, minter='id')
>>> msa
TabularMSA[DNA]
----------------------
Stats:
    sequence count: 3
    position count: 16
----------------------
ACCGTTGTA-GTAGCT
A--GTCGAA-GTACCT
AGAGTTGAAGGTATCT
>>> msa.index
Index(['seq1', 'sequence-2', '3'], dtype='object')

Now let's write the ``TabularMSA`` to file in PHYLIP format and take a look at
the output:

>>> from io import StringIO
>>> fh = StringIO()
>>> print(msa.write(fh, format='phylip').getvalue())
3 16
seq1      ACCGTTGTA- GTAGCT
sequence-2A--GTCGAA- GTACCT
3         AGAGTTGAAG GTATCT
<BLANKLINE>
>>> fh.close()

Notice that the 16-character sequences were split into two chunks, and that
each sequence appears on a single line (sequential format). Also note that each
sequence ID is padded with spaces to 10 characters in order to produce a fixed
width column.

If the index labels in a ``TabularMSA`` surpass the 10-character limit, an
error will be raised when writing:

>>> msa.index = ['seq1', 'long-sequence-2', 'seq3']
>>> fh = StringIO()
>>> msa.write(fh, format='phylip')
Traceback (most recent call last):
    ...
skbio.io._exception.PhylipFormatError: ``TabularMSA`` can only be written in \
PHYLIP format if all sequence index labels have 10 or fewer characters. Found \
sequence with index label 'long-sequence-2' that exceeds this limit. Use \
``TabularMSA.reassign_index`` to assign shorter index labels.
>>> fh.close()

One way to work around this is to assign shorter index labels. The recommended
way to do this is via ``TabularMSA.reassign_index``. For example, to reassign
default integer index labels:

>>> msa.reassign_index()
>>> msa.index
RangeIndex(start=0, stop=3, step=1)

We can now write the ``TabularMSA`` in PHYLIP format:

>>> fh = StringIO()
>>> print(msa.write(fh, format='phylip').getvalue())
3 16
0         ACCGTTGTA- GTAGCT
1         A--GTCGAA- GTACCT
2         AGAGTTGAAG GTATCT
<BLANKLINE>
>>> fh.close()

References
----------
.. [1] http://evolution.genetics.washington.edu/phylip.html
.. [2] RAxML Version 8: A tool for Phylogenetic Analysis and
   Post-Analysis of Large Phylogenies". In Bioinformatics, 2014
.. [3] http://evolution.genetics.washington.edu/phylip/doc/sequence.html
.. [4] http://www.phylo.org/tools/obsolete/phylip.html
.. [5] http://www.bioperl.org/wiki/PHYLIP_multiple_alignment_format

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.alignment import TabularMSA
from skbio.io import create_format, PhylipFormatError
from skbio.util._misc import chunk_str


phylip = create_format('phylip')


@phylip.sniffer()
def _phylip_sniffer(fh):
    # Strategy:
    #   Read the header and a single sequence; verify that the sequence length
    #   matches the header information.  Do not verify that the total number of
    #   lines matches the header information, since that would require reading
    #   the whole file.
    try:
        header = next(_line_generator(fh))
        _, seq_len = _validate_header(header)
        line = next(_line_generator(fh))
        _validate_line(line, seq_len)
    except (StopIteration, PhylipFormatError):
        return False, {}
    return True, {}


@phylip.reader(TabularMSA)
def _phylip_to_tabular_msa(fh, constructor=None):
    if constructor is None:
        raise ValueError("Must provide `constructor`.")

    seqs = []
    index = []
    for seq, ID in _parse_phylip_raw(fh):
        seqs.append(constructor(seq))
        index.append(ID)
    return TabularMSA(seqs, index=index)


@phylip.writer(TabularMSA)
def _tabular_msa_to_phylip(obj, fh):
    sequence_count = obj.shape.sequence
    if sequence_count < 1:
        raise PhylipFormatError(
            "TabularMSA can only be written in PHYLIP format if there is at "
            "least one sequence in the alignment.")

    sequence_length = obj.shape.position
    if sequence_length < 1:
        raise PhylipFormatError(
            "TabularMSA can only be written in PHYLIP format if there is at "
            "least one position in the alignment.")

    chunk_size = 10
    labels = [str(label) for label in obj.index]
    for label in labels:
        if len(label) > chunk_size:
            raise PhylipFormatError(
                "``TabularMSA`` can only be written in PHYLIP format if all "
                "sequence index labels have %d or fewer characters. Found "
                "sequence with index label '%s' that exceeds this limit. Use "
                "``TabularMSA.reassign_index`` to assign shorter index labels."
                % (chunk_size, label))

    fh.write('{0:d} {1:d}\n'.format(sequence_count, sequence_length))

    fmt = '{0:%d}{1}\n' % chunk_size
    for label, seq in zip(labels, obj):
        chunked_seq = chunk_str(str(seq), chunk_size, ' ')
        fh.write(fmt.format(label, chunked_seq))


def _validate_header(header):
    header_vals = header.split()
    try:
        n_seqs, seq_len = [int(x) for x in header_vals]
        if n_seqs < 1 or seq_len < 1:
            raise PhylipFormatError(
                'The number of sequences and the length must be positive.')
    except ValueError:
        raise PhylipFormatError(
            'Found non-header line when attempting to read the 1st record '
            '(header line should have two space-separated integers): '
            '"%s"' % header)
    return n_seqs, seq_len


def _validate_line(line, seq_len):
    if not line:
        raise PhylipFormatError("Empty lines are not allowed.")
    ID = line[:10].strip()
    seq = line[10:].replace(' ', '')
    if len(seq) != seq_len:
        raise PhylipFormatError(
            "The length of sequence %s is not %s as specified in the header."
            % (ID, seq_len))
    return (seq, ID)


def _parse_phylip_raw(fh):
    """Raw parser for PHYLIP files.

    Returns a list of raw (seq, id) values.  It is the responsibility of the
    caller to construct the correct in-memory object to hold the data.

    """
    # Note: this returns the full data instead of yielding each sequence,
    # because the header specifies the number of sequences, so the file cannot
    # be validated until it's read completely.

    # File should have a single header on the first line.
    try:
        header = next(_line_generator(fh))
    except StopIteration:
        raise PhylipFormatError("This file is empty.")
    n_seqs, seq_len = _validate_header(header)

    # All following lines should be ID+sequence. No blank lines are allowed.
    data = []
    for line in _line_generator(fh):
        data.append(_validate_line(line, seq_len))
    if len(data) != n_seqs:
        raise PhylipFormatError(
            "The number of sequences is not %s " % n_seqs +
            "as specified in the header.")
    return data


def _line_generator(fh):
    """Just remove linebreak characters and yield lines.
    """
    for line in fh:
        yield line.rstrip('\n')
