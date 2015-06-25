"""
PHYLIP multiple sequence alignment format (:mod:`skbio.io.formats.phylip`)
==========================================================================

.. currentmodule:: skbio.io.formats.phylip

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
**Has Sniffer: No**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|No    |Yes   |:mod:`skbio.alignment.Alignment`                               |
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

.. note:: scikit-bio currently only supports writing strict, sequential
   PHYLIP-formatted files from an ``skbio.alignment.Alignment``. It does not
   yet support reading PHYLIP-formatted files, nor does it support relaxed or
   interleaved PHYLIP formats.

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

.. note:: While not explicitly stated in the original PHYLIP format
   description, scikit-bio only supports writing unique sequence identifiers
   (i.e., duplicates are not allowed). Uniqueness is required because an
   ``skbio.alignment.Alignment`` cannot be created with duplicate IDs.

   scikit-bio supports the empty string (``''``) as a valid sequence ID. An
   empty ID will be padded with 10 spaces.

Sequence characters immediately follow the sequence ID. They *must* start at
the 11th character in the line, as the first 10 characters are reserved for the
sequence ID. While PHYLIP format does not explicitly restrict the set of
supported characters that may be used to represent a sequence, the original
format description [3]_ specifies the IUPAC nucleic acid lexicon for DNA or RNA
sequences, and the IUPAC protein lexicon for protein sequences. The original
PHYLIP specification uses ``-`` as a gap character, though older versions also
supported ``.``. The sequence characters may contain optional spaces (e.g., to
improve readability), and both upper and lower case characters are supported.

.. note:: scikit-bio will write a PHYLIP-formatted file even if the alignment's
   sequence characters are not valid IUPAC characters. This differs from the
   PHYLIP specification, which states that a PHYLIP-formatted file can only
   contain valid IUPAC characters. To check whether all characters are valid
   before writing, the user can call ``Alignment.is_valid()``.

   Since scikit-bio supports both ``-`` and ``.`` as gap characters (e.g., in
   ``skbio.alignment.Alignment``), both are supported when writing a
   PHYLIP-formatted file.

   When writing a PHYLIP-formatted file, scikit-bio will split up each sequence
   into chunks that are 10 characters long. Each chunk will be separated by a
   single space. The sequence will always appear on a single line (sequential
   format). It will *not* be wrapped across multiple lines. Sequences are
   chunked in this manner for improved readability, and because most example
   PHYLIP files are chunked in a similar way (e.g., see the example file
   above). Note that this chunking is not required by the PHYLIP format.

Examples
--------
Let's create an alignment with three DNA sequences of equal length:

>>> from skbio import Alignment, DNA
>>> seqs = [DNA('ACCGTTGTA-GTAGCT', metadata={'id':'seq1'}),
...         DNA('A--GTCGAA-GTACCT', metadata={'id':'sequence-2'}),
...         DNA('AGAGTTGAAGGTATCT', metadata={'id':'3'})]
>>> aln = Alignment(seqs)
>>> aln
<Alignment: n=3; mean +/- std length=16.00 +/- 0.00>

Now let's write the alignment to file in PHYLIP format, and take a look at the
output:

>>> from io import StringIO
>>> fh = StringIO()
>>> print(aln.write(fh, format='phylip').getvalue())
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

If the sequence IDs in an alignment surpass the 10-character limit, an error
will be raised when we try to write a PHYLIP file:

>>> long_id_seqs = [DNA('ACCGT', metadata={'id':'seq1'}),
...                 DNA('A--GT', metadata={'id':'long-sequence-2'}),
...                 DNA('AGAGT', metadata={'id':'seq3'})]
>>> long_id_aln = Alignment(long_id_seqs)
>>> fh = StringIO()
>>> long_id_aln.write(fh, format='phylip')
Traceback (most recent call last):
    ...
PhylipFormatError: Alignment can only be written in PHYLIP format if all \
sequence IDs have 10 or fewer characters. Found sequence with ID \
'long-sequence-2' that exceeds this limit. Use Alignment.update_ids to assign \
shorter IDs.
>>> fh.close()

One way to work around this is to update the IDs to be shorter. The recommended
way of accomplishing this is via ``Alignment.update_ids``, which provides a
flexible way of creating a new ``Alignment`` with updated IDs. For example, to
remap each of the IDs to integer-based IDs:

>>> short_id_aln, _ = long_id_aln.update_ids()
>>> short_id_aln.ids()
['1', '2', '3']

We can now write the new alignment in PHYLIP format:

>>> fh = StringIO()
>>> print(short_id_aln.write(fh, format='phylip').getvalue())
3 5
1         ACCGT
2         A--GT
3         AGAGT
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

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from skbio.alignment import Alignment
from skbio.io import create_format, PhylipFormatError
from skbio.io.formats._base import _chunk_str

phylip = create_format('phylip')

@phylip.writer(Alignment)
def _alignment_to_phylip(obj, fh):

    if obj.is_empty():
        raise PhylipFormatError(
            "Alignment can only be written in PHYLIP format if there is at "
            "least one sequence in the alignment.")

    sequence_length = obj.sequence_length()
    if sequence_length == 0:
        raise PhylipFormatError(
            "Alignment can only be written in PHYLIP format if there is at "
            "least one position in the alignment.")

    chunk_size = 10
    for id_ in obj.ids():
        if len(id_) > chunk_size:
            raise PhylipFormatError(
                "Alignment can only be written in PHYLIP format if all "
                "sequence IDs have %d or fewer characters. Found sequence "
                "with ID '%s' that exceeds this limit. Use "
                "Alignment.update_ids to assign shorter IDs." %
                (chunk_size, id_))

    sequence_count = obj.sequence_count()
    fh.write('{0:d} {1:d}\n'.format(sequence_count, sequence_length))

    fmt = '{0:%d}{1}\n' % chunk_size
    for seq in obj:
        chunked_seq = _chunk_str(str(seq), chunk_size, ' ')
        fh.write(fmt.format(seq.metadata['id'], chunked_seq))
