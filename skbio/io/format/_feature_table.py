# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

'''INSDC defined a common "feature table" specification for
DDBJ/ENA/GenBank formats It is described
`here<http://www.ebi.ac.uk/ena/WebFeat/>`_. The feature types and
their qualifiers are also described
`here<http://www.insdc.org/documents/feature_table.html>`_

This table explains how you store the metadata for each individual
interval feature. This tries to bring GenBank and GFF3 onto the same
common vocabularies to describe a interval feature.

+-----------+-----------+-----------+---------+------------------------------+
|   INSDC   |   GFF3    |    Key    |  Value  |         Description          |
|  feature  |columns or |  stored   |  type   |                              |
|   table   |attributes |           | stored  |                              |
+-----------+-----------+-----------+---------+------------------------------+
|    N/A    |  source   |  source   |   str   |the algorithm or experiment   |
|           |(column 2) |           |         |used to generate this feature |
+-----------+-----------+-----------+---------+------------------------------+
|feature key|   type    |   type    |   str   |the type of the feature       |
|           |(column 3) |           |         |                              |
+-----------+-----------+-----------+---------+------------------------------+
|    N/A    |   score   |   score   |  float  |the score of the feature      |
|           |(column 6) |           |         |                              |
+-----------+-----------+-----------+---------+------------------------------+
|    N/A    |  strand   |  strand   |   str   |the strand of the feature. +  |
|           |(column 7) |           |         |for positive strand, - for    |
|           |           |           |         |minus strand, and . for       |
|           |           |           |         |features that are not         |
|           |           |           |         |stranded. In addition, ?  can |
|           |           |           |         |be used for features whose    |
|           |           |           |         |strandedness is relevant, but |
|           |           |           |         |unknown.                      |
+-----------+-----------+-----------+---------+------------------------------+
|codon_start|   phase   |   phase   | int (0, |the offset at which the first |
|           |(column 8) |           |   1,    |complete codon of a coding    |
|           |           |           |  or 2)  |feature can be found, relative|
|           |           |           |         |to the first base of that     |
|           |           |           |         |feature. It is 0, 1, or 2 in  |
|           |           |           |         |GFF3 or 1, 2, or 3 in GenBank |
+-----------+-----------+-----------+---------+------------------------------+
|  db_xref  |  Dbxref   |  db_xref  | list of |A database cross reference    |
|           |           |           |   str   |                              |
+-----------+-----------+-----------+---------+------------------------------+
| locus_tag |    ID     | locus_tag |   str   |a submitter-supplied,         |
|           |           |           |         |systematic, stable identifier |
|           |           |           |         |for a gene and its associated |
|           |           |           |         |features, used for tracking   |
|           |           |           |         |purposes                      |
+-----------+-----------+-----------+---------+------------------------------+
|   note    |   Note    |   note    |   str   |any comment or additional     |
|           |           |           |         |information                   |
+-----------+-----------+-----------+---------+------------------------------+
|translation|    N/A    |translation|   str   |the protein sequence for CDS  |
|           |           |           |         |features                      |
+-----------+-----------+-----------+---------+------------------------------+

'''

def _vocabulary_change(format='genbank', read_in=True):
    '''Return a dict that converts between memory and output format.'''
    convert = {'phase': {'genbank': 'codon_start'},
               'db_xref': {'gff3': 'Dbxref'},
               'locus_tag': {'gff3': 'ID'},
               'note': {'gff3': 'Note'}}
    if read_in:
        return {v[format]: k for k, v in convert.items() if format in v}
    else:
        return {k: v[format] for k, v in convert.items() if format in v}


def _vocabulary_skip(format='genbank'):
    '''Return a list of keys that should be skipped when output to disk
    for the specified format.'''
    skip = {'source': ('genbank'),
            'type': ('genbank', 'gff3'),
            'translation': ('gff3'),
            'strand': ('genbank')}
    return [k for k, v in skip.items() if format in v]


def _parse_feature_table(lines):
    '''parse DDBJ/ENA/GenBank Feature Table.'''
