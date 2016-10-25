# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

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
