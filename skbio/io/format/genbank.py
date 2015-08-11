"""
Genbank format (:mod:`skbio.io.format.genbank`)
===============================================

.. currentmodule:: skbio.io.format.genbank

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
The following sections define the genbank formats.

<http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>

Feature Table Documentation:
http://www.insdc.org/files/feature_table.html
ftp://ftp.ncbi.nih.gov/genbank/docs/FTv10_3.html

    The International Nucleotide Sequence Database Collaboration (INSDC)
    between the DDBJ, EMBL, and GenBank.  These organisations all use the
    same "Feature Table" layout in their plain text flat file formats.

    However, the header and sequence sections of an EMBL file are very
    different in layout to those produced by GenBank/DDBJ.

Examples
--------

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
from future.builtins import range, zip
from six.moves import zip_longest
import re
import textwrap

import numpy as np
import pandas as pd
from datetime import datetime
from functools import partial

from skbio.io import create_format, GenbankFormatError
from skbio.io.registry import FileSentinel
from skbio.io.format._base import (_get_nth_sequence,
                                   _parse_fasta_like_header,
                                   _format_fasta_like_records, _line_generator,
                                   _too_many_blanks)
from skbio.util._misc import chunk_str
from skbio.sequence import Sequence, DNA, RNA, Protein
from pprint import pprint

genbank = create_format('genbank')

# date format in locus line of genbank record
_TIME_FORMAT = '%d-%b-%Y'
# This list is ordered
# used to read and write genbank file.
_HEADERS = ['LOCUS',
            'DEFINITION',
            'ACCESSION',
            'VERSION',
            'DBSOURCE',
            'DBLINK',
            'KEYWORDS',
            'SOURCE',
            'REFERENCE',
            'COMMENT',
            'FEATURES',
            'ORIGIN']


@genbank.sniffer()
def _genbank_sniffer(fh):
    # check the 1st line starts with 'LOCUS' and the last line ends with '//'
    if _too_many_blanks(fh, 5):
        return False, {}
    try:
        line = next(_line_generator(fh, skip_blanks=True))
    except StopIteration:
        return False, {}

    if line.startswith(_HEADERS[0]):
        for line in _line_generator(fh, skip_blanks=True):
            if line == '//':
                return True, {}
        return False, {}
    else:
        return False, {}


@genbank.reader(None)
def _genbank_to_generator(fh, constructor=None, **kwargs):
    for record in _parse_genbanks(fh):
        yield _construct(record, constructor, **kwargs)


@genbank.reader(Sequence)
def _genbank_to_biological_sequence(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_genbanks(fh), seq_num)
    return _construct(record, Sequence, **kwargs)


@genbank.reader(DNA)
def _genbank_to_dna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_genbanks(fh), seq_num)
    return _construct(record, DNA, **kwargs)


@genbank.reader(RNA)
def _genbank_to_rna(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_genbanks(fh), seq_num)
    return _construct(record, DNA, **kwargs).transcribe()


@genbank.reader(Protein)
def _genbank_to_protein(fh, seq_num=1, **kwargs):
    record = _get_nth_sequence(_parse_genbanks(fh), seq_num)
    return _construct(record, Protein, **kwargs)


@genbank.writer(None)
def _sequences_to_genbank(fh):
    pass


def _construct(record, constructor=None, **kwargs):
    seq, md, pmd = record
    if constructor is None:
        unit = md['LOCUS']['unit']
        if unit == 'bp':
            # RNA mol type has T instead of U;
            # so still read in as DNA
            constructor_ = DNA
        elif unit == 'aa':
            constructor_ = Protein
        else:
            constructor_ = Sequence
    else:
        constructor_ = constructor
    return constructor_(seq, metadata=md, positional_metadata=pmd, **kwargs)


def _parse_genbanks(fh):
    data_chunks = []
    for line in _line_generator(fh, skip_blanks=True, strip=False):
        if line.startswith('//'):
            yield _parse_single_genbank(data_chunks)
            data_chunks = []
        else:
            data_chunks.append(line)


def _parse_single_genbank(chunks):
    metadata = dict()
    metadata['REFERENCE'] = []
    sequence = ''
    # each section starts with a HEADER without indent.
    section_splitter = yield_section(lambda x: not x[0].isspace(), strip=False)
    for section in section_splitter(chunks):
        header = section[0].split(None, 1)[0]
        parser = globals().get('_parse_%s' % header.lower())
        if not parser:
            parser = _parse_section_default
        if header == 'FEATURES':
            parser = partial(
                parser, length=metadata['LOCUS']['size'])

        try:
            parsed = parser(section)
        except:
            raise GenbankFormatError(
                'Could not parse this section with %s:\n%s' %
                (parser, ''.join(section)))

        # reference can appear multiple times
        if header == 'REFERENCE':
            metadata[header].append(parsed)
        elif header == 'ORIGIN':
            sequence = parsed.upper()
        elif header == 'FEATURES':
            metadata[header] = parsed[0]
            positional_metadata = pd.concat(parsed[1], axis=1).to_sparse()
        else:
            metadata[header] = parsed
    return sequence, metadata, positional_metadata


def _parse_locus(lines):
    '''Parse the line LOCUS.

    Format:
    #    Positions  Contents
    #    ---------  --------
    #    00:06      LOCUS
    #    06:12      spaces
    #    12:??      Locus name
    #    ??:??      space
    #    ??:29      Length of sequence, right-justified
    #    29:33      space, bp/aa/rc, space
    #    33:41      molecule type (can be blank): DNA, ssDNA, dsRNA, tRNA, etc.
    #    41:42      space
    #    42:51      Blank (implies linear), linear or circular
    #    51:52      space
    #    52:55      The division code (e.g. BCT, VRL, INV)
    #    55:62      space
    #    62:73      Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)
    '''
    line = lines[0]
    pattern = (r'LOCUS'
               ' +([^\s]+)'
               ' +([0-9]+)'
               ' +(bp|aa|rc)'
               ' +(.*DNA|.*RNA)?'
               ' +(linear|circular)?'
               ' +(PRI|ROD|MAM|VRT|INV|PLN|BCT|VRL|PHG|'
               'SYN|UNA|EST|PAT|STS|GSS|HTG|HTC|ENV|CON)'
               ' +([0-9]{2}-[A-Z]{3}-[0-9]{4})')
    matches = re.match(pattern, line)
    res = dict()
    try:
        res = dict(zip(
            ['locus_name', 'size', 'unit',
             'mol_type', 'shape',
             'division', 'date'],
            matches.groups()))
    except:
        raise GenbankFormatError(
            "Could not parse the LOCUS line:\n%s" % line)

    # if res['shape'] is None:
    #     res['shape'] = 'linear'

    res['size'] = int(res['size'])
    res['date'] = datetime.strptime(res['date'], _TIME_FORMAT)
    return res


def _parse_reference(lines):
    '''Parse single REFERENCE field.
    '''
    res = {}
    # magic number 11: the non keyworded lines in REFERENCE
    # are at least indented with 11 spaces.
    feature_indent = ' ' * 11
    section_splitter = yield_section(
        lambda x: not x.startswith(feature_indent),
        skip_blanks=True, strip=False)
    for section in section_splitter(lines):
        label, data = _parse_section_default(
            section, join_delimitor=' ', return_label=True)
        res[label] = data
    return res


def _parse_source(lines):
    '''Parse SOURCE field.
    '''
    res = {}
    # magic number 11: the non keyworded lines in SOURCE
    # are at least indented with 11 spaces.
    feature_indent = ' ' * 11
    section_splitter = yield_section(
        lambda x: not x.startswith(feature_indent),
        skip_blanks=True, strip=False)
    # SOURCE line is not informative; skip it
    _, organism = list(section_splitter(lines))

    res['ORGANISM'] = organism[0].split(None, 1)[1].strip()
    res['taxonomy'] = ' '.join([i.strip() for i in organism[1:]])
    return res


def _parse_features(lines, length):
    '''Parse FEATURES field.
    '''
    features = []
    positional_metadata = []
    # skip the 1st FEATURES line
    if lines[0].startswith('FEATURES'):
        lines = lines[1:]
    # magic number 20: the lines following header of each feature
    # are at least indented with 20 spaces.
    feature_indent = ' ' * 20
    section_splitter = yield_section(
        lambda x: not x.startswith(feature_indent),
        skip_blanks=True, strip=False)
    for i, section in enumerate(section_splitter(lines)):
        feature, pmd = _parse_single_feature(section, length, i)
        features.append(feature)
        positional_metadata.append(pmd)
    return features, positional_metadata


def _parse_single_feature(lines, length, index):
    '''Parse a feature.

    Returns
    -------
    tuple
        Tuple of a dict of `metadata` and a pandas.Series of
        `positional_metadata` for the feature.

    '''
    feature = dict()
    feature['index_'] = index
    # each component of a feature starts with '/', except the 1st
    # component of location.
    section_splitter = yield_section(
        lambda x: x.startswith('/'), strip=True)
    first = True
    for section in section_splitter(lines):
        if first:
            # first section is the Location string
            first = False
            type, location = _parse_section_default(
                section, join_delimitor='', return_label=True)
            feature['type_'] = type
            feature['location'] = location
            loc, loc_pmd = _parse_loc_str(location, length)
            feature.update(loc)
        else:
            # following sections are Qualifiers
            k, v = _parse_section_default(
                section, label_delimitor='=',
                join_delimitor=' ', return_label=True)
            k = k[1:]
            # strip the quotes if it is quoted.
            # v could be empty
            if v and v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            # some Qualifiers can appear multiple times
            if k in feature:
                feature[k] = [feature[k]]
                feature[k].append(v)
            else:
                feature[k] = v
    return feature, loc_pmd


def _parse_loc_str(loc_str, length):
    '''Parse location string.

    Warning: This converts coordinates to 0-based from 1-based as
    in Genbank format.
    '''
    pmd = np.zeros(length, dtype=bool)
    res = {'rc_': False,
           'left_partial_': False,
           'right_partial_': False}
    items = re.split('[(),]+', loc_str)
    operators = ['join', 'complement']
    if 'complement' in items:
        res['rc_'] = True
    for i in items:
        i = i.strip()
        if i in operators or not i:
            continue
        elif '..' in i:  # span
            beg, end = i.split('..')
            if beg.startswith('<'):
                beg = beg[1:]
                res['left_partial_'] = True
            if end.startswith('>'):
                end = end[1:]
                res['right_partial_'] = True
            beg = int(beg)
            end = int(end)
            index = range(beg-1, end)
        elif i.isdigit():  # single base
            index = int(i) - 1
        else:
            raise GenbankFormatError(
                'Could not parse location string:\n%s' %
                loc_str)
        pmd[index] = True

    return res, pd.Series(pmd)


def _parse_origin(lines):
    '''Parse the ORIGIN section for sequence.
    '''
    sequence = []
    for line in lines:
        if line.startswith('ORIGIN'):
            continue
        # remove the number at the beg of each line
        items = line.split()
        sequence.append(''.join(items[1:]))
    return ''.join(sequence)


def _parse_section_default(
        lines, label_delimitor=None, join_delimitor=' ', return_label=False):
    '''Parse sections in default way.

    Do 2 things:
        1. split first line with label_delimitor for label
        2. join all the lines into one str with join_delimitor.
    '''
    data = []
    first = True
    label = None
    for line in lines:
        if first:
            try:
                items = line.split(label_delimitor, 1)
            except:
                GenbankFormatError('Could not split the line:\n%s', line)
            if len(items) == 2:
                label, section = items
            else:
                label = items[0]
                section = ""
            data.append(section)
            first = False
        else:
            data.append(line)
    data = join_delimitor.join(i.strip() for i in data)
    if return_label:
        return label, data
    else:
        return data


def yield_section(is_another_section, **kwargs):
    '''Returns function that returns successive sections from file.

    Parameters
    ----------
    is_another_section : callable
        It takes a string as input and return a boolean indicating
        a new section starts.
    kwargs : dict, optional
        Keyword arguments will be passed to `_line_generator`.

    Returns
    -------
    generator
        A generator to yield section one by one.
    '''
    def parser(lines):
        curr = []
        for line in _line_generator(lines, **kwargs):
            # if we find another, return the previous section
            if is_another_section(line):
                if curr:
                    yield curr
                    curr = []
            curr.append(line)
        # don't forget to return the last section in the file
        if curr:
            yield curr
    return parser
