#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from unittest import TestCase, main

from numpy import array

from skbio.parse.sequences.factory import load
from skbio.parse.sequences.iterator import FastaIterator


# copied from maths/stats/ordination/tests/test_ordination.py
def get_data_path(fn):
    """Return path to filename `fn` in the data folder."""
    path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(path, 'data', fn)
    return data_path


class SequenceLoadTests(TestCase):
    def setUp(self):
        self.fna1 = get_data_path('fna1.fasta')
        self.fna1gz = get_data_path('fna1.fna.gz')
        self.fq1 = get_data_path('fq1.fq')
        self.fq1gz = get_data_path('fq1.fastq.gz')
        self.qual1 = get_data_path('fna1.qual')
        self.noext = get_data_path('noextensionfasta')

    def test_single_files(self):
        """load should handle a single file, and can be gzipped"""
        it = load(self.fna1)
        obs = [rec.copy() for rec in it]
        exp = [{'Sequence': 'ATGC', 'SequenceID': 's1',
                'QualID': None, 'Qual': None},
               {'Sequence': 'AATTGG', 'SequenceID': 's2',
                'QualID': None, 'Qual': None}]
        self.assertEqual(obs, exp)
        it = load(self.fq1, phred_offset=64)
        obs = [rec.copy() for rec in it]
        exp = [{'Sequence': 'ATGC', 'SequenceID': 's1',
                'QualID': 's1', 'Qual': array([40, 40, 40, 40])},
               {'Sequence': 'AATTGG', 'SequenceID': 's2',
                'QualID': 's2', 'Qual': array([39, 39, 39, 39, 40, 40])}]
        for o, e in zip(obs, exp):
            self.assertEqual(o['Sequence'], e['Sequence'])
            self.assertEqual(o['SequenceID'], e['SequenceID'])
            self.assertEqual(o['QualID'], e['QualID'])
            self.assertTrue((o['Qual'] == e['Qual']).all())

        it = load(self.fna1gz)
        obs = [rec.copy() for rec in it]
        exp = [{'Sequence': 'ATGC', 'SequenceID': 's1',
                'QualID': None, 'Qual': None},
               {'Sequence': 'AATTGG', 'SequenceID': 's2',
                'QualID': None, 'Qual': None}]
        self.assertEqual(obs, exp)

        it = load(self.fq1gz, phred_offset=64)
        obs = [rec.copy() for rec in it]
        exp = [{'Sequence': 'ATGC', 'SequenceID': 's1',
                'QualID': 's1', 'Qual': array([40, 40, 40, 40])},
               {'Sequence': 'AATTGG', 'SequenceID': 's2',
                'QualID': 's2', 'Qual': array([39, 39, 39, 39, 40, 40])}]
        for o, e in zip(obs, exp):
            self.assertEqual(o['Sequence'], e['Sequence'])
            self.assertEqual(o['SequenceID'], e['SequenceID'])
            self.assertEqual(o['QualID'], e['QualID'])
            self.assertTrue((o['Qual'] == e['Qual']).all())

    def test_multiple_files(self):
        """load should handle multiple files of different types"""
        it = load([self.fq1, self.fna1], phred_offset=64)
        obs = [rec.copy() for rec in it]
        exp = [{'Sequence': 'ATGC', 'SequenceID': 's1',
                'QualID': 's1', 'Qual': array([40, 40, 40, 40])},
               {'Sequence': 'AATTGG', 'SequenceID': 's2',
                'QualID': 's2', 'Qual': array([39, 39, 39, 39, 40, 40])},
               {'Sequence': 'ATGC', 'SequenceID': 's1',
                'QualID': None, 'Qual': None},
               {'Sequence': 'AATTGG', 'SequenceID': 's2',
                'QualID': None, 'Qual': None}]

        o = obs[0]
        e = exp[0]
        self.assertEqual(o['Sequence'], e['Sequence'])
        self.assertEqual(o['SequenceID'], e['SequenceID'])
        self.assertEqual(o['QualID'], e['QualID'])
        self.assertTrue((o['Qual'] == e['Qual']).all())

        o = obs[1]
        e = exp[1]
        self.assertEqual(o['Sequence'], e['Sequence'])
        self.assertEqual(o['SequenceID'], e['SequenceID'])
        self.assertEqual(o['QualID'], e['QualID'])
        self.assertTrue((o['Qual'] == e['Qual']).all())

        o = obs[2]
        e = exp[2]
        self.assertEqual(o['Sequence'], e['Sequence'])
        self.assertEqual(o['SequenceID'], e['SequenceID'])
        self.assertEqual(o['QualID'], e['QualID'])
        self.assertEqual(o['Qual'], e['Qual'])

        o = obs[3]
        e = exp[3]
        self.assertEqual(o['Sequence'], e['Sequence'])
        self.assertEqual(o['SequenceID'], e['SequenceID'])
        self.assertEqual(o['QualID'], e['QualID'])
        self.assertEqual(o['Qual'], e['Qual'])

    def test_transform(self):
        """load should pass transform methods to the iterators"""
        def rev_f(st):
            st['Sequence'] = st['Sequence'][::-1]
            st['Qual'] = st['Qual'][::-1] if st['Qual'] is not None else None

        it = load([self.fq1gz, self.fna1], transform=rev_f, phred_offset=64)
        obs = [rec.copy() for rec in it]
        exp = [{'Sequence': 'CGTA', 'SequenceID': 's1',
                'QualID': 's1', 'Qual': array([40, 40, 40, 40])},
               {'Sequence': 'GGTTAA', 'SequenceID': 's2',
                'QualID': 's2', 'Qual': array([40, 40, 39, 39, 39, 39])},
               {'Sequence': 'CGTA', 'SequenceID': 's1',
                'QualID': None, 'Qual': None},
               {'Sequence': 'GGTTAA', 'SequenceID': 's2',
                'QualID': None, 'Qual': None}]

        o = obs[0]
        e = exp[0]
        self.assertEqual(o['Sequence'], e['Sequence'])
        self.assertEqual(o['SequenceID'], e['SequenceID'])
        self.assertEqual(o['QualID'], e['QualID'])
        self.assertTrue((o['Qual'] == e['Qual']).all())

        o = obs[1]
        e = exp[1]
        self.assertEqual(o['Sequence'], e['Sequence'])
        self.assertEqual(o['SequenceID'], e['SequenceID'])
        self.assertEqual(o['QualID'], e['QualID'])
        self.assertTrue((o['Qual'] == e['Qual']).all())

        o = obs[2]
        e = exp[2]
        self.assertEqual(o['Sequence'], e['Sequence'])
        self.assertEqual(o['SequenceID'], e['SequenceID'])
        self.assertEqual(o['QualID'], e['QualID'])
        self.assertEqual(o['Qual'], e['Qual'])

        o = obs[3]
        e = exp[3]
        self.assertEqual(o['Sequence'], e['Sequence'])
        self.assertEqual(o['SequenceID'], e['SequenceID'])
        self.assertEqual(o['QualID'], e['QualID'])
        self.assertEqual(o['Qual'], e['Qual'])

    def test_force_constructor(self):
        it = load([self.noext], constructor=FastaIterator)
        obs = [rec.copy() for rec in it]
        exp = [{'Sequence': 'AATTGG', 'SequenceID': 'seq1',
                'Qual': None, 'QualID': None},
               {'Sequence': 'ATATA', 'SequenceID': 'seq2',
                'Qual': None, 'QualID': None}]
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
