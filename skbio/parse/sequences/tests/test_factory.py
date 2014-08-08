#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from numpy import array

from skbio import FastaIterator
from skbio.parse.sequences import load
from skbio.parse.sequences.factory import (
    _open_or_none, _is_single_iterator_type)
from skbio.util.testing import get_data_path


class SequenceLoadTests(TestCase):
    def setUp(self):
        self.fna1 = get_data_path('fna1.fasta')
        self.fna1gz = get_data_path('fna1.fna.gz')
        self.fq1 = get_data_path('fq1.fq')
        self.fq1gz = get_data_path('fq1.fastq.gz')
        self.qual1 = get_data_path('fna1.qual')
        self.noext = get_data_path('noextensionfasta')
        self.qs1 = get_data_path('qs1.qseq')
        self.qs1gz = get_data_path('qs1.qseq.gz')

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

        it = load(self.qs1, phred_offset=64)
        obs = [rec.copy() for rec in it]
        exp = [{'Qual': array([2, 27, 27, 27]),
                'QualID': 'CRESSIA_242:1:2204:1453:1918#0/1',
                'Sequence': 'TTAA',
                'SequenceID': 'CRESSIA_242:1:2204:1453:1918#0/1'},
               {'Qual': array([2, 2, 2, 2],),
                'QualID': 'CRESSIA_242:1:2204:1491:1920#0/1',
                'Sequence': 'AAAA',
                'SequenceID': 'CRESSIA_242:1:2204:1491:1920#0/1'}]
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

        it = load(self.qs1gz, phred_offset=64)
        obs = [rec.copy() for rec in it]
        exp = [{'Qual': array([2, 27, 27, 27]),
                'QualID': 'CRESSIA_242:1:2204:1453:1918#0/1',
                'Sequence': 'TTAA',
                'SequenceID': 'CRESSIA_242:1:2204:1453:1918#0/1'},
               {'Qual': array([2, 2, 2, 2],),
                'QualID': 'CRESSIA_242:1:2204:1491:1920#0/1',
                'Sequence': 'AAAA',
                'SequenceID': 'CRESSIA_242:1:2204:1491:1920#0/1'}]
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

    def test_no_seqs(self):
        for null in ('', [], (), None):
            with self.assertRaises(ValueError):
                load(null)

    def test_unknown_filetype(self):
        with self.assertRaises(IOError):
            load('seqs.mpeg')

    def test_file_path_does_not_exist(self):
        with self.assertRaises(IOError):
            load('this-seqs-file-had-better-not-exist-or-this-test-will-'
                 'fail.fna')

    def test_multiple_types_fasta_fastq_qual(self):
        with self.assertRaises(ValueError):
            load([self.fna1, self.fq1], qual=self.qual1)

    def test_open_or_none_no_opener(self):
        obs = _open_or_none(None, self.fna1)
        self.assertTrue(obs is None)

    def test_open_or_none_opener_error(self):
        def bogus_opener(f):
            raise IOError('hahaha')

        with self.assertRaises(IOError):
            _open_or_none(bogus_opener, self.fna1)

    def test_is_single_iterator_type_null_case(self):
        self.assertTrue(_is_single_iterator_type([]))


if __name__ == '__main__':
    main()
