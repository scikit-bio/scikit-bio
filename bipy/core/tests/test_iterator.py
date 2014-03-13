#!/usr/bin/env python

from StringIO import StringIO
from numpy import arange, array
from bipy.util.unit_test import TestCase, main
from bipy.core.iterator import (SequenceIterator, SequenceRecord,
        FastaIterator, FastqIterator)


class SeqIterTests(TestCase):
    def setUp(self):
        self.seq_ok = SequenceRecord('foo',
                                     'AATTGGCC',
                                     None,
                                     None)

        self.seqqual_ok = SequenceRecord('foo',
                                         'AATTGGCC',
                                         'foo',
                                         arange(8))

        self.seq_bad = SequenceRecord('foo',
                                      'AATT  GGCC',
                                      None,
                                      None)

        self.seqqual_bad_id = SequenceRecord('foo',
                                             'AATTGGCC',
                                             'bar',
                                             arange(8))

        self.seqqual_bad_qual = SequenceRecord('foo',
                                               'AATTGGCC',
                                               'foo',
                                               arange(5))


        def rev_f(st):
            st['Sequence']= st['Sequence'][::-1]
            st['Qual'] = st['Qual'][::-1] if st['Qual'] is not None else None

        self.rev_f = rev_f

    def test_validate_ids_true(self):
        wk = SequenceIterator(['aattgg'], valid_id=True)

        wk.initialize_state(self.seq_ok)
        wk.validate_ids()
        self.assertFalse(wk.failed)

        wk.initialize_state(self.seqqual_ok)
        wk.validate_ids()
        self.assertFalse(wk.failed)

        wk.initialize_state(self.seqqual_bad_id)
        wk.validate_ids()
        self.assertTrue(wk.failed)

    def test_validate_ids_false(self):
        wk = SequenceIterator(['aattgg'], valid_id=False)

        wk.initialize_state(self.seq_ok)
        wk.validate_ids()
        self.assertFalse(wk.failed)

        wk.initialize_state(self.seqqual_ok)
        wk.validate_ids()
        self.assertFalse(wk.failed)

        wk.initialize_state(self.seqqual_bad_id)
        wk.validate_ids()
        self.assertFalse(wk.failed)

    def test_validate_lengths_true(self):
        wk = SequenceIterator(['aattgg'], valid_length=True)

        wk.initialize_state(self.seq_ok)
        wk.valid_lengths()
        self.assertFalse(wk.failed)

        wk.initialize_state(self.seqqual_ok)
        wk.valid_lengths()
        self.assertFalse(wk.failed)

        wk.initialize_state(self.seqqual_bad_qual)
        wk.valid_lengths()
        self.assertTrue(wk.failed)

    def test_validate_lengths_false(self):
        wk = SequenceIterator(['aattgg'], valid_length=False)

        wk.initialize_state(self.seq_ok)
        wk.valid_lengths()
        self.assertFalse(wk.failed)

        wk.initialize_state(self.seqqual_ok)
        wk.valid_lengths()
        self.assertFalse(wk.failed)

        wk.initialize_state(self.seqqual_bad_qual)
        wk.valid_lengths()
        self.assertFalse(wk.failed)

    def test_transform(self):
        wk = SequenceIterator(['aattgg'], transform=self.rev_f)

        wk.initialize_state(self.seqqual_ok)
        self.assertEqual(wk.state['Sequence'], self.seqqual_ok.Sequence)
        wk.transform()
        self.assertEqual(wk.state['Sequence'], self.seqqual_ok.Sequence[::-1])
        self.assertTrue((wk.state['Qual'] == self.seqqual_ok.Qual[::-1]).all())


class FastaTests(TestCase):
    def setUp(self):
        self.fastas = [StringIO(fasta1), StringIO(fasta2), StringIO(fasta3)]
        self.quals = [StringIO(qual1), StringIO(qual2), StringIO(qual3)]

        self.bad_qual_val = [StringIO(qual1), StringIO(qual_bad_val),
                             StringIO(qual3)]
        self.bad_qual_id = [StringIO(qual1), StringIO(qual_bad_id),
                            StringIO(qual3)]

    def test_fasta_gen(self):
        wk = FastaIterator(seq=self.fastas)
        gen = wk()

        exp1 = {'SequenceID': '1', 'Sequence': 'aattggcc', 'Qual': None,
                'QualID': None}
        exp2 = {'SequenceID': '2', 'Sequence': 'aattaatt', 'Qual': None,
                'QualID': None}
        exp3 = {'SequenceID': '3', 'Sequence': 'atat', 'Qual': None,
                'QualID': None}
        exp4 = {'SequenceID': '4', 'Sequence': 'attatt', 'Qual': None,
                'QualID': None}
        exp5 = {'SequenceID': '5', 'Sequence': 'ggccc', 'Qual': None,
                'QualID': None}

        obs1 = gen.next()
        self.assertEqual(obs1, exp1)
        self.assertFalse(wk.failed)

        obs2 = gen.next()
        self.assertEqual(obs2, exp2)
        self.assertFalse(wk.failed)

        obs3 = gen.next()
        self.assertEqual(obs3, exp3)
        self.assertFalse(wk.failed)

        obs4 = gen.next()
        self.assertEqual(obs4, exp4)
        self.assertFalse(wk.failed)

        obs5 = gen.next()
        self.assertEqual(obs5, exp5)
        self.assertFalse(wk.failed)

    def test_fasta_qual(self):
        wk = FastaIterator(seq=self.fastas, qual=self.quals)
        gen = wk()

        exp1 = {'SequenceID': '1', 'Sequence': 'aattggcc',
                'Qual': arange(1, 9), 'QualID': '1'}
        exp2 = {'SequenceID': '2', 'Sequence': 'aattaatt', 'QualID': '2',
                'Qual': arange(1, 9)[::-1]}
        exp3 = {'SequenceID': '3', 'Sequence': 'atat', 'Qual': arange(1, 5),
                'QualID': '3'}
        exp4 = {'SequenceID': '4', 'Sequence': 'attatt', 'Qual': arange(1, 7),
                'QualID': '4'}
        exp5 = {'SequenceID': '5', 'Sequence': 'ggccc', 'Qual': arange(1, 6),
                'QualID': '5'}

        obs1 = gen.next()
        self.assertTrue((obs1['Qual'] == exp1['Qual']).all())
        obs1.pop('Qual'); exp1.pop('Qual')
        self.assertEqual(obs1, exp1)
        self.assertFalse(wk.failed)

        obs2 = gen.next()
        self.assertTrue((obs2['Qual'] == exp2['Qual']).all())
        obs2.pop('Qual'); exp2.pop('Qual')
        self.assertEqual(obs2, exp2)
        self.assertFalse(wk.failed)

        obs3 = gen.next()
        self.assertTrue((obs3['Qual'] == exp3['Qual']).all())
        obs3.pop('Qual'); exp3.pop('Qual')
        self.assertEqual(obs3, exp3)
        self.assertFalse(wk.failed)

        obs4 = gen.next()
        self.assertTrue((obs4['Qual'] == exp4['Qual']).all())
        obs4.pop('Qual'); exp4.pop('Qual')
        self.assertEqual(obs4, exp4)
        self.assertFalse(wk.failed)

        obs5 = gen.next()
        self.assertTrue((obs5['Qual'] == exp5['Qual']).all())
        obs5.pop('Qual'); exp5.pop('Qual')
        self.assertEqual(obs5, exp5)
        self.assertFalse(wk.failed)

    def test_fasta_badqual_val(self):
        wk = FastaIterator(seq=self.fastas, qual=self.bad_qual_val)
        gen = wk()

        # default behavior is to sliently ignore
        exp_ids = ['1','2','4','5']
        obs_ids = [r['SequenceID'] for r in gen]

        self.assertEqual(obs_ids, exp_ids)

    def test_fasta_badqual_id(self):
        wk = FastaIterator(seq=self.fastas, qual=self.bad_qual_id)
        gen = wk()

        # default behavior is to sliently ignore
        exp_ids = ['1','2','4','5']
        obs_ids = [r['SequenceID'] for r in gen]

        self.assertEqual(obs_ids, exp_ids)


class FastqTests(TestCase):
    def setUp(self):
        self.fastqs = [StringIO(fastq1), StringIO(fastq2)]

    def test_fastq_gen(self):
        wk = FastqIterator(seq=self.fastqs)
        gen = wk()

        exp1 = {'SequenceID': '1', 'Sequence': 'atat', 'QualID': '1',
                'Qual': array([1, 2, 3, 4])}
        exp2 = {'SequenceID': '2', 'Sequence': 'atgc', 'QualID': '2',
                'Qual': array([2, 3, 4, 5])}
        exp3 = {'SequenceID': '3', 'Sequence': 'taa', 'QualID': '3',
                'Qual': array([5, 6, 7])}

        obs1 = gen.next()
        self.assertTrue((obs1['Qual'] == exp1['Qual']).all())
        obs1.pop('Qual'); exp1.pop('Qual')
        self.assertEqual(obs1, exp1)

        obs2 = gen.next()
        self.assertTrue((obs2['Qual'] == exp2['Qual']).all())
        obs2.pop('Qual'); exp2.pop('Qual')
        self.assertEqual(obs2, exp2)

        obs3 = gen.next()
        self.assertTrue((obs3['Qual'] == exp3['Qual']).all())
        obs3.pop('Qual'); exp3.pop('Qual')
        self.assertEqual(obs3, exp3)

fasta1 = """>1
aattggcc
>2
aattaatt
"""

fasta2 = """>3
atat
"""

fasta3 = """>4
attatt
>5
ggccc
"""

qual1 = """>1
1 2 3 4 5 6 7 8
>2
8 7 6 5 4 3 2 1
"""

qual2 = """>3
1 2 3 4
"""

qual3 = """>4
1 2 3 4 5 6
>5
1 2 3 4 5
"""

qual_bad_val = """>3
1 2
"""

qual_bad_id = """>asdasd
1 2 3 4
"""

fastq1 = """@1
atat
+
ABCD
@2
atgc
+
BCDE
"""

fastq2 = """@3
taa
+
EFG
"""

if __name__ == '__main__':
    main()
