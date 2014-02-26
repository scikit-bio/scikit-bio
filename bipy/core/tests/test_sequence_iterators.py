#!/usr/bin/env python

from StringIO import StringIO
from numpy import arange, array
from bipy.util.unit_test import TestCase, main
from bipy.core.sequence_iterators import (SequenceIterator, SequenceRecord,
        FastaIterator, FastqIterator)


class SeqIterTests(TestCase):
    def setUp(self):
        self.seq_ok = SequenceRecord('foo',
                                     'AATTGGCC',
                                     None,
                                     None,
                                     None,
                                     None)

        self.seqqual_ok = SequenceRecord('foo',
                                         'AATTGGCC',
                                         'foo',
                                         arange(8),
                                         None,
                                         None)

        self.seqqualbc_ok = SequenceRecord('foo',
                                           'AATTGGCC',
                                           'foo',
                                           arange(8),
                                           'foo',
                                           'aattcc')

        self.seq_bad = SequenceRecord('foo',
                                      'AATT  GGCC',
                                      None,
                                      None,
                                      None,
                                      None)

        self.seqqual_bad_id = SequenceRecord('foo',
                                             'AATTGGCC',
                                             'bar',
                                             arange(8),
                                             None,
                                             None)

        self.seqqual_bad_qual = SequenceRecord('foo',
                                               'AATTGGCC',
                                               'foo',
                                               arange(5),
                                               None,
                                               None)

        self.seqqualbc_bad_id = SequenceRecord('foo',
                                               'AATTGGCC',
                                               'foo',
                                               arange(8),
                                               'bar',
                                               'aattgg')

        self.seqbc_bad_id = SequenceRecord('foo',
                                           'AATTGGCC',
                                           None,
                                           None,
                                           'bar',
                                           'aattgg')

        self.rev_f = lambda x: x[::-1]

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

        wk.initialize_state(self.seqbc_bad_id)
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

        wk.initialize_state(self.seqbc_bad_id)
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

    def test_transform_seq(self):
        wk = SequenceIterator(['aattgg'], transform_seq=self.rev_f)

        wk.initialize_state(self.seq_ok)
        self.assertEqual(wk.state['Sequence'], self.seq_ok.Sequence)
        wk.transform_seq()
        self.assertEqual(wk.state['Sequence'], self.seq_ok.Sequence[::-1])

    def test_transform_qual(self):
        wk = SequenceIterator(['aattgg'], transform_qual=self.rev_f)

        wk.initialize_state(self.seqqual_ok)
        self.assertTrue((wk.state['Qual'] == self.seqqual_ok.Qual).all())
        wk.transform_qual()
        self.assertTrue((wk.state['Qual'] == self.seqqual_ok.Qual[::-1]).all())

    def test_transform_bc(self):
        wk = SequenceIterator(['aattgg'], transform_bc=self.rev_f)

        wk.initialize_state(self.seqqualbc_ok)
        self.assertEqual(wk.state['Barcode'], self.seqqualbc_ok.Barcode)
        wk.transform_bc()
        self.assertEqual(wk.state['Barcode'], self.seqqualbc_ok.Barcode[::-1])

class FastaTests(TestCase):
    def setUp(self):
        self.fastas = [StringIO(fasta1), StringIO(fasta2), StringIO(fasta3)]
        self.quals = [StringIO(qual1), StringIO(qual2), StringIO(qual3)]
        self.bcs = [StringIO(bc)]

        self.bad_qual_val = [StringIO(qual1), StringIO(qual_bad_val),
                             StringIO(qual3)]
        self.bad_qual_id = [StringIO(qual1), StringIO(qual_bad_id),
                            StringIO(qual3)]

        self.bad_bc = [StringIO(bc_bad)]

    def test_fasta_gen(self):
        wk = FastaIterator(seq=self.fastas)
        gen = wk()

        exp1 = {'SequenceID': '1', 'Sequence': 'aattggcc', 'Qual': None,
                'QualID': None, 'BarcodeID': None, 'Barcode': None}
        exp2 = {'SequenceID': '2', 'Sequence': 'aattaatt', 'Qual': None,
                'QualID': None, 'BarcodeID': None, 'Barcode': None}
        exp3 = {'SequenceID': '3', 'Sequence': 'atat', 'Qual': None,
                'QualID': None, 'BarcodeID': None, 'Barcode': None}
        exp4 = {'SequenceID': '4', 'Sequence': 'attatt', 'Qual': None,
                'QualID': None, 'BarcodeID': None, 'Barcode': None}
        exp5 = {'SequenceID': '5', 'Sequence': 'ggccc', 'Qual': None,
                'QualID': None, 'BarcodeID': None, 'Barcode': None}

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
                'Qual': arange(1, 9), 'QualID': '1', 'BarcodeID': None,
                'Barcode': None}
        exp2 = {'SequenceID': '2', 'Sequence': 'aattaatt', 'QualID': '2',
                'Qual': arange(1, 9)[::-1], 'BarcodeID': None, 'Barcode': None}
        exp3 = {'SequenceID': '3', 'Sequence': 'atat', 'Qual': arange(1, 5),
                'QualID': '3', 'BarcodeID': None, 'Barcode': None}
        exp4 = {'SequenceID': '4', 'Sequence': 'attatt', 'Qual': arange(1, 7),
                'QualID': '4', 'BarcodeID': None, 'Barcode': None}
        exp5 = {'SequenceID': '5', 'Sequence': 'ggccc', 'Qual': arange(1, 6),
                'QualID': '5', 'BarcodeID': None, 'Barcode': None}

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

    def test_fasta_qual_bc(self):
        wk = FastaIterator(seq=self.fastas, qual=self.quals, bc=self.bcs)
        gen = wk()

        exp1 = {'SequenceID': '1', 'Sequence': 'aattggcc',
                'Qual': arange(1, 9), 'QualID': '1', 'BarcodeID': '1',
                'Barcode': 'at'}
        exp2 = {'SequenceID': '2', 'Sequence': 'aattaatt', 'QualID': '2',
                'Qual': arange(1, 9)[::-1], 'BarcodeID': '2',
                'Barcode': 'atat'}
        exp3 = {'SequenceID': '3', 'Sequence': 'atat', 'Qual': arange(1, 5),
                'QualID': '3', 'BarcodeID': '3', 'Barcode': 'atatat'}
        exp4 = {'SequenceID': '4', 'Sequence': 'attatt', 'Qual': arange(1, 7),
                'QualID': '4', 'BarcodeID': '4', 'Barcode': 'ata'}
        exp5 = {'SequenceID': '5', 'Sequence': 'ggccc', 'Qual': arange(1, 6),
                'QualID': '5', 'BarcodeID': '5', 'Barcode': 'gtgt'}

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

    def test_fasta_bc(self):
        wk = FastaIterator(seq=self.fastas, bc=self.bcs)
        gen = wk()

        exp1 = {'SequenceID': '1', 'Sequence': 'aattggcc',
                'Qual': None, 'QualID': None, 'BarcodeID': '1',
                'Barcode': 'at'}
        exp2 = {'SequenceID': '2', 'Sequence': 'aattaatt', 'QualID': None,
                'Qual': None, 'BarcodeID': '2',
                'Barcode': 'atat'}
        exp3 = {'SequenceID': '3', 'Sequence': 'atat', 'Qual': None,
                'QualID': None, 'BarcodeID': '3', 'Barcode': 'atatat'}
        exp4 = {'SequenceID': '4', 'Sequence': 'attatt', 'Qual': None,
                'QualID': None, 'BarcodeID': '4', 'Barcode': 'ata'}
        exp5 = {'SequenceID': '5', 'Sequence': 'ggccc', 'Qual': None,
                'QualID': None, 'BarcodeID': '5', 'Barcode': 'gtgt'}

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

    def test_fasta_bad_bc(self):
        wk = FastaIterator(seq=self.fastas, bc=self.bad_bc)
        gen = wk()

        # default behavior is to sliently ignore
        exp_ids = ['1','2','4','5']
        obs_ids = [r['SequenceID'] for r in gen]

        self.assertEqual(obs_ids, exp_ids)

class FastqTests(TestCase):
    def setUp(self):
        self.fastqs = [StringIO(fastq1), StringIO(fastq2)]
        self.bcs = [StringIO(fq_bc1), StringIO(fq_bc2), StringIO(fq_bc3)]

    def test_fastq_gen(self):
        wk = FastqIterator(seq=self.fastqs)
        gen = wk()

        exp1 = {'SequenceID': '1', 'Sequence': 'atat', 'QualID': '1',
                'Qual': array([1, 2, 3, 4]), 'BarcodeID': None,
                'Barcode': None}
        exp2 = {'SequenceID': '2', 'Sequence': 'atgc', 'QualID': '2',
                'Qual': array([2, 3, 4, 5]), 'BarcodeID': None,
                'Barcode': None}
        exp3 = {'SequenceID': '3', 'Sequence': 'taa', 'QualID': '3',
                'Qual': array([5, 6, 7]), 'BarcodeID': None,
                'Barcode': None}

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

    def test_fastq_bc_gen(self):
        wk = FastqIterator(seq=self.fastqs, bc=self.bcs)
        gen = wk()

        exp1 = {'SequenceID': '1', 'Sequence': 'atat', 'QualID': '1',
                'Qual': array([1, 2, 3, 4]), 'BarcodeID': '1',
                'Barcode': 'at'}
        exp2 = {'SequenceID': '2', 'Sequence': 'atgc', 'QualID': '2',
                'Qual': array([2, 3, 4, 5]), 'BarcodeID': '2',
                'Barcode': 'atg'}
        exp3 = {'SequenceID': '3', 'Sequence': 'taa', 'QualID': '3',
                'Qual': array([5, 6, 7]), 'BarcodeID': '3',
                'Barcode': 'gg'}

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

bc = """>1
at
>2
atat
>3
atatat
>4
ata
>5
gtgt
"""

bc_bad = """>1
at
>2
atat
>badid
atatat
>4
ata
>5
gtgt
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

fq_bc1 = """@1
at
+
AA
"""

fq_bc2 = """@2
atg
+
AA
"""

fq_bc3 = """@3
gg
+
AA
"""

if __name__ == '__main__':
    main()
