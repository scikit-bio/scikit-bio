# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import six

import unittest

from skbio import TabularMSA, Protein, DNA, RNA
from skbio.io import StockholmFormatError
from skbio.io.format.stockholm import (_stockholm_to_tabular_msa,
                                       _stockholm_sniffer)
from skbio.util import get_data_path


class TestStockholmSniffer(unittest.TestCase):
    def setUp(self):
        self.positives = [get_data_path(e) for e in [
            'stockholm_extensive',
            'stockholm_minimal',
            'stockholm_rna',
            'stockholm_runon_gf',
            'stockholm_duplicate_sequence_names',
            'stockholm_duplicate_gr',
            'stockholm_duplicate_gc',
            'stockholm_invalid_nonexistent_gr',
            'stockholm_invalid_nonexistent_gs',
            'stockholm_no_data',
            'stockholm_blank_lines',
            'stockholm_differing_gc_data_length',
            'stockholm_differing_gr_data_length',
            'stockholm_differing_seq_lengths',
            'stockholm_duplicate_sequence_names',
            'stockholm_duplicate_tree_ids',
            'stockholm_extensive_mixed',
            'stockholm_invalid_data_type',
            'stockholm_malformed_gf_line',
            'stockholm_malformed_gs_line',
            'stockholm_malformed_gr_line',
            'stockholm_malformed_gc_line',
            'stockholm_malformed_data_line',
            'stockholm_metadata_only',
            'stockholm_multiple_msa',
            'stockholm_multiple_trees',
            'stockholm_runon_gs',
            'stockholm_single_tree_with_id',
            'stockholm_single_tree_without_id',
            'stockholm_whitespace_only_lines'
            ]]

        self.negatives = [get_data_path(e) for e in [
            'stockholm_missing_header',
            'stockholm_missing_footer',
            'empty',
            'whitespace_only'
            ]]

    def test_positives(self):
        for fp in self.positives:
            self.assertEqual(_stockholm_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negatives:
            self.assertEqual(_stockholm_sniffer(fp), (False, {}))


class TestStockholmReader(unittest.TestCase):
    def test_stockholm_extensive(self):
        fp = get_data_path('stockholm_extensive')
        msa = _stockholm_to_tabular_msa(fp, constructor=Protein)
        exp = TabularMSA([Protein('MTCRAQLIAVPRASSLAE..AIACAQKM....'
                                  'RVSRVPVYERS',
                                  positional_metadata={'SA': list('9998877564'
                                                                  '53524252..'
                                                                  '55152525..'
                                                                  '..36463774'
                                                                  '777')}),
                          Protein('EVMLTDIPRLHINDPIMK..GFGMVINN....'
                                  '..GFVCVENDE',
                                  metadata={'OS': 'Bacillus subtilis'},
                                  positional_metadata={'SS': list('CCCCCCCHHHH'
                                                                  'HHHHHHH..HE'
                                                                  'EEEEEE....E'
                                                                  'EEEEEE'
                                                                  'EEEH')}),
                          Protein('EVMLTDIPRLHINDPIMK..GFGMVINN...'
                                  '...GFVCVENDE',
                                  positional_metadata={'AS': list('___________'
                                                                  '_____*_____'
                                                                  '___________'
                                                                  '________'
                                                                  '__'),
                                                       'IN': list('___________'
                                                                  '_1_________'
                                                                  '_____2_____'
                                                                  '_____0_'
                                                                  '___')})],
                         metadata={'ID': 'CBS', 'AC': 'PF00571',
                                   'AU': 'Bateman A', 'SQ': '67'},
                         positional_metadata={'SS_cons': list('CCCCCHHHHHHHH'
                                                              'HHHHH..EEEEEE'
                                                              'EE....EEEEEEE'
                                                              'EEEH')},
                         index=['O83071/192-246', 'O31698/88-139',
                                'O31699/88-139'])
        self.assertEqual(msa, exp)

    def test_stockholm_extensive_mixed(self):
        fp = get_data_path('stockholm_extensive_mixed')
        msa = _stockholm_to_tabular_msa(fp, constructor=Protein)
        exp = TabularMSA([Protein('MTCRAQLIAVPRASSLAE..AIACAQKM....'
                                  'RVSRVPVYERS',
                                  positional_metadata={'SA': list('9998877564'
                                                                  '53524252..'
                                                                  '55152525..'
                                                                  '..36463774'
                                                                  '777')}),
                          Protein('EVMLTDIPRLHINDPIMK..GFGMVINN....'
                                  '..GFVCVENDE',
                                  metadata={'OS': 'Bacillus subtilis'},
                                  positional_metadata={'SS': list('CCCCCCCHHHH'
                                                                  'HHHHHHH..HE'
                                                                  'EEEEEE....E'
                                                                  'EEEEEE'
                                                                  'EEEH')}),
                          Protein('EVMLTDIPRLHINDPIMK..GFGMVINN...'
                                  '...GFVCVENDE',
                                  positional_metadata={'AS': list('___________'
                                                                  '_____*_____'
                                                                  '___________'
                                                                  '________'
                                                                  '__'),
                                                       'IN': list('___________'
                                                                  '_1_________'
                                                                  '_____2_____'
                                                                  '_____0_'
                                                                  '___')})],
                         metadata={'ID': 'CBS', 'AC': 'PF00571',
                                   'AU': 'Bateman A', 'SQ': '67'},
                         positional_metadata={'SS_cons': list('CCCCCHHHHHHHH'
                                                              'HHHHH..EEEEEE'
                                                              'EE....EEEEEEE'
                                                              'EEEH')},
                         index=['O83071/192-246', 'O31698/88-139',
                                'O31699/88-139'])
        self.assertEqual(msa, exp)

    def test_stockholm_minimal(self):
        fp = get_data_path('stockholm_minimal')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([DNA('TGTGTCGCAGTTGTCGTTTG')], index=['0235244'])
        self.assertEqual(msa, exp)

    def test_stockholm_rna(self):
        fp = get_data_path('stockholm_rna')
        msa = _stockholm_to_tabular_msa(fp, constructor=RNA)
        exp = TabularMSA([RNA('AAGGGUUAUUUAUAUACUUU'),
                          RNA('UGCUAAGAGUGGGGAUGAUU'),
                          RNA('GCCACAACCGAUUAGAUAGA'),
                          RNA('UUAGAAACCGAUGGACCGAA')],
                         metadata={'AC': 'G2134T23', 'ID': 'ARD'},
                         positional_metadata=(
                         {'AC_cons': list('GGGACUGGACAUCUAUUCAG')}),
                         index=['RTC2231', 'RTF2124', 'RTH3322', 'RTB1512'])
        self.assertEqual(msa, exp)

    def test_stockholm_runon_gf(self):
        fp = get_data_path('stockholm_runon_gf')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([DNA('ACTGGTTCAATG')],
                         metadata={'CC': 'CBS domains are small intracellular'
                                         ' modules mostly found in 2 or four '
                                         'copies within a protein.'},
                         index=['GG1344'])
        self.assertEqual(msa, exp)

    def test_stockholm_runon_gs(self):
        fp = get_data_path('stockholm_runon_gs')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([DNA('ATCGTTCAGTG',
                              metadata={'AL': 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'})],
                         index=['seq1'])
        self.assertEqual(msa, exp)

    def test_stockholm_metadata_only(self):
        fp = get_data_path('stockholm_metadata_only')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([], metadata={'NM': 'Kestrel Gorlick',
                                       'DT': 'February 5th, 2016'})
        self.assertEqual(msa, exp)

    def test_stockholm_no_data(self):
        fp = get_data_path('stockholm_no_data')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([])
        self.assertEqual(msa, exp)

    def test_stockholm_with_blank_lines(self):
        fp = get_data_path('stockholm_blank_lines')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([], metadata={'AL': 'ABCD', 'NM': '1234'})
        self.assertEqual(msa, exp)

    def test_stockholm_with_whitespace_only_lines(self):
        fp = get_data_path('stockholm_whitespace_only_lines')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([], metadata={'AL': 'ABCD', 'NM': '1234'})
        self.assertEqual(msa, exp)

    def test_stockholm_single_tree_without_id(self):
        fp = get_data_path('stockholm_single_tree_without_id')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([], metadata={'NH': 'ABCD'})
        self.assertEqual(msa, exp)

    def test_stockholm_single_tree_with_id(self):
        fp = get_data_path('stockholm_single_tree_with_id')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([], metadata={'NH': {'tree1': 'ABCD'}})
        self.assertEqual(msa, exp)

    def test_stockholm_multiple_trees(self):
        fp = get_data_path('stockholm_multiple_trees')
        msa = _stockholm_to_tabular_msa(fp, constructor=DNA)
        exp = TabularMSA([], metadata={'NH': {'tree1': 'ABCD',
                                              'tree2': 'EFGH',
                                              'tree3': 'IJKL'}})
        self.assertEqual(msa, exp)

    def test_multiple_msa_file(self):
        fp = get_data_path('stockholm_multiple_msa')
        msa = _stockholm_to_tabular_msa(fp, constructor=RNA)
        exp = TabularMSA([RNA('AAGGGUUAUUUAUAUACUUU'),
                          RNA('UGCUAAGAGUGGGGAUGAUU'),
                          RNA('GCCACAACCGAUUAGAUAGA'),
                          RNA('UUAGAAACCGAUGGACCGAA')],
                         metadata={'AC': 'G2134T23', 'ID': 'ARD'},
                         positional_metadata=(
                         {'AC_cons': list('GGGACUGGACAUCUAUUCAG')}),
                         index=['RTC2231', 'RTF2124', 'RTH3322', 'RTB1512'])
        self.assertEqual(msa, exp)


    def test_stockholm_duplicate_tree_id_error(self):
        fp = get_data_path('stockholm_duplicate_tree_ids')
        with self.assertRaisesRegex(StockholmFormatError,
                                    'Tree.*tree1.*in file.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_nonexistent_gr_error(self):
        fp = get_data_path('stockholm_invalid_nonexistent_gr')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   'Markup line references.*RL1355.'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_nonexistent_gs_error(self):
        fp = get_data_path('stockholm_invalid_nonexistent_gs')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   'Markup line references.*AC14.'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_duplicate_sequence_names_error(self):
        fp = get_data_path('stockholm_duplicate_sequence_names')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   'Found multiple data.*ASR132.*supported by '
                                   'the reader.'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_duplicate_gr_error(self):
        fp = get_data_path('stockholm_duplicate_gr')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   'Found duplicate GR.*OS.*LFDR3.*supported '
                                   'by the reader.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_duplicate_gc_error(self):
        fp = get_data_path('stockholm_duplicate_gc')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   'Found duplicate GC.*SS_cons.*supported '
                                   'by the reader.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_empty_file_error(self):
        fp = get_data_path('empty')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   'File is empty.'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_missing_header_error(self):
        fp = get_data_path('stockholm_missing_header')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   r'File missing.*`# STOCKHOLM 1.0\\n`.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_missing_footer_error(self):
        fp = get_data_path('stockholm_missing_footer')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   'Final line.*only `//`.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_data_type_error(self):
        fp = get_data_path('stockholm_invalid_data_type')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                   "Unrecognized.*'#=GZ'"):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_malformed_gf_line_error(self):
        fp = get_data_path('stockholm_malformed_gf_line')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                    'Line only contains 2.*must contain 3.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_malformed_gs_line_error(self):
        fp = get_data_path('stockholm_malformed_gs_line')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                    'Line only contains 3.*must contain 4.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_malformed_gr_line_error(self):
        fp = get_data_path('stockholm_malformed_gr_line')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                    'Line only contains 2.*must contain 4.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_malformed_gc_line_error(self):
        fp = get_data_path('stockholm_malformed_gc_line')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                    'Line only contains 2.*must contain 3.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_malformed_data_line_error(self):
        fp = get_data_path('stockholm_malformed_data_line')
        with six.assertRaisesRegex(self, StockholmFormatError,
                                    'Line only contains 1.*must contain 2.'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

    def test_differing_sequence_lengths_error(self):
        fp = get_data_path('stockholm_differing_seq_lengths')
        with six.assertRaisesRegex(self, ValueError,
                                   'Each sequence.*11 != 10'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_differing_data_lengths_gr_error(self):
        fp = get_data_path('stockholm_differing_gr_data_length')
        with six.assertRaisesRegex(self, ValueError,
                                   'Number.*7.*(8).'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_differing_data_lengths_gc_error(self):
        fp = get_data_path('stockholm_differing_gc_data_length')
        with six.assertRaisesRegex(self, ValueError,
                                   'Number.*12.*(10).'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_no_constructor_error(self):
        fp = get_data_path('empty')
        with six.assertRaisesRegex(self, ValueError, 'Must.*parameter.'):
            _stockholm_to_tabular_msa(fp)

    def test_unsupported_constructor_error(self):
        fp = get_data_path('empty')
        with six.assertRaisesRegex(self, TypeError,
                                   '`constructor`.*`IUPACSequence`'):
            _stockholm_to_tabular_msa(fp, constructor=TabularMSA)

if __name__ == '__main__':
    unittest.main()
