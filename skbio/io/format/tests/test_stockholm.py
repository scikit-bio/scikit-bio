# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from six import assertRaisesRegex

import unittest

from skbio import TabularMSA
from skbio.io import StockholmFormatError
from skbio.io.format.stockholm import (_stockholm_to_tabular_msa,
                                       _stockholm_sniffer)
from skbio import Protein, DNA, RNA
from skbio.util import get_data_path


class TestStockholmSniffer(unittest.TestCase):
    def setUp(self):
        self.positives = [get_data_path(e) for e in [
            'stockholm_extensive',
            ]]

        self.negatives = [get_data_path(e) for e in [
            'stockholm_data_only',
            'empty',
            ]]

    def test_positives(self):
        for fp in self.positives:
            self.assertEqual(_stockholm_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negatives:
            self.assertEqual(_stockholm_sniffer(fp), (False, {}))


class TestStockholmReader(unittest.TestCase):
    def test_stockholm_valid_extensive(self):
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

    def test_nonexistent_data_errors(self):
        fp = get_data_path('stockholm_invalid_nonexistent_gr')
        with assertRaisesRegex(self, StockholmFormatError,
                               'Markup line references.*RL1355.'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)
        fp = get_data_path('stockholm_invalid_nonexistent_gs')
        with assertRaisesRegex(self, StockholmFormatError,
                               'Markup line references.*AC14.'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_duplicate_data_error(self):
        fp = get_data_path('stockholm_duplicate_data')
        with assertRaisesRegex(self, StockholmFormatError,
                               'Found multiple data.*ASR132'):
            _stockholm_to_tabular_msa(fp, constructor=RNA)

    def test_no_data_error(self):
        fp = get_data_path('empty')
        with assertRaisesRegex(self, StockholmFormatError, 'No data present'):
            _stockholm_to_tabular_msa(fp)
        fp = get_data_path('stockholm_no_data')
        with assertRaisesRegex(self, StockholmFormatError, 'No data present'):
            _stockholm_to_tabular_msa(fp)

    def test_duplicate_label_errors(self):
        fp = get_data_path('stockholm_duplicate_gr')
        with assertRaisesRegex(self, StockholmFormatError,
                               'Found duplicate GR.*OS.*LFDR3'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)
        fp = get_data_path('stockholm_duplicate_gc')
        with assertRaisesRegex(self, StockholmFormatError,
                               'Found duplicate GC.*SS_cons'):
            _stockholm_to_tabular_msa(fp, constructor=DNA)

if __name__ == '__main__':
    unittest.main()
