# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from io import StringIO

import pandas as pd
import numpy as np

from skbio.util import get_data_path, assert_data_frame_almost_equal
from skbio.util._testing import _data_frame_to_default_int_type
from skbio.io.format.taxdump import _taxdump_to_data_frame


class TestTaxdumpReader(unittest.TestCase):
    def test_nodes_default(self):
        # subset of a real NCBI taxonomy nodes.dmp file
        fp = get_data_path('taxdump_nodes.dmp')
        obs = _taxdump_to_data_frame(fp, scheme='nodes')
        exp = pd.DataFrame([
            [1, 1, 'no rank', np.nan, 8, False, 1, False,
             0, False, False, False, np.nan],
            [2, 131567, 'superkingdom', np.nan, 0, False, 11, False,
             0, False, False, False, np.nan],
            [543, 91347, 'family', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [548, 570, 'species', 'KA', 0, True, 11, True,
             0, True, True, False, np.nan],
            [561, 543, 'genus', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [562, 561, 'species', 'EC', 0, True, 11, True,
             0, True, True, False, np.nan],
            [570, 543, 'genus', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [620, 543, 'genus', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [622, 620, 'species', 'SD', 0, True, 11, True,
             0, True, True, False, np.nan],
            [766, 28211, 'order', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [1224, 2, 'phylum', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [1236, 1224, 'class', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [28211, 1224, 'class', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [91347, 1236, 'order', np.nan, 0, True, 11, True,
             0, True, False, False, np.nan],
            [118884, 1236, 'no rank', np.nan, 0, True, 11, True,
             0, True, True, False, np.nan],
            [126792, 36549, 'species', 'PP', 0, True, 11, True,
             0, True, True, False, np.nan],
            [131567, 1, 'no rank', np.nan, 8, True, 1, True,
             0, True, True, False, np.nan],
            [585056, 562, 'no rank', np.nan, 0, True, 11, True,
             0, True, True, False, np.nan],
            [1038927, 562, 'no rank', np.nan, 0, True, 11, True,
             0, True, True, False, np.nan],
            [2580236, 488338, 'species', 'SE', 7, True, 11, True,
             0, True, False, False, np.nan]], columns=[
                 'tax_id', 'parent_tax_id', 'rank', 'embl_code',
                 'division_id', 'inherited_div_flag', 'genetic_code_id',
                 'inherited_GC_flag', 'mitochondrial_genetic_code_id',
                 'inherited_MGC_flag', 'GenBank_hidden_flag',
                 'hidden_subtree_root_flag', 'comments']).set_index('tax_id')
        exp['comments'] = exp['comments'].astype('O')
        _data_frame_to_default_int_type(exp)
        assert_data_frame_almost_equal(obs, exp)

    def test_names_default(self):
        # subset of a real NCBI taxonomy names.dmp file
        fp = get_data_path('taxdump_names.dmp')
        obs = _taxdump_to_data_frame(fp, scheme='names')
        exp = pd.DataFrame([
            [1, 'root', np.nan, 'scientific name'],
            [2, 'Bacteria', 'Bacteria <bacteria>', 'scientific name'],
            [2, 'eubacteria', np.nan, 'genbank common name'],
            [543, 'Enterobacteriaceae', np.nan, 'scientific name'],
            [548, 'Klebsiella aerogenes', np.nan, 'scientific name'],
            [561, 'Escherichia', np.nan, 'scientific name'],
            [562, '"Bacillus coli" Migula 1895', np.nan, 'authority'],
            [562, 'Escherichia coli', np.nan, 'scientific name'],
            [562, 'Escherichia/Shigella coli', np.nan, 'equivalent name'],
            [570, 'Donovania', np.nan, 'synonym'],
            [570, 'Klebsiella', np.nan, 'scientific name'],
            [620, 'Shigella', np.nan, 'scientific name'],
            [622, 'Shigella dysenteriae', np.nan, 'scientific name'],
            [766, 'Rickettsiales', np.nan, 'scientific name'],
            [1224, 'Proteobacteria', np.nan, 'scientific name'],
            [1236, 'Gammaproteobacteria', np.nan, 'scientific name'],
            [28211, 'Alphaproteobacteria', np.nan, 'scientific name'],
            [91347, 'Enterobacterales', np.nan, 'scientific name'],
            [118884, 'unclassified Gammaproteobacteria', np.nan,
             'scientific name'],
            [126792, 'Plasmid pPY113', np.nan, 'scientific name'],
            [131567, 'cellular organisms', np.nan, 'scientific name'],
            [585056, 'Escherichia coli UMN026', np.nan, 'scientific name'],
            [1038927, 'Escherichia coli O104:H4', np.nan, 'scientific name'],
            [2580236, 'synthetic Escherichia coli Syn61', np.nan,
             'scientific name']],
             columns=['tax_id', 'name_txt', 'unique_name',
                      'name_class']).set_index('tax_id')
        assert_data_frame_almost_equal(obs, exp)

    def test_nodes_slim(self):
        fp = get_data_path('taxdump_nodes.dmp')
        obs = _taxdump_to_data_frame(fp, scheme='nodes_slim')
        exp = pd.DataFrame([
            [1,       1,      'no rank'],
            [2,       131567, 'superkingdom'],
            [543,     91347,  'family'],
            [548,     570,    'species'],
            [561,     543,    'genus'],
            [562,     561,    'species'],
            [570,     543,    'genus'],
            [620,     543,    'genus'],
            [622,     620,    'species'],
            [766,     28211,  'order'],
            [1224,    2,      'phylum'],
            [1236,    1224,   'class'],
            [28211,   1224,   'class'],
            [91347,   1236,   'order'],
            [118884,  1236,   'no rank'],
            [126792,  36549,  'species'],
            [131567,  1,      'no rank'],
            [585056,  562,    'no rank'],
            [1038927, 562,    'no rank'],
            [2580236, 488338, 'species']],
            columns=['tax_id', 'parent_tax_id', 'rank']).set_index('tax_id')
        _data_frame_to_default_int_type(exp)
        assert_data_frame_almost_equal(obs, exp)

    def test_custom_scheme(self):
        fs = StringIO('\n'.join(map('\t|\t'.join, [
            ('a', 'a'),
            ('b', 'a'),
            ('c', 'a')
        ])))
        obs = _taxdump_to_data_frame(fs, scheme={'self': str, 'parent': str})
        exp = pd.DataFrame([
            ['a', 'a'],
            ['b', 'a'],
            ['c', 'a']],
            columns=['self', 'parent']).set_index('self')
        assert_data_frame_almost_equal(obs, exp)

    def test_invalid_scheme(self):
        fp = get_data_path('taxdump_names.dmp')
        with self.assertRaises(ValueError) as ctx:
            _taxdump_to_data_frame(fp, scheme='hello')
        self.assertEqual(str(ctx.exception),
                         'Invalid taxdump column scheme: "hello".')

    def test_invalid_id(self):
        fs = StringIO('\n'.join(map('\t|\t'.join, [
            ('1', '2', 'family'),
            ('3', '4', 'genus'),
            ('x', '6', 'species'),  # 'x' is not a number
        ])))
        with self.assertRaises(ValueError) as ctx:
            _taxdump_to_data_frame(fs, scheme='nodes_slim')
        self.assertEqual(str(ctx.exception),
                         'Invalid taxdump file format.')


if __name__ == '__main__':
    unittest.main()
