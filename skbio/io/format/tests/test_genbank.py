from __future__ import absolute_import, division, print_function
from future.builtins import map, range, zip
import six

import io
import numpy as np
import numpy.testing as nptest
from unittest import TestCase, main
from functools import partial
from datetime import datetime

from skbio import Protein, DNA, RNA
from skbio.util import get_data_path
from skbio.io import GenbankFormatError
from skbio.io.format.genbank import (
    _genbank_to_generator, _genbank_to_biological_sequence,
    _parse_genbanks, _parse_single_genbank,
    _parse_locus, _parse_reference,
    _parse_source,  _parse_origin,
    _parse_features, _parse_single_feature, _parse_loc_str,
    _parse_section_default,
    yield_section)


class SnifferTests(TestCase):
    pass


class ReaderTests(TestCase):
    def setUp(self):
        self.valid = []

    def test_parse_reference(self):
        lines = '''REFERENCE   1  (bases 1 to 154478)
  AUTHORS   Sato,S., Nakamura,Y., Kaneko,T., and Tabata,S.
  TITLE     Complete structure of the chloroplast genome of
            Arabidopsis thaliana
  JOURNAL   DNA Res. 6 (5), 283-290 (1999)
   PUBMED   10574454'''.split('\n')

        exp = {'AUTHORS': 'Sato,S., Nakamura,Y., Kaneko,T., and Tabata,S.',
               'JOURNAL': 'DNA Res. 6 (5), 283-290 (1999)',
               'PUBMED': '10574454',
               'REFERENCE': '1  (bases 1 to 154478)',
               'TITLE': ('Complete structure of the chloroplast genome of'
                         ' Arabidopsis thaliana')}
        self.assertEqual(_parse_reference(lines), exp)

    def test_parse_locus(self):
        lines = [
            ['LOCUS       NC_005816               9609 bp'
             '    DNA     circular CON 07-FEB-2015'],
            ['LOCUS       SCU49845     5028 bp'
             '    DNA             PLN       21-JUN-1999'],
            ['LOCUS       NP_001832                360 aa'
             '            linear   PRI 18-DEC-2001']]
        expects = [
            {'division': 'CON', 'mol_type': DNA, 'shape': 'circular',
             'locus_name': 'NC_005816', 'date': datetime(2015, 2, 7, 0, 0),
             'unit': 'bp', 'size': 9609},
            {'division': 'PLN', 'mol_type': DNA, 'shape': None,
             'locus_name': 'SCU49845', 'date': datetime(1999, 6, 21, 0, 0),
             'unit': 'bp', 'size': 5028},
            {'division': 'PRI', 'mol_type': Protein, 'shape': 'linear',
             'locus_name': 'NP_001832', 'date': datetime(2001, 12, 18, 0, 0),
             'unit': 'aa', 'size': 360}]

        for line, exp in zip(lines, expects):
            self.assertEqual(_parse_locus(line), exp)

    def test_parse_locus_invalid(self):
        lines = [
            # missing unit
            'LOCUS       NC_005816               9609 '
            '    DNA     circular CON 07-FEB-2015',
            # missing division
            'LOCUS       SCU49845     5028 bp'
            '    DNA                    21-JUN-1999',
            # wrong date format
            'LOCUS       NP_001832                360 aa'
            '            linear   PRI 2001-12-18']
        for line in lines:
            with self.assertRaises(GenbankFormatError):
                print(_parse_locus(line))

    def test_parse_single_feature(self):
        lines='''     CDS             join(1713..1891,2322..2438,2538..2633,2801..2843,
                     2918..3073,3167..3247,3874..3972,4082..4309)
                     /gene="CCT"
                     /EC_number="2.7.7.15"
                     /codon_start=1
                     /product="CTP:phosphocholine cytidylyltransferase"
                     /protein_id="AAD45922.1"
                     /db_xref="GI:5640001"
                     /translation="MSNVIGDRTEDGLSTAAAASGSTAVQSSPPTDRPVRVYADGIYD
                     LFHFGHARSLEQAKLAFPNNTYLLVGCCNDETTHKYKGRTVMTAEERYESLRHCKWVD
                     EVIPDAPWVVNQEFLDKHQIDYVAHDSLPYADSSGAGKDVYEFVKKVGRFKETQRTEG
                     ISTSDIIMRIVKDYNQYVMRNLDRGYSREDLGVSFVKEKRLRVNMRLKKLQERVKEQQ
                     ERVGEKIQTVKMLRNEWVENADRWVAGFLEIFEEGCHKMGTAIVDSIQERLMRQKSAE
                     RLENGQDDDTDDQFYEEYFDHDMGSDDDEDEKFYDEEEVKEEETEKTVMTDAKDNK"'''.split('\n')
        print(_parse_single_feature(lines, 300, 0))

    def test_parse_features(self):
        lines = '''FEATURES             Location/Qualifiers
     source          1..5485
                     /organism="Arabidopsis thaliana"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:3702"
                     /ecotype="Col-0"
     gene            1..4637
                     /gene="CCT"
     regulatory      1..1602
                     /regulatory_class="promoter"
                     /gene="CCT"
     regulatory      1554..1560
                     /regulatory_class="TATA_box"
                     /gene="CCT"
     mRNA            join(1603..1891,2322..2438,2538..2633,2801..2843,
                     2918..3073,3167..3247,3874..3972,4082..4637)
                     /gene="CCT"
                     /product="CTP:phosphocholine cytidylyltransferase"
     5'UTR           1603..1712
                     /gene="CCT"
     CDS             join(1713..1891,2322..2438,2538..2633,2801..2843,
                     2918..3073,3167..3247,3874..3972,4082..4309)
                     /gene="CCT"
                     /EC_number="2.7.7.15"
                     /codon_start=1
                     /product="CTP:phosphocholine cytidylyltransferase"
                     /protein_id="AAD45922.1"
                     /db_xref="GI:5640001"
                     /translation="MSNVIGDRTEDGLSTAAAASGSTAVQSSPPTDRPVRVYADGIYD
                     LFHFGHARSLEQAKLAFPNNTYLLVGCCNDETTHKYKGRTVMTAEERYESLRHCKWVD
                     EVIPDAPWVVNQEFLDKHQIDYVAHDSLPYADSSGAGKDVYEFVKKVGRFKETQRTEG
                     ISTSDIIMRIVKDYNQYVMRNLDRGYSREDLGVSFVKEKRLRVNMRLKKLQERVKEQQ
                     ERVGEKIQTVKMLRNEWVENADRWVAGFLEIFEEGCHKMGTAIVDSIQERLMRQKSAE
                     RLENGQDDDTDDQFYEEYFDHDMGSDDDEDEKFYDEEEVKEEETEKTVMTDAKDNK"
     3'UTR           4310..4637
                     /gene="CCT"'''.split('\n')
        pprint(_parse_features(lines, 5485))


    def test_parse_loc_str(self):
        length = 12

        examples = [
            '9',  # a single base in the presented sequence
            '3..8',
            '<3..8',
            '1..>8',
            'complement(3..8)',
            'complement(join(3..5,7..9))',
            'join(3..5,7..9)']

        expects = [
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': True, 'rc_': False},
             np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)),
            ({'right_partial_': True, 'left_partial_': False, 'rc_': False},
             np.array([1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': True},
             np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': True},
             np.array([0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0], dtype=bool)),
            ({'right_partial_': False, 'left_partial_': False, 'rc_': False},
             np.array([0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0], dtype=bool))]
        for example, expect in zip(examples, expects):
            parsed = _parse_loc_str(example, length)
            self.assertDictEqual(parsed[0], expect[0])
            nptest.assert_equal(parsed[1], expect[1])

    def test_parse_genbanks(self):
        fp = get_data_path('genbank_misc_gene_features')
        with open(fp) as fh:
            a = list(_parse_genbanks(fh))[0]
            pprint(a)
            import ipdb; ipdb.set_trace()

    def test_genbank_to_generator(self):
        fp = get_data_path('genbank_misc_gene_features')
        a = _genbank_to_generator(fp)
        x = list(a)[0]
        import ipdb; ipdb.set_trace()

    def test_genbank_to_biological_sequence(self):
        fp = get_data_path('genbank_multi_records')
        print(fp)
        a = _genbank_to_biological_sequence(fp, rec_num=2)
        import ipdb; ipdb.set_trace()


class WriterTests(TestCase):
    pass


if __name__ == '__main__':
    main()
