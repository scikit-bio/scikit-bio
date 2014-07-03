r"""
E-utilities automation module (:mod:`skbio.db.ncbi`)
====================================================

.. currentmodule:: skbio.db.ncbi

EUtils [1]_. is a web service offered by the NCBI to access the sequence,
literature and other databases by a special format of URLs. This module offers
an interface to construct the URLs and retrieve the results in text format.

Classes
-------

.. autosummary::
   :toctree: generated/

   EFetch
   ELink
   ESearch
   EUtils

References
----------
.. [1] http://www.ncbi.nih.gov/entrez/eutils

"""
from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from future.utils.six import StringIO

from skbio.db.ncbi import (EUtils, ESearch, EFetch, ELink,
                           _esearch_result_parser, _elink_result_parser,
                           _get_primary_ids, _ids_to_taxon_ids,
                           _taxon_lineage_extractor, _taxon_ids_to_lineages,
                           _taxon_ids_to_names,
                           _taxon_ids_to_names_and_lineages,
                           _get_unique_lineages, _get_unique_taxa,
                           _parse_taxonomy_using_elementtree_xml_parse)

# py3k compatibility
try:
    from string import strip
except ImportError:
    strip = str.strip

class EUtilsTests(TestCase):
    def test_simple_get(self):
        g = EUtils(db='protein', rettype='gp')
        result = g['NP_003320'].read()
        self.assertTrue(result.startswith('LOCUS'))
        self.assertTrue('NP_003320' in result)

    def test_get_slice(self):
        g = EUtils(db='protein', rettype='gp', retmax=1)
        result = g['NP_003320':'NP_003322'].read()
        lines = result.splitlines()
        is_locus = lambda x: x.startswith('LOCUS')
        loci = filter(is_locus, lines)
        self.assertEqual(len(loci), 3)

        # EUtils access of a slice should work, while limiting
        # the esearch term length
        g = EUtils(db='protein', rettype='gp', retmax=1, url_limit=2)
        result = g['NP_003320':'NP_003322'].read()
        lines = result.splitlines()
        is_locus = lambda x: x.startswith('LOCUS')
        loci = filter(is_locus, lines)
        self.assertEqual(len(loci), 3)

    def test_get_list(self):
        g = EUtils(db='protein', rettype='gp')
        result = g['NP_003320', 'NP_003321', 'NP_003322'].read()
        lines = result.splitlines()
        is_locus = lambda x: x.startswith('LOCUS')
        loci = filter(is_locus, lines)
        self.assertEqual(len(loci), 3)

        # EUtils access of a slice should work, while limiting
        # the esearch term length
        g = EUtils(db='protein', rettype='gp', url_limit=2)
        result = g['NP_003320', 'NP_003321', 'NP_003322'].read()
        lines = result.splitlines()
        is_locus = lambda x: x.startswith('LOCUS')
        loci = filter(is_locus, lines)
        self.assertEqual(len(loci), 3)

    def test_get_from_taxonomy_db(self):
        # note: this is more fragile than the nucleotide databases
        g = EUtils(db='taxonomy', rettype='xml', retmode='xml')
        ids = '9606[taxid] OR 28901[taxid]'
        fh = StringIO()
        fh.write(g[ids].read())
        fh.seek(0)
        data = _parse_taxonomy_using_elementtree_xml_parse(fh)
        result = sorted([item['ScientificName'] for item in data])
        self.assertEqual(result, ['Homo sapiens', 'Salmonella enterica'])

    def test_query(self):
        g = EUtils(db='protein', rettype='gi', retmax=100)
        result = g['homo[organism] AND erf1[ti]'].read().splitlines()
        self.assertTrue('5499721' in result)  # gi of human eRF1
        # note: successfully retrieved 841,821 ids on a query for 'rrna',
        # although it took about 40 min so not in the tests. RK 1/3/07.

    def test_query_retmax(self):
        g = EUtils(db='protein', rettype='gi', retmax=3, DEBUG=False)
        result = g['homo[organism] AND myh7'].read().splitlines()
        self.assertTrue(len(result) > 1)
        self.assertTrue('83304912' in result)  # gi of human myh7

    def test_query_max_recs(self):
        g = EUtils(db='protein', rettype='gi', max_recs=5, DEBUG=False,
                   retmax=100)
        result = g['homo[organism] AND myh7'].read().splitlines()
        self.assertEqual(len(result), 5)

    def test_query_max_recs_gt_retmax(self):
        g = EUtils(db='protein', rettype='gi', max_recs=5, DEBUG=False,
                   retmax=3)
        result = g['homo[organism] AND myh7'].read().splitlines()
        self.assertEqual(len(result), 5)


class ESearchTests(TestCase):
    def test_simple_search(self):
        s = ESearch(db='protein', rettype='gi', retmax=1000,
                    term='homo[organism] AND myh7')
        result = s.read()
        parsed = _esearch_result_parser(result)

        # gi of human cardiac beta myh7
        self.assertTrue('83304912' in parsed.IdList)


class ELinkTests(TestCase):
    def test_simple_elink(self):
        l = ELink(db='taxonomy', dbfrom='protein', id='83304912')
        result = l.read()
        parsed = _elink_result_parser(result)
        self.assertEqual(parsed, ['9606'])  # human sequence

    def test_multiple_elink(self):
        l = ELink(db='taxonomy', dbfrom='protein',
                  id='83304912 115496169 119586556 111309484')
        result = l.read()
        parsed = _elink_result_parser(result)
        self.assertEqual(sorted(parsed), ['10090', '9606'])


class EFetchTests(TestCase):
    def test_simple_efetch(self):
        f = EFetch(db='protein', rettype='fasta', retmode='text',
                   id='111309484')
        result = f.read().splitlines()
        self.assertTrue(result[0].startswith(b'>'))
        self.assertTrue(result[1].startswith(b'madaemaafg'.upper()))


class NcbiTests(TestCase):
    def setUp(self):
        self.maxDiff = None
        self.mouse_taxonomy = ['cellular organisms', 'Eukaryota',
                               'Opisthokonta', 'Metazoa', 'Eumetazoa',
                               'Bilateria', 'Deuterostomia', 'Chordata',
                               'Craniata', 'Vertebrata', 'Gnathostomata',
                               'Teleostomi', 'Euteleostomi', 'Sarcopterygii',
                               'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota',
                               'Mammalia', 'Theria', 'Eutheria',
                               'Boreoeutheria', 'Euarchontoglires', 'Glires',
                               'Rodentia', 'Sciurognathi', 'Muroidea',
                               'Muridae', 'Murinae', 'Mus', 'Mus']

        self.human_taxonomy = ['cellular organisms', 'Eukaryota',
                               'Opisthokonta', 'Metazoa', 'Eumetazoa',
                               'Bilateria', 'Deuterostomia', 'Chordata',
                               'Craniata', 'Vertebrata', 'Gnathostomata',
                               'Teleostomi', 'Euteleostomi', 'Sarcopterygii',
                               'Dipnotetrapodomorpha', 'Tetrapoda', 'Amniota',
                               'Mammalia', 'Theria', 'Eutheria',
                               'Boreoeutheria', 'Euarchontoglires', 'Primates',
                               'Haplorrhini', 'Simiiformes', 'Catarrhini',
                               'Hominoidea', 'Hominidae', 'Homininae', 'Homo']

    def test_get_primary_ids(self):
        res = _get_primary_ids('homo[orgn] AND myh7[ti]', retmax=5, max_recs=7)
        self.assertEqual(len(res), 7)
        res = _get_primary_ids('homo[orgn] AND myh7[ti]', retmax=5, max_recs=2)
        self.assertEqual(len(res), 2)
        res = _get_primary_ids('homo[orgn] AND myh7[ti]', retmax=100)

        # this may have to be updated as versions change, it previously was
        # 115496168. See this site and the revision history
        # http://www.ncbi.nlm.nih.gov/nuccore/NM_000257.3
        self.assertTrue('657940851' in res)

    def test_ids_to_taxon_ids(self):
        ids = ['83304912', '115496169', '119586556', '111309484']
        result = _ids_to_taxon_ids(ids, db='protein')
        self.assertEqual(sorted(result), ['10090', '9606'])

    def test_taxon_lineage_extractor(self):
        lines = """ignore
        <Lineage>xxx;yyy</Lineage>
        ignore

        <Lineage>aaa;bbb</Lineage>
        """
        self.assertEqual(list(_taxon_lineage_extractor(lines.splitlines())),
                         [['xxx', 'yyy'], ['aaa', 'bbb']])

    def test_parse_taxonomy_using_elementtree_xml_parse(self):
        g = EUtils(db='taxonomy', rettype='xml', retmode='xml')
        ids = '28901[taxid]'
        fh = StringIO()
        fh.write(g[ids].read())
        fh.seek(0)
        data = _parse_taxonomy_using_elementtree_xml_parse(fh)[0]
        obs = (data['Lineage'], data['TaxId'], data['ScientificName'],
               data['Rank'])
        exp = (
            'cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria'
            '; Enterobacteriales; Enterobacteriaceae; Salmonella',
            '28901',
            'Salmonella enterica',
            'species')
        self.assertEqual(obs, exp)

    def test_taxon_ids_to_lineages(self):
        taxon_ids = ['10090', '9606']
        result = [self.mouse_taxonomy, self.human_taxonomy]

        if hasattr(self, 'assertItemsEqual'):
            self.assertItemsEqual(list(_taxon_ids_to_lineages(taxon_ids)),
                                  result)
        else:
            self.assertCountEqual(list(_taxon_ids_to_lineages(taxon_ids)),
                                  result)

    def test_taxon_ids_to_names(self):
        taxon_ids = ['10090', '9606']
        result = set(['Mus musculus', 'Homo sapiens'])
        self.assertEqual(set(_taxon_ids_to_names(taxon_ids)), result)

    def test_taxon_ids_to_names_and_lineages(self):
        taxon_ids = ['10090', '9606']
        exp = [('10090', 'Mus musculus', '; '.join(self.mouse_taxonomy)),
               ('9606', 'Homo sapiens', '; '.join(self.human_taxonomy))]
        obs = list(_taxon_ids_to_names_and_lineages(taxon_ids))

        if hasattr(self, 'assertItemsEqual'):
            self.assertItemsEqual(obs, exp)
        else:
            self.assertCountEqual(obs, exp)

    def test_get_unique_lineages(self):
        result = _get_unique_lineages('angiotensin[ti] AND rodents[orgn]')

        self.assertTrue(tuple(self.mouse_taxonomy) in result)
        self.assertTrue(len(result) > 2)

    def test_get_unique_taxa(self):
        result = _get_unique_taxa('angiotensin[ti] AND primate[orgn]')
        self.assertTrue('Homo sapiens' in result)
        self.assertTrue(len(result) > 2)


if __name__ == '__main__':
    main()
