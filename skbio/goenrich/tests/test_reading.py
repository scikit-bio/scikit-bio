import unittest

import goenrich

class TestRead(unittest.TestCase):
    def test_ontology(self):
        G = goenrich.obo.ontology('db/go-basic.obo')

    def test_goa(self):
        background = goenrich.read.goa('db/gene_association.goa_ref_human.gz')

    def test_gene2go(self):
        background = goenrich.read.gene2go('db/gene2go.gz')

if __name__ == '__main__':
    unittest.main()
