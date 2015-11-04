import unittest
import subprocess
import goenrich

class TestGoenrich(unittest.TestCase):

    def test_analysis_and_export(self):
        O = goenrich.obo.ontology('db/go-basic.obo')
        gene2go = goenrich.read.gene2go('db/gene2go.gz')
        values = {k: set(v) for k,v in gene2go.groupby('GO_ID')['GeneID']}
        background_attribute = 'gene2go'
        goenrich.enrich.propagate(O, values, background_attribute)
        query = gene2go['GeneID'].unique()[:20]
        try:
            import pygraphviz
            goenrich.enrich.analyze(O, query, background_attribute, gvfile='test.dot')
            subprocess.check_call(['dot', '-Tpng', 'test.dot', '-o', 'test.png'])
            subprocess.check_call(['rm', 'test.dot', 'test.png'])
        except ImportError:
            goenrich.enrich.analyze(O, query, background_attribute)
            print('pygraphviz could not be imported')

if __name__ == '__main__':
    unittest.main()
