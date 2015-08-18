Files in this directory are the QIIME 1.9.1 "tiny test" files. These data were developed by @gregcaporaso, who gave permission to reproduce them in scikit-bio.

If you have a [QIIME 1.9.1 base installation](http://install.qiime.org), the raw input files in this directory can be obtained by running:

```bash
python -c "from qiime.test import write_test_data; write_test_data('.')"
biom convert -i biom --to-tsv -o otu-table.tsv
```

After converting to tsv, the following OTUs are removed because they are not present in the tree (they're not 16S sequences, so can't be aligned with PyNAST): ``None1``, ``None10``, ``None6``, and ``None2``. The ``not16S.1`` sample is also removed because, after removing those OTUs, it has a total count of 0. This boundary case is tested directly in the ``unifrac_*`` and ``faith_pd`` tests.

Then, in the python interpreter, we midpoint root the tree (since this is a QIIME 1.9.1 installation, this step is performed with scikit-bio 0.2.3):

```python
from skbio import TreeNode
t = TreeNode.read('./tree')
t = t.root_at_midpoint()
t.write('tree', format='newick')
```

The output files (alpha diversity values and beta diversity distance matrices) can then be obtained by running:

```bash
alpha_diversity.py -i biom -t tree -m PD_whole_tree -o pd.txt
beta_diversity.py -m weighted_unifrac,unweighted_unifrac,weighted_normalized_unifrac -i biom -t tree -o o
```
