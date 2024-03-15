r"""Taxdump format (:mod:`skbio.io.format.taxdump`)
===============================================

.. currentmodule:: skbio.io.format.taxdump

The NCBI Taxonomy database dump (``taxdump``) format stores information of
organism names, classifications and other properties. It is a tabular format
with a delimiter: ``<tab><pipe><tab>`` between columns, and a line end
``<tab><pipe>`` after all columns. The file name usually ends with .dmp.

Format Support
--------------
**Has Sniffer: No**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`pandas.DataFrame`                                        |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------

The NCBI taxonomy database [1]_ [2]_ hosts organism names and classifications.
It has a web portal [3]_ and an FTP download server [4]_. It is also accessible
using E-utilities [5]_. The database is being updated daily, and an archive is
generated every month. The data release has the file name ``taxdump``. It
consists of multiple .dmp files. These files serve different purposes, but they
follow a common format pattern:

- It is a tabular format.
- Column delimiter is ``<tab><pipe><tab>``.
- Line end is ``<tab><pipe>``.
- The first column is a numeric identifier, which usually represent taxa (i.e.,
  "TaxID"), but can also be genetic codes, citations or other entries.

The two most important files of the data release are ``nodes.dmp`` and
``names.dmp``. They store the hierarchical structure of the classification
system (i.e., taxonomy) and the names of organisms, respectively. They can be
used to construct the taxonomy tree of organisms.

The definition of columns of each .dmp file type are taken from [6]_ and [7]_.

``nodes.dmp``
^^^^^^^^^^^^^
+----------------+-------------------------------------+
|Name            |Description                          |
+================+=====================================+
|tax_id          |node id in GenBank taxonomy database |
+----------------+-------------------------------------+
|parent tax_id   |parent node id in GenBank taxonomy   |
|                |database                             |
+----------------+-------------------------------------+
|rank            |rank of this node (superkingdom,     |
|                |kingdom, ...)                        |
+----------------+-------------------------------------+
|embl code       |locus-name prefix; not unique        |
+----------------+-------------------------------------+
|division id     |see division.dmp file                |
+----------------+-------------------------------------+
|inherited div   |1 if node inherits division from     |
|flag (1 or 0)   |parent                               |
+----------------+-------------------------------------+
|genetic code id |see gencode.dmp file                 |
+----------------+-------------------------------------+
|inherited GC    |1 if node inherits genetic code from |
|flag (1 or 0)   |parent                               |
+----------------+-------------------------------------+
|mitochondrial   |see gencode.dmp file                 |
|genetic code id |                                     |
+----------------+-------------------------------------+
|inherited MGC   |1 if node inherits mitochondrial     |
|flag (1 or 0)   |gencode from parent                  |
+----------------+-------------------------------------+
|GenBank hidden  |1 if name is suppressed in GenBank   |
|flag (1 or 0)   |entry lineage                        |
+----------------+-------------------------------------+
|hidden subtree  |1 if this subtree has no sequence    |
|root flag       |data yet                             |
|(1 or 0)        |                                     |
+----------------+-------------------------------------+
|comments        |free-text comments and citations     |
+----------------+-------------------------------------+

Since 2018, NCBI releases "new taxonomy files" [8]_ (``new_taxdump``). The new
``nodes.dmp`` format is compatible with the classical format, plus five extra
columns after all aforementioned columns.

+----------------+-------------------------------------+
|Name            |Description                          |
+================+=====================================+
|plastid genetic |see gencode.dmp file                 |
|code id         |                                     |
+----------------+-------------------------------------+
|inherited PGC   |1 if node inherits plastid gencode   |
|flag (1 or 0)   |from parent                          |
+----------------+-------------------------------------+
|specified       |1 if species in the node's lineage   |
|species         |has formal name                      |
+----------------+-------------------------------------+
|hydrogenosome   |see gencode.dmp file                 |
|genetic code id |                                     |
+----------------+-------------------------------------+
|inherited HGC   |1 if node inherits hydrogenosome     |
|flag (1 or 0)   |gencode from parent                  |
+----------------+-------------------------------------+

``names.dmp``
^^^^^^^^^^^^^
+----------------+-------------------------------------+
|Name            |Description                          |
+================+=====================================+
|tax_id          |the id of node associated with this  |
|                |name                                 |
+----------------+-------------------------------------+
|name_txt        |name itself                          |
+----------------+-------------------------------------+
|unique name     |the unique variant of this name if   |
|                |name not unique                      |
+----------------+-------------------------------------+
|name class      |(synonym, common name, ...)          |
+----------------+-------------------------------------+

``division.dmp``
^^^^^^^^^^^^^^^^
+----------------+-------------------------------------+
|Name            |Description                          |
+================+=====================================+
|division id     |taxonomy database division id        |
+----------------+-------------------------------------+
|division cde    |GenBank division code (three         |
|                |characters)                          |
+----------------+-------------------------------------+
|division name   |e.g. BCT, PLN, VRT, MAM, PRI...      |
+----------------+-------------------------------------+
|comments        |                                     |
+----------------+-------------------------------------+

``gencode.dmp``
^^^^^^^^^^^^^^^
+----------------+-------------------------------------+
|Name            |Description                          |
+================+=====================================+
|genetic code id |GenBank genetic code id              |
+----------------+-------------------------------------+
|abbreviation    |genetic code name abbreviation       |
+----------------+-------------------------------------+
|name            |genetic code name                    |
+----------------+-------------------------------------+
|cde             |translation table for this genetic   |
|                |code                                 |
+----------------+-------------------------------------+
|starts          |start codons for this genetic code   |
+----------------+-------------------------------------+

Other types of .dmp files are currently not supported by scikit-bio. However,
the user may customize column definitions in using this utility. See below for
details.

Format Parameters
-----------------
The following format parameters are available in ``taxdump`` format:

- ``scheme``: The column definition scheme name of the input .dmp file.
  Available options are listed below. Alternatively, one can provide a custom
  scheme as defined in a name-to-data type dictionary.

  1. ``nodes``: The classical ``nodes.dmp`` scheme. It is also compatible with
     new ``nodes.dmp`` format, in which case only the columns defined by the
     classical format will be read.

  2. ``nodes_new``: The new ``nodes.dmp`` scheme.

  3. ``nodes_slim``: Only the first three columns: tax_id, parent_tax_id and
     rank, which are the minimum required information for constructing the
     taxonomy tree. It can be applied to both classical and new ``nodes.dmp``
     files. It can also handle custom files which only contains these three
     columns.

  4. ``names``: The ``names.dmp`` scheme.

  5. ``division``: The ``division.dmp`` scheme.

  6. ``gencode``: The ``gencode.dmp`` scheme.

.. note:: scikit-bio will read columns from leftmost till the number of columns
   defined in the scheme. Extra columns will be cropped.

Examples
--------
>>> from io import StringIO
>>> import skbio.io
>>> import pandas as pd
>>> fs = '\n'.join([
...     '1\t|\t1\t|\tno rank\t|',
...     '2\t|\t131567\t|\tsuperkingdom\t|',
...     '6\t|\t335928\t|\tgenus\t|'
... ])
>>> fh = StringIO(fs)

Read the file into a ``pd.DataFrame`` and specify that the "nodes_slim" scheme
should be used:

>>> df = skbio.io.read(fh, format="taxdump", into=pd.DataFrame,
...                    scheme="nodes_slim")
>>> df # doctest: +NORMALIZE_WHITESPACE
        parent_tax_id          rank
tax_id
1                   1       no rank
2              131567  superkingdom
6              335928         genus

References
----------
.. [1] Federhen, S. (2012). The NCBI taxonomy database. Nucleic acids
   research, 40(D1), D136-D143.
.. [2] Schoch, C. L., Ciufo, S., Domrachev, M., Hotton, C. L., Kannan, S.,
   Khovanskaya, R., ... & Karsch-Mizrachi, I. (2020). NCBI Taxonomy: a
   comprehensive update on curation, resources and tools. Database, 2020.
.. [3] https://www.ncbi.nlm.nih.gov/taxonomy
.. [4] https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
.. [5] Kans, J. (2022). Entrez direct: E-utilities on the UNIX command line.
   In Entrez Programming Utilities Help [Internet]. National Center for
   Biotechnology Information (US).
.. [6] https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt
.. [7] https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/taxdump_readme.txt
.. [8] https://ncbiinsights.ncbi.nlm.nih.gov/2018/02/22/new-taxonomy-files-
       available-with-lineage-type-and-host-information/


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from skbio.io import create_format


taxdump = create_format("taxdump")

_taxdump_column_schemes = {
    "nodes_slim": {"tax_id": int, "parent_tax_id": int, "rank": str},
    "nodes": {
        "tax_id": int,
        "parent_tax_id": int,
        "rank": str,
        "embl_code": str,
        "division_id": int,
        "inherited_div_flag": bool,
        "genetic_code_id": int,
        "inherited_GC_flag": bool,
        "mitochondrial_genetic_code_id": int,
        "inherited_MGC_flag": bool,
        "GenBank_hidden_flag": bool,
        "hidden_subtree_root_flag": bool,
        "comments": str,
    },
    "names": {"tax_id": int, "name_txt": str, "unique_name": str, "name_class": str},
    "division": {
        "division_id": int,
        "division_cde": str,
        "division_name": str,
        "comments": str,
    },
    "gencode": {
        "genetic_code_id": int,
        "abbreviation": str,
        "name": str,
        "cde": str,
        "starts": str,
    },
}

_taxdump_column_schemes["nodes_new"] = dict(
    _taxdump_column_schemes["nodes"],
    **{
        "plastid_genetic_code_id": bool,
        "inherited_PGC_flag": bool,
        "specified_species": bool,
        "hydrogenosome_genetic_code_id": int,
        "inherited_HGC_flag": bool,
    },
)


@taxdump.reader(pd.DataFrame, monkey_patch=False)
def _taxdump_to_data_frame(fh, scheme):
    """Read a taxdump file into a data frame.

    Parameters
    ----------
    fh : file handle
        Input taxdump file
    scheme : str
        Name of column scheme

    Returns
    -------
    pd.DataFrame
        Parsed table

    """
    if isinstance(scheme, str):
        if scheme not in _taxdump_column_schemes:
            raise ValueError(f'Invalid taxdump column scheme: "{scheme}".')
        scheme = _taxdump_column_schemes[scheme]
    names = list(scheme.keys())
    try:
        return pd.read_csv(
            fh,
            sep="\t\\|(?:\t|$)",
            engine="python",
            index_col=0,
            names=names,
            dtype=scheme,
            usecols=range(len(names)),
        )
    except ValueError:
        raise ValueError("Invalid taxdump file format.")
