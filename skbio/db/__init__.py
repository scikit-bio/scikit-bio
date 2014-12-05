r"""
Exploring online databases (:mod:`skbio.db`)
============================================

.. currentmodule:: skbio.db

EUtils is a web service offered by the NCBI to access the sequence, literature
and other databases by a special format of URLs. scikit-bio offers an interface
to construct the URLs and retrieve the results in text format.

From Pubmed
-----------

Retrieving PubMed records from NCBI by PubMed ID
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The process for getting PubMed records by PubMed ID (PMID) is very similar to
that for getting sequences. Basically, you just need to pass in the unique id
associated with the article. For example, if you want to get the reference to
the original PyCogent paper to see how far we've come since then, you can do
this:


>>> from skbio.db.ncbi import EFetch
>>> ef = EFetch(id='17708774', db='pubmed', rettype='brief')
>>> print(ef.read()) # doctest: +ELLIPSIS
<BLANKLINE>
1. Genome Biol. 2007;8(8):R171.
<BLANKLINE>
PyCogent: a toolkit for making sense from sequen...
<BLANKLINE>
<BLANKLINE>

If you want more information, there are other rettypes, e.g.

>>> ef = EFetch(id='17708774', db='pubmed', rettype='citation')
>>> print(ef.read()) # doctest: +ELLIPSIS
<BLANKLINE>
1. Genome Biol. 2007;8(8):R171.
<BLANKLINE>
PyCogent: a toolkit for making sense from sequence.
<BLANKLINE>
Knight R(1), Maxwell P, Birmingham A, Carnes J, Caporaso...
<BLANKLINE>
PMCID: PMC2375001
PMID: 17708774  [PubMed - indexed for MEDLINE]
<BLANKLINE>
<BLANKLINE>


Similarly, if you want something more machine-readable (but quite a lot less
human-readable), you can specify XML in the retmode:

>>> ef = EFetch(id='17708774', db='pubmed', rettype='citation', retmode='xml')
>>> cite = ef.read()
>>> for line in cite.splitlines()[:5]:
...     print(line) # doctest: +ELLIPSIS
...
<?xml version="1.0"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st...
<PubmedArticleSet>
<PubmedArticle>
    <MedlineCitation Owner="NLM" Status="MEDLINE">

Only a partial example is shown as there are quite a few lines. As with
sequences, you can retrieve multiple accessions at once.

Searching for records using EUtils
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Getting records by their primary identifiers is very nice if you actually have
the primary identifiers, but what if you don't? For example, what if you want
to do a search based on a keyword, or have a genbank accession number rather
than a gi, or want to get a range of records?

Fortunately, the more general EUtils class allows this kind of complex workflow
with relatively little intervention. For example, if you want to search for
articles that mention PyCogent:

>>> from skbio.db.ncbi import EUtils
>>> eu = EUtils(db='pubmed', rettype='brief')
>>> res = eu['PyCogent']
>>> print(res.read()) # doctest: +ELLIPSIS
<BLANKLINE>
1. J Appl Crystallogr. 2011 Apr 1;44(Pt 2):424-428...
2. RNA. 2008 Mar;14(3):410-6. doi: 10.1261/rna.881...
3. Genome Biol. 2007;8(8):R...

Perhaps you want only the ones with PyCogent in the title, in which case you
can use any qualifier that NCBI supports:

>>> res = eu['PyCogent[ti]']
>>> print(res.read()) # doctest: +ELLIPSIS
<BLANKLINE>
1. J Appl Crystallogr. 2011 Apr 1;44(Pt 2):424-4...
2. Genome Biol. 2007;8(8):R1...

The NCBI-supported list of field qualifiers, and lots of documentation
generally on how to do pubmed queries, is `here
<http://www.ncbi.nlm.nih.gov/bookshelf/br.fcgi?book=helppubmed&part=pubmedhelp>`_.

One especially useful feature is the ability to get a list of primary
identifiers matching a query. You do this by setting ``rettype='uilist'`` (not
idlist any more, so again you may need to update old code examples). For
example:

>>> eu = EUtils(db='pubmed', rettype='uilist')
>>> res = eu['PyCogent']
>>> print(res.read())
22479120
18230758
17708774
<BLANKLINE>

This is especially useful when you want to do a bunch of queries (whether for
journal articles, as shown here, or for sequences), combine the results, then
download the actual unique records only once. You could of course do this with
an incredibly complex single query, but good luck debugging that query...


For sequences
-------------

Fetching FASTA or Genpept sequences from NCBI using EFetch with GI's
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you already have a list of GI's (the numeric identifiers that are used by
GenBank internally as identifers), your job is very easy: you just need to use
``EFetch`` to retrieve the corresponding records. In general, this works for
any case where the identifiers you have are the primary keys, e.g. for PubMed
you use the PubMed ID (see example below).

Here is an example of getting the nucleotide record that corresponds to one
particular id, in this case id # 459567 (chosen arbitrarily). The record
arrives as a file-like object that can be read; in this case, we are looking at
each line and printing the first 40 characters.

>>> from skbio.db.ncbi import EFetch
>>> ef = EFetch(id='459567', rettype='fasta')
>>> lines = ef.read().splitlines()
>>> for line in lines:
...     print(line[:40])
...
>gi|459567|dbj|D28543.1|HPCNS5PC Hepatit
GAGCACGACATCTACCAATGTTGCCAACTGAACCCAGAGG
GGCTTTACCTTGGTGGTCCCATGTTTAACTCGCGAGGTCA
CGGGGTTCTTCCAACCAGCATGGGCAATACCCTCACATGT
GCAGGCCTCACCAATTCTGACATGTTGGTTTGCGGAGATG
TC
<BLANKLINE>

Similarly, if your id refers to a protein record, you can get that by setting
the ``rettype`` to 'gp'.

>>> genpept = EFetch(id='1234567,459567', rettype='gp').read()

You'll probably notice that the lines look suspiciously like FASTA-format
records. This is in fact true: the ``rettype`` parameter controls what type of
record you get back. For example, if we do the same search with
``rettype='brief'``, we get.

>>> from skbio.db.ncbi import EFetch
>>> ef = EFetch(id='459567', rettype='brief')
>>> lines = ef.read().splitlines()
>>> for line in lines:
...     print(line) # doctest: +SKIP
...
D28543 Hepatitis C virus... [gi:459567]

The current ``rettypes`` (as of this writing on 4/14/2010) for the 'core' NCBI
databases are native, fasta, gb, gp, gbwithparts, gbc, gpc, est, gss, seqid,
acc, ft. Formerly, but not currently, 'genbank' was a synonym for 'gb' and
'genpept' was a synonym for 'gp': however, these synonyms no longer work and
need to be fixed if you encounter them in old code. For more information check
NCBI's `format documentation
<http://www.ncbi.nlm.nih.gov/corehtml/query/static/efetchseq_help.html>`_.

Note that there are two separate concepts: rettype and retmode. rettype
controls what kind of data you'll get, and retmode controls how the data will
be formatted.

For example:

>>> from skbio.db.ncbi import EFetch
>>> ef = EFetch(id='459567', rettype='fasta', retmode='text')
>>> lines = ef.read().splitlines()
>>> for line in lines:
...     print(line[:40])
...
>gi|459567|dbj|D28543.1|HPCNS5PC Hepatit
GAGCACGACATCTACCAATGTTGCCAACTGAACCCAGAGG
GGCTTTACCTTGGTGGTCCCATGTTTAACTCGCGAGGTCA
CGGGGTTCTTCCAACCAGCATGGGCAATACCCTCACATGT
GCAGGCCTCACCAATTCTGACATGTTGGTTTGCGGAGATG
TC
<BLANKLINE>
>>> ef = EFetch(id='459567', rettype='fasta', retmode='html')
>>> lines = ef.read().splitlines()
>>> for line in lines:
...     print(line[:40]) # doctest: +SKIP
...
>gi|459567|dbj|D28543.1|HPCNS5PC Hepatit
GAGCACGACATCTACCAATGTTGCCAACTGAACCCAGAGG
GGCTTTACCTTGGTGGTCCCATGTTTAACTCGCGAGGTCA
CGGGGTTCTTCCAACCAGCATGGGCAATACCCTCACATGT
GCAGGCCTCACCAATTCTGACATGTTGGTTTGCGGAGATG
TC
<BLANKLINE>
>>> ef = EFetch(id='459567', rettype='fasta', retmode='xml')
>>> lines = ef.read().splitlines()
>>> for line in lines:
...     print(line[:40])
...
<?xml version="1.0"?>
 <!DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI
 <TSeqSet>
<TSeq>
  <TSeq_seqtype value="nucleotide"/>
  <TSeq_gi>459567</TSeq_gi>
  <TSeq_accver>D28543.1</TSeq_accver>
  <TSeq_taxid>11103</TSeq_taxid>
  <TSeq_orgname>Hepatitis C virus</TSeq_
  <TSeq_defline>Hepatitis C virus gene f
  <TSeq_length>282</TSeq_length>
  <TSeq_sequence>GAGCACGACATCTACCAATGTTG
</TSeq>
<BLANKLINE>
</TSeqSet>
<BLANKLINE>

You'll notice that the second case is some funny-looking html. Thanks, NCBI!
This is not our fault, please don't file a bug report. To figure out whether
something is actually surprising behavior at NCBI, you can always capture the
command-line and run it in a web browser. You can do this by calling str() on
the ``ef``, or by printing it. For example:

>>> print(ef) # doctest: +ELLIPSIS
http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi...

If you paste the resulting string into your web browser and you get the same
incorrect result that you get using PyCogent, you know that you should direct
your support requests NCBI's way. If you want to use your own email address
instead of leaving it as the default (the module developer), you can do that
just by passing it in as a parameter. For example, in the unlikely event that I
want NCBI to contact me instead of Mike if something goes wrong with my script,
I can achieve that as follows:

>>> ef = EFetch(id='459567', rettype='fasta', retmode='xml',
...             email='rob@spot.colorado.edu')
>>> print(ef) # doctest: +ELLIPSIS
http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi...

You can also select multiple ids (pass in as comma-delimited list):

>>> ef = EFetch(id='459567,459568', rettype='brief')
>>> ef.read() # doctest: +SKIP
'D28543 Hepatitis C virus... [gi:459567]\nBAA05896 NS5 protein [Hepa...8]'

Retrieving GenPept files from NCBI via Eutils
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We query for just one accession to illustrate the process. A more general query
can be executed by replacing ``'BAB52044`` with ``'"lysyl tRNA-synthetase"[ti]
AND bacteria[orgn]'`` in the snippet below.

>>> from skbio.db.ncbi import EUtils
>>> e = EUtils(numseqs=100, db='protein', rettype='gp')
>>> result = e['BAB52044']
>>> print(result.read()) # doctest: +ELLIPSIS
LOCUS       BAB52044                 548 aa            linear   BCT 16-MAY-2009
DEFINITION  lysyl tRNA synthetase [Mesorhizobium loti MAFF303099].
ACCESSION   BAB52044
VERSION     BAB52044.1  GI:14025444...

The EUtils modules are generic, so additional databases like OMIM can be
accessed using similar mechanisms.

Retrieving PubMed abstracts from NCBI via EUtils:

>>> from skbio.db.ncbi import EUtils
>>> e = EUtils(db='pubmed',rettype='brief')
>>> result = e['Simon Easteal AND Von Bing Yap'].read()
>>> print(result) # doctest: +ELLIPSIS
<BLANKLINE>
1. Mol Biol Evol. 2010 Mar;27(3):726-34. doi: 10.1093/molbev/msp232...
<BLANKLINE>
2. BMC Bioinformatics. 2008 Dec 19;9:550. doi: 10.1186/1471-2105-9...
<BLANKLINE>

Retrieving PubMed abstracts via PMID:

>>> from skbio.db.ncbi import EUtils
>>> e = EUtils(db='pubmed',rettype='abstract')
>>> result = e['14983078'].read()
>>> print(result) # doctest: +ELLIPSIS
<BLANKLINE>
1. Protein Eng. 2003 Dec;16(12):979-85.
<BLANKLINE>
Loops In Proteins (LIP)--a comprehensive loop database for...

"""

from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from numpy.testing import Tester
test = Tester().test
