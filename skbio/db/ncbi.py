r"""
E-utilities automation module (:mod:`skbio.db.ncbi`)
====================================================

.. currentmodule:: skbio.db.ncbi

This module provides a few core objects and functions to automate the retrieval
of results using the Entrez Programming Utilities [1]_.

Classes
-------

.. autosummary::
   :toctree: generated/

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

from urllib import urlopen, urlretrieve
from time import sleep
from xml.dom.minidom import parseString
from xml.etree.ElementTree import parse
from StringIO import StringIO
from string import strip

from skbio.core.exception import QueryNotFoundError
from skbio.db.util import (UrlGetter, expand_slice,
                           make_lists_of_expanded_slices_of_set_size,
                           make_lists_of_accessions_of_set_size)
from skbio.parse.record_finder import DelimitedRecordFinder


# eutils_base='http://eutils.ncbi.nlm.nih.gov/entrez/eutils'
eutils_base = 'http://www.ncbi.nlm.nih.gov/entrez/eutils'

# EUtils requires a tool and and email address
default_tool_string = 'scikit-bio'
default_email_address = 'yoshiki.vazquezbaeza@colorado.edu'

# databases last updated 7/22/05
valid_databases = dict.fromkeys(["pubmed", "protein", "nucleotide",
                                 "structure", "genome", "books",
                                 "cancerchromosomes", "cdd", "domains", "gene",
                                 "genomeprj", "gensat", "geo", "gds",
                                 "homologene", "journals", "mesh",
                                 "ncbisearch", "nlmcatalog", "omim", "pmc",
                                 "popset", "probe", "pcassay", "pccompound",
                                 "pcsubstance", "snp", "taxonomy", "unigene",
                                 "unists"])

# rettypes last updated 7/22/05 somehow, I don't think we'll be writing parsers
# for all these...
# WARNING BY RK 4/13/09: THESE RETTYPES ARE HIGHLY MISLEADING AND NO LONGER
# WORK. See this URL for the list of "official" rettypes, which is highly
# incomplete and has some important omissions (e.g. rettype 'gi' is missing but
# is the "official" replacement for 'GiList'):
# http://eutils.ncbi.nlm.nih.gov/entrez/query/static/efetchseq_help.html In
# particular, use gb or gp for GenBank or GenPept, use gi for GiList, use fasta
# for FASTA, and several other changes.  Until we get a complete accounting of
# what all the changes are, treat the rettypes below with extreme caution and
# experiment in the interpreter.
rettypes = {}
rettypes['pubmed'] = ('DocSum Brief Abstract Citation MEDLINE XML uilist '
                      'ExternalLink ASN1 pubmed_pubmed pubmed_pubmed_refs '
                      'pubmed_books_refs pubmed_cancerchromosomes pubmed_cdd '
                      'pubmed_domains pubmed_gds pubmed_gene pubmed_gene_rif '
                      'pubmed_genome pubmed_genomeprj pubmed_gensat pubmed_geo'
                      ' pubmed_homologene pubmed_nucleotide pubmed_omim '
                      ' pubmed_pcassay pubmed_pccompound '
                      'pubmed_pccompound_mesh '
                      'pubmed_pcsubstance pubmed_pcsubstance_mesh pubmed_pmc '
                      'pubmed_pmc_refs pubmed_popset '
                      'pubmed_probe pubmed_protein '
                      'pubmed_snp pubmed_structure '
                      'pubmed_unigene pubmed_unists')

rettypes['protein'] = ('DocSum ASN1 FASTA XML GenPept GiList graph fasta_xml '
                      'igp_xml gpc_xml ExternalLink protein_protein '
                      'protein_cdd protein_domains protein_gene protein_genome'
                      ' protein_genomeprj protein_homologene '
                      'protein_nucleotide protein_nucleotide_mgc protein_omim'
                      ' protein_pcassay protein_pccompound protein_pcsubstance'
                      ' protein_pmc protein_popset protein_pubmed protein_snp'
                      ' protein_snp_genegenotype protein_structure '
                      'protein_taxonomy protein_unigene')

rettypes['nucleotide'] = ('DocSum ASN1 FASTA XML GenBank GiList graph '
                          'fasta_xml gb_xml gbc_xml ExternalLink '
                          'nucleotide_comp_nucleotide nucleotide_nucleotide '
                          'nucleotide_nucleotide_comp '
                          'nucleotide_nucleotide_mrna nucleotide_comp_genome '
                          'nucleotide_gene nucleotide_genome '
                          'nucleotide_genome_samespecies nucleotide_gensat '
                          'nucleotide_geo nucleotide_homologene '
                          'nucleotide_mrna_genome nucleotide_omim '
                          'nucleotide_pcassay nucleotide_pccompound '
                          'nucleotide_pcsubstance nucleotide_pmc '
                          'nucleotide_popset nucleotide_probe '
                          'nucleotide_protein nucleotide_pubmed nucleotide_snp'
                          ' nucleotide_snp_genegenotype nucleotide_structure'
                          ' nucleotide_taxonomy nucleotide_unigene '
                          'nucleotide_unists')

rettypes['structure'] = ('DocSum Brief Structure Summary uilist ExternalLink'
                         ' structure_domains structure_genome '
                         'structure_nucleotide structure_omim '
                         'structure_pcassay structure_pccompound '
                         'structure_pcsubstance structure_pmc '
                         'structure_protein structure_pubmed structure_snp '
                         'structure_taxonomy')

rettypes['genome'] = ('DocSum ASN1 GenBank XML ExternalLink genome_genomeprj '
                      'genome_nucleotide genome_nucleotide_comp '
                      'genome_nucleotide_mrna genome_nucleotide_samespecies '
                      'genome_omim genome_pmc genome_protein genome_pubmed '
                      'genome_structure genome_taxonomy')

rettypes['books'] = ('DocSum Brief Books books_gene books_omim books_pmc_refs '
                     'books_pubmed_refs')

rettypes['cancerchromosomes'] = ('DocSum SkyCghDetails SkyCghCommon '
                                 'SkyCghCommonVerbose '
                                 'cancerchromosomes_cancerchromosomes_casecell'
                                 ' cancerchromosomes_cancerchromosomes_cellca'
                                 'se cancerchromosomes_cancerchromosomes_cytoc'
                                 'gh cancerchromosomes_cancerchromosomes_cytoc'
                                 'lincgh cancerchromosomes_cancerchromosomes_'
                                 'cytoclinsky cancerchromosomes_cancerchromos'
                                 'omes_cytodiagcgh cancerchromosomes_cancerch'
                                 'romosomes_cytodiagsky cancerchromosomes_can'
                                 'cerchromosomes_cytosky cancerchromosomes_ca'
                                 'ncerchromosomes_diag cancerchromosomes_canc'
                                 'erchromosomes_textual cancerchromosomes_pmc'
                                 ' cancerchromosomes_pubmed')

rettypes['cdd'] = ('DocSum Brief uilist cdd_cdd_fused cdd_cdd_related '
                   'cdd_gene cdd_homologene cdd_pmc cdd_protein cdd_pubmed '
                   'cdd_taxonomy')

rettypes['domains'] = ('DocSum Brief uilist domains_domains_new domains_pmc '
                       'domains_protein domains_pubmed domains_structure '
                       'domains_taxonomy')

rettypes['gene'] = ('Default DocSum Brief ASN.1 XML Graphics gene_table uilist'
                    ' ExternalLink gene_books gene_cdd gene_gensat gene_geo '
                    'gene_homologene gene_nucleotide gene_nucleotide_mgc '
                    'gene_omim gene_pmc gene_probe gene_protein gene_pubmed '
                    'gene_pubmed_rif gene_snp gene_snp_genegenotype '
                    'gene_taxonomy gene_unigene gene_unists')

rettypes['genomeprj'] = ('DocSum Brief Overview genomeprj_genomeprj '
                         'genomeprj_genome genomeprj_nucleotide '
                         'genomeprj_nucleotide_mrna genomeprj_nucleotide_orga'
                         'nella genomeprj_nucleotide_wgs genomeprj_pmc '
                         'genomeprj_popset genomeprj_protein genomeprj_pubmed '
                         'genomeprj_taxonomy')

rettypes['gensat'] = ('Group Detail DocSum Brief gensat_gensat gensat_gene '
                      'gensat_geo gensat_nucleotide gensat_pmc gensat_pubmed '
                      'gensat_taxonomy gensat_unigene')

rettypes['geo'] = ('DocSum Brief ExternalLink geo_geo_homologs geo_geo_prof '
                   'geo_geo_seq geo_gds geo_gene geo_gensat geo_homologene '
                   'geo_nucleotide geo_omim geo_pmc geo_pubmed geo_taxonomy '
                   'geo_unigene')

rettypes['gds'] = ('DocSum Brief gds_gds gds_geo gds_pmc gds_pubmed '
                   'gds_taxonomy')

rettypes['homologene'] = ('DocSum Brief HomoloGene AlignmentScores '
                          'MultipleAlignment ASN1 XML FASTA '
                          'homologene_homologene homologene_cdd '
                          'homologene_gene homologene_geo homologene_nucleoti'
                          'de homologene_omim homologene_pmc homologene_prote'
                          'in homologene_pubmed homologene_snp homologene_snp'
                          '_genegenotype homologene_taxonomy '
                          'homologene_unigene')

rettypes['journals'] = ('DocSum full journals_PubMed journals_Protein '
                        'journals_Nucleotide journals_Genome journals_Popset '
                        'journals_PMC journals_nlmcatalog')

rettypes['mesh'] = 'Full DocSum Brief mesh_PubMed'

rettypes['ncbisearch'] = 'DocSum Brief Home+Page+View ncbisearch_ncbisearch'

rettypes['nlmcatalog'] = 'Brief DocSum XML Expanded Full Subject ExternalLink'

rettypes['omim'] = ('DocSum Detailed Synopsis Variants ASN1 XML ExternalLink '
                    'omim_omim omim_books omim_gene omim_genome omim_geo '
                    'omim_homologene omim_nucleotide omim_pmc omim_protein '
                    'omim_pubmed omim_snp omim_snp_genegenotype omim_structure'
                    ' omim_unigene omim_unists')

rettypes['pmc'] = ('DocSum Brief XML TxTree pmc_books_refs '
                   'pmc_cancerchromosomes pmc_cdd pmc_domains pmc_gds pmc_gene'
                   ' pmc_genome pmc_genomeprj pmc_gensat pmc_geo '
                   'pmc_homologene pmc_nucleotide pmc_omim pmc_pccompound '
                   'pmc_pcsubstance pmc_popset pmc_protein pmc_pubmed '
                   'pmc_refs_pubmed pmc_snp pmc_structure pmc_taxonomy '
                   'pmc_unists')

rettypes['popset'] = ('DocSum PS ASN1 XML GiList ExternalLink TxTree '
                      'popset_genomeprj popset_nucleotide popset_protein '
                      'popset_pubmed popset_taxonomy')

rettypes['probe'] = ('DocSum Brief ASN1 XML Probe probe_probe probe_gene '
                     'probe_nucleotide probe_pubmed probe_taxonomy')

rettypes['pcassay'] = ('DocSum Brief uilist pcassay_nucleotide '
                       'pcassay_pccompound pcassay_pccompound_active '
                       'pcassay_pccompound_inactive pcassay_pcsubstance '
                       'pcassay_pcsubstance_active '
                       'pcassay_pcsubstance_inactive pcassay_protein '
                       'pcassay_pubmed pcassay_structure')

rettypes['pccompound'] = ('Brief DocSum PROP SYNONYMS pc_fetch '
                          'pccompound_pccompound_pulldown '
                          'pccompound_pccompound_sameanytautomer_pulldown '
                          'pccompound_pccompound_sameconnectivity_pulldown '
                          'pccompound_pccompound_sameisotopic_pulldown '
                          'pccompound_pccompound_samestereochem_pulldown '
                          'pccompound_nucleotide pccompound_pcassay '
                          'pccompound_pcassay_active '
                          'pccompound_pcassay_inactive pccompound_pcsubstance '
                          'pccompound_pmc pccompound_protein pccompound_pubmed'
                          ' pccompound_pubmed_mesh pccompound_structure')

rettypes['pcsubstance'] = ('Brief DocSum PROP SYNONYMS pc_fetch IDLIST '
                           'pcsubstance_pcsubstance_pulldown '
                           'pcsubstance_pcsubstance_same_pulldown '
                           'pcsubstance_pcsubstance_sameanytautomer_pulldown '
                           'pcsubstance_pcsubstance_sameconnectivity_pulldow '
                           'pcsubstance_pcsubstance_sameisotopic_pulldown '
                           'pcsubstance_pcsubstance_samestereochem_pulldown '
                           'pcsubstance_mesh pcsubstance_nucleotide '
                           'pcsubstance_pcassay pcsubstance_pcassay_active '
                           'pcsubstance_pcassay_inactive '
                           'pcsubstance_pccompound pcsubstance_pmc '
                           'pcsubstance_protein pcsubstance_pubmed '
                           'pcsubstance_pubmed_mesh pcsubstance_structure')

rettypes['snp'] = ('DocSum Brief FLT ASN1 XML FASTA RSR ssexemplar CHR FREQXML'
                   ' GENB GEN GENXML DocSet Batch uilist GbExp ExternalLink '
                   'MergeStatus snp_snp_genegenotype snp_gene snp_homologene '
                   'snp_nucleotide snp_omim snp_pmc snp_protein snp_pubmed '
                   'snp_structure snp_taxonomy snp_unigene snp_unists')

rettypes['taxonomy'] = ('DocSum Brief TxUidList TxInfo XML TxTree ExternalLink'
                        ' taxonomy_protein taxonomy_nucleotide '
                        'taxonomy_structure taxonomy_genome taxonomy_gene '
                        'taxonomy_cdd taxonomy_domains taxonomy_gds '
                        'taxonomy_genomeprj taxonomy_gensat '
                        'taxonomy_homologene taxonomy_pmc taxonomy_popset '
                        'taxonomy_probe taxonomy_pubmed taxonomy_snp '
                        'taxonomy_unigene taxonomy_unists')

rettypes['unigene'] = ('DocSum Brief ExternalLink unigene_unigene '
                       'unigene_unigene_expression unigene_unigene_homologous '
                       'unigene_gene unigene_gensat unigene_geo '
                       'unigene_homologene unigene_nucleotide '
                       'unigene_nucleotide_mgc unigene_omim unigene_protein '
                       'unigene_pubmed unigene_snp unigene_snp_genegenotype '
                       'unigene_taxonomy unigene_unists')

rettypes['unists'] = ('DocSum Brief ExternalLink unists_gene unists_nucleotide'
                      ' unists_omim unists_pmc unists_pubmed unists_snp '
                      'unists_taxonomy unists_unigene')

# convert into dict of known rettypes for efficient lookups -- don't want to
# scan list every time.
for key, val in rettypes.items():
    rettypes[key] = dict.fromkeys(val.split())


class ESearch(UrlGetter):
    """Performs an ESearch, getting a list of ids from an arbitrary query."""
    PrintedFields = dict.fromkeys(['db', 'usehistory', 'term', 'retmax',
                                   'retstart', 'tool', 'email'])
    Defaults = {'db': 'nucleotide', 'usehistory': 'y', 'retmax': 1000,
                'tool': default_tool_string, 'email': default_email_address}
    BaseUrl = eutils_base + '/esearch.fcgi?'


class EFetch(UrlGetter):
    """Retrieves a list of primary ids.

    WARNING: retmax (the maximum return value) is only 3 by default, so you
    will only get 3 records (this is for testing purposes). You will probably
    want to increase for real searches.
    """
    PrintedFields = dict.fromkeys(['db',
                                   'rettype',
                                   'retmode',
                                   'query_key',
                                   'WebEnv',
                                   'retmax',
                                   'retstart',
                                   'id',
                                   'tool',
                                   'email'])
    Defaults = {'retmode': 'text', 'rettype': 'fasta', 'db': 'nucleotide',
                'retstart': 0, 'retmax': 100, 'tool': default_tool_string,
                'email': default_email_address}
    BaseUrl = eutils_base + '/efetch.fcgi?'


class ELink(UrlGetter):
    """Retrieves a list of ids from one db that link to another db."""
    PrintedFields = dict.fromkeys(['db',
                                   'id',
                                   'reldate',
                                   'mindate',
                                   'maxdate',
                                   'datetype',
                                   'term',
                                   'retmode',
                                   'db',
                                   'dbfrom',
                                   'WebEnv',
                                   'query_key',
                                   'holding',
                                   'cmd',
                                   'tool',
                                   'email'])
    Defaults = {'tool': default_tool_string, 'email': default_email_address}
    BaseUrl = eutils_base + '/elink.fcgi?'


class ESearchResult(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __str__(self):
        return str(self.__dict__)


def id_list_constructor(id_list_node):
    """Takes an id_list xml node and converts it into list of ids as strings"""
    return [str_constructor(n) for n in id_list_node.childNodes
            if n.nodeType != n.TEXT_NODE]


def int_constructor(node):
    """Makes an int out of node's first textnode child."""
    return int(node.firstChild.data)


def str_constructor(node):
    """Makes an str out of node's first textnode child."""
    return str(node.firstChild.data)

# the following are the only keys we explicitly handle now:
#(note difference in capitalization from parameters passed in)
esearch_constructors = {
    'Count': int_constructor,
    'RetMax': int_constructor,
    'RetStart': int_constructor,
    'QueryKey': int_constructor,
    'WebEnv': str_constructor,
    'IdList': id_list_constructor}


def esearch_result_parser(result_as_string):
    """Parses an ESearch result. Returns ESearchResult object."""
    if '414 Request-URI Too Large' in result_as_string:
        raise ValueError("Tried to pass too large an URI:\n" +
                         result_as_string)
    doc = parseString(result_as_string)
    # assume one query result -- may need to fix
    query = doc.childNodes[-1]
    result = {}
    for n in query.childNodes:
        # skip top-level text nodes
        if n.nodeType == n.TEXT_NODE:
            continue
        name = str(n.tagName)  # who cares about unicode anyway...
        if name in esearch_constructors:
            result[name] = esearch_constructors[name](n)
        else:  # just keep the data if we don't know what it is
            result[name] = n.toxml()
    return ESearchResult(**result)


def elink_result_parser(text):
    """Gets the linked ids out of a single ELink result.

    Does not use the XML parser because of problems with long results.
    Only handles cases where there is a single set of links between
    databases.
    """
    result = []
    in_links = False
    for line in text.splitlines():
        if '<LinkName>' in line:
            in_links = True
        elif in_links and ('<Id>' in line):
            try:
                # expect line of form <Id>xxxx</Id>: want xxxx
                result.append(line.split('>', 1)[1].rsplit('<', 1)[0])
            except (IndexError, TypeError):
                pass
        elif '</LinkSetDb>' in line:  # end of block
            break
    return result


class EUtils(object):

    """Retrieves records from NCBI using EUtils."""

    def __init__(
            self,
            filename=None,
            wait=0.5,
            retmax=100,
            url_limit=400,
            DEBUG=False,
            max_recs=None,
            **kwargs):
        self.__dict__.update(kwargs)
        self.filename = filename
        self.wait = wait
        self.retstart = 0  # was originally set to 1
        self.DEBUG = DEBUG
        self.retmax = retmax
        self.url_limit = url_limit  # limits url esearch term size
        self.max_recs = max_recs
        # adjust retmax if max_recs is set: no point getting more records
        if max_recs is not None and max_recs < retmax:
            self.retmax = max_recs

    def __getitem__(self, query):
        """Gets an query from NCBI. Assumes lists are lists of accessions.

        Returns a handle to the result (either in memory or file on disk).
        WARNING: result is not guaranteed to contain any data.
        """
        # check if it's a slice
        if isinstance(query, slice):
            #query = expand_slice(query)
            queries = make_lists_of_expanded_slices_of_set_size(query)
            return self.grab_data(queries)

        # check if it's a list -- if so, delimit with ' '
        if isinstance(query, list) or isinstance(query, tuple):
            #query = ' '.join(map(str, query))
            queries = make_lists_of_accessions_of_set_size(query)
            return self.grab_data(queries)

        # most likey a general set of search terms
        # e.g. '9606[taxid] OR 28901[taxid]' . So just return.
        return self.grab_data([query])

    def grab_data(self, queries):
        """Iterates through list of search terms and combines results.
            -queries : list of lists of accession lists / query items

        This will mostly only apply whe the user wants to download 1000s of
        sequences via accessions. This will superced the GenBank url
        length limit. So, we break up the accession list into sets of 400
        terms per list.

        WARNING: if you _really_ have more than 300-400 terms similar to:
            'angiotensin[ti] AND rodents[orgn]'
            The results will not be what you want anyway due do the
            limitations of the esearch url length at GenBank. You'll just end
            up returning sets of results from the broken up word based search
            terms.
        """
        # figure out where to put the data
        if self.filename:
            result = open(self.filename, 'w')
        else:
            result = StringIO()

        for query in queries:
            self.term = query
            search_query = ESearch(**self.__dict__)
            # don't want the ids, just want to post search
            search_query.retmax = 0
            if self.DEBUG:
                print('SEARCH QUERY:')
                print(str(search_query))
            cookie = search_query.read()
            if self.DEBUG:
                print('COOKIE:')
                print(repr(cookie))
            search_result = esearch_result_parser(cookie)
            if self.DEBUG:
                print('SEARCH RESULT:')
                print(search_result)
            try:
                self.query_key = search_result.QueryKey
                self.WebEnv = search_result.WebEnv
            except AttributeError:
                # The query_key and/or WebEnv not Found!
                # GenBank occiasionally does not return these when user
                # attempts to only fetch data by Accession or UID. So we just
                # move on to extract UID list directly from the search result
                try:
                    self.id = ','.join(search_result.IdList)
                except AttributeError:
                    raise QueryNotFoundError("WebEnv or query_key not Found! "
                        "Query %s returned no results.\nURL was:\n%s" % (
                        repr(query), str(search_query)))

            count = search_result.Count

            # wrap the fetch in a loop so we get all the results
            fetch_query = EFetch(**self.__dict__)
            curr_rec = 0
            # check if we need to get additional ids

            if self.max_recs:  # cut off at max_recs if set
                count = min(count, self.max_recs)
                retmax = min(self.retmax, self.max_recs)
            else:
                retmax = self.retmax

            while curr_rec < count:
                # do the fetch
                if count - curr_rec < self.retmax:
                    fetch_query.retmax = count - curr_rec
                fetch_query.retstart = curr_rec
                if self.DEBUG:
                    print('FETCH QUERY')
                    print('CURR REC:', curr_rec, 'COUNT:', count)
                    print(str(fetch_query))
                # return the result of the fetch
                curr = fetch_query.read()
                result.write(curr)
                if not curr.endswith('\n'):
                    result.write('\n')
                curr_rec += retmax
                sleep(self.wait)
            # clean up after retrieval
        if self.filename:
            result.close()
            return open(self.filename, 'r')
        else:
            result.seek(0)
            return result

# The following are convenience wrappers for some of the above functionality

def _get_primary_ids(term, retmax=100, max_recs=None, **kwargs):
    """Gets primary ids from query."""
    search_result = None
    records_got = 0
    if max_recs:
        retmax = min(retmax, max_recs)
    search_query = ESearch(term=term, retmax=retmax, **kwargs)
    while True:
        cookie = search_query.read()
        if search_result is None:
            search_result = esearch_result_parser(cookie)
        else:
            search_result.IdList.extend(esearch_result_parser(cookie).IdList)
        # set the query key and WebEnv
        search_query.query_key = search_result.QueryKey
        search_query.WebEnv = search_result.WebEnv

        # if more results than retmax, keep adding results
        if max_recs:
            recs_to_get = min(max_recs, search_result.Count)
        else:
            recs_to_get = search_result.Count
        records_got += retmax
        if records_got >= recs_to_get:
            break
        elif recs_to_get - records_got < retmax:
            search_query.retmax = recs_to_get - records_got
        search_query.retstart = records_got
    return search_result.IdList


def _ids_to_taxon_ids(ids, db='nucleotide'):
    """Converts primary ids to taxon ids"""
    link = ELink(id=' '.join(ids), db='taxonomy', dbfrom=db, DEBUG=True)
    return elink_result_parser(link.read())


def _get_between_tags(line):
    """"Returns portion of line between xml tags."""
    return line.split('>', 1)[1].rsplit('<', 1)[0]


def _taxon_lineage_extractor(lines):
    """Extracts lineage from taxonomy record lines, not incl. species."""
    for line in lines:
        if '<Lineage>' in line:
            # expect line of form <Lineage>xxxx</Lineage> where xxxx semicolon-
            # delimited
            between_tags = line.split('>', 1)[1].rsplit('<', 1)[0]
            yield map(strip, between_tags.split(';'))

taxon_record_finder = DelimitedRecordFinder('</Taxon>', constructor=None,
                                            strict=False)


def _get_taxid_name_lineage(rec):
    """Returns taxon id, name, and lineage from single xml taxon record."""
    tax_tag = '  <TaxId>'
    name_tag = '  <ScientificName>'
    lineage_tag = '  <Lineage>'
    taxid = name = lineage = None
    for line in rec:
        if line.startswith(tax_tag):
            taxid = _get_between_tags(line)
        elif line.startswith(name_tag):
            name = _get_between_tags(line)
        elif line.startswith(lineage_tag):
            lineage = map(strip, _get_between_tags(line).split(';'))
    return taxid, name, lineage


def _get_taxa_names_lineages(lines):
    """Extracts taxon, name and lineage from each entry in an XML record."""
    empty_result = (None, None, None)
    for rec in taxon_record_finder(lines):
        curr = _get_taxid_name_lineage(rec)
        if curr != empty_result:
            yield curr


def _parse_taxonomy_using_elementtree_xml_parse(search_result):
    """Returns upper level XML taxonomy information from GenBank.
        search_result: StringIO object

    Returns list of all results in the form of:
        [{result_01},{result_02},{result_03}]
    For each dict the key and values would be:
        key,value = xml label, e.g. [{'Lineage':'Bacteria; Proteobacteria...',
                    'TaxId':'28901', 'ScientificName':'Salmonella enterica'},
                    {...}...]
    """
    xml_data = parse(search_result)
    xml_data_root = xml_data.getroot()

    tax_info_list = ['']

    l = []
    for individual_result in xml_data_root:
        children = list(individual_result)
        d = {}
        for child in children:
            key = child.tag
            value = child.text.strip()
            # We only want to retain the upper-level taxonomy information
            # from the xml parser and ignore all the rest of the information.
            # May revisit this in the future so that we can extract
            # 'GeneticCode', 'GCId', 'GCName', etc... <-- These values at this
            # level have whitespace, so we just ignore. Must traverse deeper to
            # obtain this information. Again, may implement in the future if
            # needed
            if value == '':
                continue
            else:
                d[key] = value
        l.append(d)
    return l


def _taxon_ids_to_names_and_lineages(ids, retmax=1000):
    """Yields taxon id, name and lineage for a set of taxon ids."""
    e = EUtils(db='taxonomy', rettype='xml', retmode='xml', retmax=retmax,
               DEBUG=False)
    fids = _fix_taxon_ids(ids)
    result = StringIO()
    result.write(e[fids].read())
    result.seek(0)
    data = _parse_taxonomy_using_elementtree_xml_parse(result)
    return [(i['TaxId'], i['ScientificName'], i['Lineage'])for i in data]


def _taxon_ids_to_lineages(ids, retmax=1000):
    """Returns full taxonomy (excluding species) from set of taxon ids.

    WARNING: Resulting lineages aren't in the same order as input. Use
    taxon_ids_to_name_and_lineage if you need the names and/or lineages
    associated with the specific ids.
    """
    ids = _fix_taxon_ids(ids)
    e = EUtils(db='taxonomy', rettype='xml', retmode='xml', retmax=retmax,
               DEBUG=False)
    result = e[ids].read().splitlines()
    return _taxon_lineage_extractor(result)


def _taxon_ids_to_names(ids, retmax=1000):
    """Returns names (e.g. species) from set of taxon ids.

    WARNING: Resulting lineages aren't in the same order as input. Use
    taxon_ids_to_name_and_lineage if you need the names and/or lineages
    associated with the specific ids.
    """
    e = EUtils(db='taxonomy', rettype='xml', retmode='xml', retmax=retmax,
               DEBUG=False)
    transformed_ids = _fix_taxon_ids(ids)
    h = StringIO()
    h.write(e[transformed_ids].read())
    h.seek(0)
    result = _parse_taxonomy_using_elementtree_xml_parse(h)
    return [i['ScientificName'] for i in result]


def _fix_taxon_ids(ids):
    """Fixes list of taxonomy ids by adding [taxid] to each.

    Need to add taxid field restriction to each id because NCBI broke taxon
    id search around 3/07 and has no plans to fix it.
    """
    if isinstance(ids, str):
        if not ids.endswith('[taxid]'):
            ids += '[taxid]'
        transformed_ids = ids
    else:
        transformed_ids = []
        for i in ids:
            if not i.endswith('[taxid]'):
                i = i.strip() + '[taxid]'
                transformed_ids.append(i)
    transformed_ids = ' OR '.join(transformed_ids)
    return transformed_ids


def _get_unique_lineages(query, db='protein'):
    """Gets the unique lineages directly from a query."""
    return set(map(tuple, _taxon_ids_to_lineages(_ids_to_taxon_ids(
        _get_primary_ids(query, db=db), db=db))))


def _get_unique_taxa(query, db='protein'):
    """Gets the unique lineages directly from a query."""
    return set(
        _taxon_ids_to_names(_ids_to_taxon_ids(_get_primary_ids(query, db=db),
                                            db=db)))

if __name__ == '__main__':
    main()
