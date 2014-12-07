r"""
E-utilities core objects (:mod:`skbio.db.ncbi`)
===============================================

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

from time import sleep
from xml.dom.minidom import parseString
from xml.etree.ElementTree import parse

from six import StringIO

from ._exception import QueryNotFoundError
from skbio.db.base import (URLGetter,
                           make_lists_of_expanded_slices_of_set_size,
                           make_lists_of_accessions_of_set_size)
from skbio.parse.record_finder import DelimitedRecordFinder


EUTILS_BASE_URL = 'http://www.ncbi.nlm.nih.gov/entrez/eutils'

# EUtils requires a tool and and email address
DEFAULT_TOOL_STRING = 'scikit-bio'
DEFAULT_EMAIL_ADDRESS = 'foo@bar.com'

VALID_DATABASES = {"pubmed", "protein", "nucleotide", "structure", "genome",
                   "books", "cancerchromosomes", "cdd", "domains", "gene",
                   "genomeprj", "gensat", "geo", "gds", "homologene",
                   "journals", "mesh", "ncbisearch", "nlmcatalog", "omim",
                   "pmc", "popset", "probe", "pcassay", "pccompound",
                   "pcsubstance", "snp", "taxonomy", "unigene", "unists"}

# Updated on December 6, 2014 from:
# http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_
# of__retmode_and/
RETTYPES = {}
RETTYPES['bioproject'] = {'docsum', 'uilist', 'xml'}
RETTYPES['biosample'] = {'docsum', 'uilist', 'full'}
RETTYPES['biosystems'] = {'docsum', 'uilist', 'xml'}
RETTYPES['gds'] = {'docsum', 'uilist', 'summary'}
RETTYPES['gene'] = {'docsum', 'uilist', 'gene_table', ''}
RETTYPES['homologene'] = {'docsum', 'uilist', 'alignmentscores', 'fasta',
                          'homologene'}
RETTYPES['mesh'] = {'docsum', 'uilist', 'full'}
RETTYPES['nlmcatalog'] = {'docsum', 'uilist'}
RETTYPES['nuccore'] = {'docsum', 'uilist', 'native', 'acc', 'fasta', 'seqid',
                       'gb', 'gbc', 'ft', 'gbwithparts', 'fasta_cds_na'}
RETTYPES['nucest'] = {'docsum', 'uilist', 'native', 'acc', 'fasta', 'seqid',
                      'gb', 'gbc', 'est'}
RETTYPES['nucgss'] = {'docsum', 'uilist', 'native', 'acc', 'fasta', 'seqid',
                      'gb', 'gbc', 'gss'}
RETTYPES['protein'] = {'docsum', 'uilist', 'native', 'acc', 'fasta', 'seqid',
                       'ft', 'gp', 'gp', 'gpc', 'ipg'}
RETTYPES['popset'] = {'docsum', 'uilist', 'native', 'acc', 'fasta', 'seqid',
                      'gb', 'gbc'}
RETTYPES['pmc'] = {'docsum', 'uilist', 'medline'}
RETTYPES['pubmed'] = {'docsum', 'uilist', 'medline', 'uilist', 'abstract'}
RETTYPES['snp'] = {'docsum', 'uilist', 'flt', 'fasta', 'rsr', 'ssexemplar',
                   'chr', 'genxml', 'docset'}
RETTYPES['sra'] = {'docsum', 'uilist', 'full'}
RETTYPES['taxonomy'] = {'docsum', 'uilist'}


class ESearch(URLGetter):

    """URLGetter subclass that uses the `esearch` operator

    Performs an `esearch` operation, getting a list of ids from an arbitrary
    query.

    References
    ----------
    .. [1] http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch

    See Also
    --------
    skbio.db.base.URLGetter

    """
    printed_fields = {'db', 'usehistory', 'term', 'retmax', 'retstart', 'tool',
                      'email'}
    defaults = {'db': 'nucleotide', 'usehistory': 'y', 'retmax': 1000,
                'tool': DEFAULT_TOOL_STRING, 'email': DEFAULT_EMAIL_ADDRESS}
    base_url = EUTILS_BASE_URL + '/esearch.fcgi?'


class EFetch(URLGetter):

    """URLGetter subclass that uses the `efetch` operator

    Performs an `efetch` operation, getting a list of ids from an arbitrary
    query.

    References
    ----------
    .. [1] http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

    See Also
    --------
    skbio.db.base.URLGetter

    """
    printed_fields = {'db', 'rettype', 'retmode', 'query_key', 'WebEnv',
                      'retmax', 'retstart', 'id', 'tool', 'email'}
    defaults = {'retmode': 'text', 'rettype': 'fasta', 'db': 'nucleotide',
                'retstart': 0, 'retmax': 100, 'tool': DEFAULT_TOOL_STRING,
                'email': DEFAULT_EMAIL_ADDRESS}
    base_url = EUTILS_BASE_URL + '/efetch.fcgi?'


class ELink(URLGetter):

    """URLGetter subclass that uses the `elink` operator

    Performs an `elink` operation, getting a list of ids from an arbitrary
    query (operates on databases that are linked together).

    References
    ----------
    .. [1] http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ELink

    See Also
    --------
    skbio.db.base.URLGetter

    """
    printed_fields = {'db', 'id', 'reldate', 'mindate', 'maxdate', 'datetype',
                      'term', 'retmode', 'db', 'dbfrom', 'WebEnv', 'query_key',
                      'holding', 'cmd', 'tool', 'email'}
    defaults = {'tool': DEFAULT_TOOL_STRING,
                'email': DEFAULT_EMAIL_ADDRESS}
    base_url = EUTILS_BASE_URL + '/elink.fcgi?'


class ESearchResult(object):

    """Container object for results retrieved on EUtils

    Notes
    -----
    This class is intended to be used internally ESearch.

    """

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __str__(self):
        return str(self.__dict__)


def _id_list_constructor(id_list_node):
    """Takes an id_list xml node and converts it into list of ids as strings"""
    return [_str_constructor(n) for n in id_list_node.childNodes
            if n.nodeType != n.TEXT_NODE]


def _int_constructor(node):
    """Makes an int out of node's first textnode child."""
    return int(node.firstChild.data)


def _str_constructor(node):
    """Makes an str out of node's first textnode child."""
    return node.firstChild.data

# the following are the only keys we explicitly handle now:
# (note difference in capitalization from parameters passed in)
ESEARCH_CONSTRUCTORS = {
    'Count': _int_constructor,
    'RetMax': _int_constructor,
    'RetStart': _int_constructor,
    'QueryKey': _int_constructor,
    'WebEnv': _str_constructor,
    'IdList': _id_list_constructor}


def _esearch_result_parser(result_as_string):
    """Parses an ESearch result. Returns ESearchResult object

    Parameters
    ----------
    result_as_string : str
        the string to be parsed

    Returns
    -------
    ESearchResult
        object representation of the parsed string
    """

    if '414 Request-URI Too Large' in result_as_string:
        raise ValueError("Tried to pass too large an URI:\n", result_as_string)

    doc = parseString(result_as_string)

    # assume one query result -- may need to fix
    query = doc.childNodes[-1]
    result = {}

    for n in query.childNodes:

        # skip top-level text nodes
        if n.nodeType == n.TEXT_NODE:
            continue
        name = n.tagName
        if name in ESEARCH_CONSTRUCTORS:
            result[name] = ESEARCH_CONSTRUCTORS[name](n)
        else:  # just keep the data if we don't know what it is
            result[name] = n.toxml()
    return ESearchResult(**result)


def _elink_result_parser(text):
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

    """Retrieves records from NCBI using the Entrez utilities

    Parameters
    ----------
    fp : str, optional
        The file path where results are written to. Defaults to ``None``.
    wait : float, optional
        Seconds to wait between fetching operations. Defaults to `0.5`.
    retmax : int, optional
        Maximum number of records per query. Defaults to `100`.
    url_limit : int, optional
        Limits the URL ESearch term size. Defaults to `400`.
    verbose : bool, optional
        Sets whether or not debugging information will be printed to standard
        output. Defaults to ``False``.
    max_recs : int, optional
        Total maximum number of records to retrieve. Defaults to None i. e. all
        records will be returned.

    See Also
    --------
    EFetch
    ELink
    ESearch
    skbio.db.base.URLGetter

    Examples
    --------
    >>> from skbio.db.ncbi import EUtils
    >>> eu = EUtils(db='pubmed', rettype='brief')
    >>> res = eu['PyCogent']
    >>> print(res.read()) # doctest: +ELLIPSIS
    <BLANKLINE>
    1...[PubMed]
    <BLANKLINE>
    <BLANKLINE>
    2...[PubMed - indexed for MEDLINE]
    <BLANKLINE>
    <BLANKLINE>
    3...[PubMed - indexed for MEDLINE]
    <BLANKLINE>
    <BLANKLINE>


    References
    ----------
    .. [1] http://www.ncbi.nlm.nih.gov/books/NBK25499/
    """

    def __init__(self, fp=None, wait=0.5, retmax=100, url_limit=400,
                 verbose=False, max_recs=None, **kwargs):
        self.__dict__.update(kwargs)
        self.fp = fp
        self.wait = wait
        self.retmax = retmax
        self.url_limit = url_limit
        self.max_recs = max_recs
        self.verbose = verbose
        self.term = None

        # adjust retmax if max_recs is set: no point getting more records
        if max_recs is not None and max_recs < retmax:
            self.retmax = max_recs

    def __getitem__(self, query):
        """Gets a query from NCBI. Assumes lists are lists of accessions.

        Returns
        -------
        file-like
            Returns a handle to the result (either in memory or file on disk).

        Notes
        -----
        The result is not guaranteed to contain any data.

        ..shownumpydoc
        """
        # check if it's a slice
        if isinstance(query, slice):
            queries = make_lists_of_expanded_slices_of_set_size(query)
            return self.grab_data(queries)

        # check if it's a list -- if so, delimit with ' '
        if isinstance(query, list) or isinstance(query, tuple):
            queries = make_lists_of_accessions_of_set_size(query)
            return self.grab_data(queries)

        # most likey a general set of search terms
        # e.g. '9606[taxid] OR 28901[taxid]' . So just return.
        return self.grab_data([query])

    def grab_data(self, queries):
        """Iterates through list of search terms and combines results.

        This will mostly only apply when the user wants to download thosands of
        sequences via accessions. This will supersede the GenBank URL length
        limit. So, we break up the accession list into sets of `400` terms per
        list.

        Parameters
        ----------
        queries : list of list objects
            list of lists of accession lists/query items

        Notes
        -----

        If you *really* have more than 300-400 terms similar to a query. The
        results will not be what you want anyway due do the limitations of the
        `esearch` URL length at GenBank. You'll just end up returning sets of
        results from the broken up word based search terms.
        """
        # figure out where to put the data
        if self.fp:
            result = open(self.fp, 'w')
        else:
            result = StringIO()

        for query in queries:
            self.term = query
            search_query = ESearch(**self.__dict__)

            # don't want the ids, just want to post search
            search_query.retmax = 0
            if self.verbose:
                print('SEARCH QUERY:')
                print(str(search_query))
            cookie = search_query.read()
            if self.verbose:
                print('COOKIE:')
                print(repr(cookie))
            search_result = _esearch_result_parser(cookie)
            if self.verbose:
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
                    raise QueryNotFoundError(
                        "WebEnv or query_key not Found! Query %s returned no "
                        "results.\nURL was:\n%s" %
                        (repr(query), str(search_query)))

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
                if self.verbose:
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
        if self.fp:
            result.close()
            return open(self.fp, 'r')
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
            search_result = _esearch_result_parser(cookie)
        else:
            search_result.IdList.extend(_esearch_result_parser(cookie).IdList)
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
    return _elink_result_parser(link.read())


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
            yield [element.strip() for element in between_tags.split(';')]

TAXON_RECORD_FINDER = DelimitedRecordFinder('</Taxon>', constructor=None,
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
            lineage = [e.strip() for e in _get_between_tags(line).split(';')]
    return taxid, name, lineage


def _get_taxa_names_lineages(lines):
    """Extracts taxon, name and lineage from each entry in an XML record."""
    empty_result = (None, None, None)
    for rec in TAXON_RECORD_FINDER(lines):
        curr = _get_taxid_name_lineage(rec)
        if curr != empty_result:
            yield curr


def _parse_taxonomy_using_elementtree_xml_parse(search_result):
    """Returns upper level XML taxonomy information from GenBank.
        search_result: StringIO object

    Returns
    -------
    list
        Returns list of all results in the form of:
        `[{result_01},{result_02},{result_03}]`.
        For each dict the key and values would be:
        key,value = xml label, e.g. [{'Lineage':'Bacteria; Proteobacteria...',
        'TaxId':'28901', 'ScientificName':'Salmonella enterica'}, {...}...]`
    """
    xml_data = parse(search_result)
    xml_data_root = xml_data.getroot()

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

    Notes
    -----

    Resulting lineages aren't in the same order as input. Use
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

    Notes
    -----

    Resulting lineages aren't in the same order as input. Use
    `taxon_ids_to_name_and_lineage` if you need the names and/or lineages
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
    primary_ids = _get_primary_ids(query, db=db)
    txon_ids = _ids_to_taxon_ids(primary_ids, db=db)
    return set([tuple(e) for e in _taxon_ids_to_lineages(txon_ids)])


def _get_unique_taxa(query, db='protein'):
    """Gets the unique lineages directly from a query."""
    return set(
        _taxon_ids_to_names(_ids_to_taxon_ids(_get_primary_ids(query, db=db),
                                              db=db)))
