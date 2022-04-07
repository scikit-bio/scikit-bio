import json
import time
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from urllib.error import HTTPError
from xml.etree import ElementTree

# Identification parameters used by NCBI
tool = None
email = None
api_key = None
# Time tracker to ensure that an excessive amount of calls per second is not issued
_latest_request_time = 0


class EntrezError(Exception):
    """Raised when the Entrez database throws an error in response to a request."""
    def __init__(self, msg):
        msg = json.loads(msg)
        super().__init__(msg["error"])


def einfo(db=None,
          version=None,
          retmode="xml",
          parse_mode="utf8",
          max_tries=2,
          interval=5):
    """Wrapper for Entrez's EInfo function.

    Functions:

    - Provides a list of the names of all valid Entrez databases
    - Provides statistics for a single database, including lists of indexing fields and
      available link names

    Parameters
    ----------
    db : str, optional
        Target database about which to gather statistics. Value must be a valid Entrez database
        name. If no `db` parameter is provided, einfo will return a list of the names of all valid
        Entrez databases.
    version : float, optional
        Used to specify version 2.0 EInfo XML. The only supported value is '2.0'.
    retmode : str, optional
        Retrieval type. Determines the format of the returned output. The default value is 'xml' for
        EInfo XML, but 'json' is also supported to return output in JSON format.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``).
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.

    Examples
    --------
    Return a list of all Entrez database names parsed into a python ``dict``:

    >>> einfo_result = einfo(retmode="json", parse_mode="json")

    Return version 2.0 statistics for Entrez Protein:

    >>> einfo_result = einfo("protein", version=2)

    Notes
    -----
    To find more details on how to use the EInfo function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EInfo
    """
    parameters = {
        "db": db,
        "version": version,
        "retmode": retmode
    }
    return _make_request("einfo.fcgi",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def esearch(db,
            term,
            usehistory=None,
            webenv=None,
            query_key=None,
            retstart=0,
            retmax=20,
            rettype="uilist",
            retmode="xml",
            sort=None,
            field=None,
            idtype=None,
            datetype=None,
            reldate=None,
            mindate=None,
            maxdate=None,
            parse_mode="utf8",
            max_tries=2,
            interval=5):
    """Wrapper for Entrez's ESearch function.

    Functions:

    - Provides a list of UIDs matching a text query
    - Posts the results of a search on the History server
    - Downloads all UIDs from a dataset stored on the History server
    - Combines or limits UID datasets stored on the History server
    - Sorts sets of UIDs

    Parameters
    ----------
    db : str
        Database to search. Value must be a valid Entrez database name.
    term : str
        Entrez text query. All special characters must be URL encoded. Spaces may be replaced by '+'
        signs. See the PubMed or Entrez help for information about search field descriptions and
        tags. Search fields and tags are database specific.
    usehistory : str, optional
        When `usehistory` is set to 'y', ESearch will post the UIDs resulting from the search
        operation onto the History server so that they can be used directly in a subsequent
        E-utility call. Also, `usehistory` must be set to 'y' for ESearch to interpret query key
        values included in `term` or to accept a `webenv` as input.
    webenv : str, optional
        Web environment string returned from a previous ESearch, EPost or ELink call. When provided,
        ESearch will post the results of the search operation to this pre-existing WebEnv, thereby
        appending the results to the existing environment. In addition, providing `webenv` allows
        query keys to be used in term so that previous search sets can be combined or limited. As
        described above, if `webenv` is used, `usehistory` must be set to 'y'.
    query_key : int, optional
        Integer query key returned by a previous ESearch, EPost or ELink call. When provided,
        ESearch will find the intersection of the set specified by `query_key` and the set retrieved
        by the query in `term` (i.e. joins the two with AND). For `query_key` to function, `webenv`
        must be assigned an existing WebEnv string and `usehistory` must be set to 'y'.
    retstart : int, optional
        Sequential index of the first UID in the retrieved set to be shown in the XML output
        (default=0, corresponding to the first record of the entire set). This parameter can be
        used in conjunction with `retmax` to download an arbitrary subset of UIDs retrieved from a
        search.
    retmax : int, optional
        Total number of UIDs from the retrieved set to be shown in the XML output (default=20). By
        default, ESearch only includes the first 20 UIDs retrieved in the XML output. If `usehistory`
        is set to 'y', the remainder of the retrieved set will be stored on the History server;
        otherwise these UIDs are lost. Increasing `retmax` allows more of the retrieved UIDs to be
        included in the XML output, up to a maximum of 100,000 records. To retrieve more than
        100,000 UIDs, submit multiple ESearch requests while incrementing the value of `retstart`.
    rettype : str, optional
        Retrieval type. There are two allowed values for ESearch: 'uilist' (default), which displays
        the standard XML output, and 'count', which displays only the <Count> tag.
    retmode : str, optional
        Retrieval type. Determines the format of the returned output. The default value is 'xml' for
        ESearch XML, but 'json' is also supported to return output in JSON format.
    sort : str, optional
        Specifies the method used to sort UIDs in the ESearch output. The available values vary by
        database (`db`) and may be found in the Display Settings menu on an Entrez search results
        page. If `usehistory` is set to 'y', the UIDs are loaded onto the History Server in the
        specified sort order and will be retrieved in that order by ESummary or EFetch. Example
        values are 'relevance' and 'name' for Gene and 'first+author' and 'pub+date' for PubMed.
        Users should be aware that the default value of `sort` varies from one database to another,
        and that the default value used by ESearch for a given database may differ from that used on
        NCBI web search pages.
    field : str, optional
        Search field. If used, the entire search term will be limited to the specified Entrez field.
    idtype : str, optional
        Specifies the type of identifier to return for sequence databases (nuccore, nucest, nucgss,
        popset, protein). By default, ESearch returns GI numbers in its output. If `idtype` is set
        to 'acc', ESearch will return accession.version identifiers rather than GI numbers.
    datetype : str, optional
        Type of date used to limit a search. The allowed values vary between Entrez databases, but
        common values are 'mdat' (modification date), 'pdat' (publication date) and 'edat'
        (Entrez date). Generally an Entrez database will have only two allowed values for `datetype`.
    reldate : int, optional
        When `reldate` is set to an integer n, the search returns only those items that have a date
        specified by `datetype` within the last n days.
    mindate : str, optional
        Minimum date of the range used to limit a search result by the date specified by `datetype`.
        The general date format is YYYY/MM/DD, and these variants are also allowed: YYYY, YYYY/MM.
        Musty be used along `maxdate`.
    maxdate : str, optional
        Maximum date of the range used to limit a search result by the date specified by `datetype`.
        The general date format is YYYY/MM/DD, and these variants are also allowed: YYYY, YYYY/MM.
        Musty be used along `mindate`.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``). This parameter is most
        useful when combined with `retmode`.
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.

    Examples
    --------
    Search in PubMed with the term cancer for abstracts that have an Entrez date within the last 60
    days; retrieve the first 100 PMIDs and translations; post the results on the History server and
    return a `webenv` and `query_key`:

    >>> search_result = esearch("pubmed", "cancer", reldate=60, datetype="edat", retmax=100, usehistory='y')

    Search in Nucleotide for all tRNAs:

    >>> search_result = esearch("nucleotide", "biomol+trna[prop]")

    Search in Protein for a molecular weight range:

    >>> search_result = esearch("protein", "70000:90000[molecular+weight]")

    Notes
    -----
    To find more details on how to use the ESearch function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch

    """
    if type(mindate) != type(maxdate):
        raise ValueError("'mindate' and 'maxdate' parameters must be either both defined or both"
                         "undefined")
    if webenv is not None and usehistory is None:
        raise ValueError("to use 'webenv' the 'usehistory' parameter must be set to 'y'")
    if query_key is not None and webenv is None:
        raise ValueError("to use 'query_key' a value must be provided for 'webenv' and the "
                         "'usehistory' parameter must be set to 'y'")

    parameters = {
        "db": db,
        "term": term,
        "usehistory": usehistory,
        "webenv": webenv,
        "query_key": query_key,
        "retstart": retstart,
        "retmax": retmax,
        "rettype": rettype,
        "retmode": retmode,
        "sort": sort,
        "field": field,
        "idtype": idtype,
        "datetype": datetype,
        "reldate": reldate,
        "mindate": mindate,
        "maxdate": maxdate
    }
    return _make_request("esearch.fcgi",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def epost(db,
          id,
          webenv=None,
          parse_mode="utf8",
          max_tries=2,
          interval=5):
    """Wrapper for Entrez's EPost function.

    Functions:

    - Uploads a list of UIDs to the Entrez History server
    - Appends a list of UIDs to an existing set of UID lists attached to a Web Environment

    Parameters
    ----------
    db : str
        Database containing the UIDs in the input list. The value must be a valid Entrez database
        name.
    id : str or int or list of (str or int), optional
        UID list. Either a single UID, a ``str`` containing UIDs delimited by commas or an iterable
        object of multiple UIDs may be provided. All of the UIDs must be from the database specified
        by `db`. For sequence databases (nuccore, nucest, nucgss, popset, protein), the UID list may
        be a mixed list of GI numbers and accession.version identifiers. Note: When using
        accession.version identifiers, there is a conversion step that takes place that causes large
        lists of identifiers to time out, even when using POST. Therefore, we recommend batching
        these types of requests in sizes of about 500 UIDs or less, to avoid retrieving only a
        partial amount of records from your original POST input list.
    webenv : str, optional
        Web Environment. If provided, this parameter specifies the Web Environment that will receive
        the UID list sent by post. EPost will create a new query key associated with that Web
        Environment. Usually this WebEnv value is obtained from the output of a previous ESearch,
        EPost or ELink call. If no `webenv` parameter is provided, EPost will create a new Web
        Environment and post the UID list to query key 1.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``).
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.

    Examples
    --------
    Post records to PubMed:

    >>> epost_result = epost("pubmed", [11237011, 12466850])

    Notes
    -----
    To find more details on how to use the EPost function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EPost
    """
    parameters = {
        "db": db,
        "id": id,
        "WebEnv": webenv
    }
    return _make_request("epost.fcgi",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def esummary(db,
             id=None,
             query_key=None,
             webenv=None,
             retstart=0,
             retmax=None,
             retmode="xml",
             version=None,
             parse_mode="utf8",
             max_tries=2,
             interval=5):
    """Wrapper for Entrez's ESummary function.

    Functions:

    - Returns document summaries (DocSums) for a list of input UIDs
    - Returns DocSums for a set of UIDs stored on the Entrez History server

    Parameters
    ----------
    db : str
        Database from which to retrieve DocSums. The value must be a valid Entrez database name.
    id : str or int or list of (str or int), optional
        UID list. Required when input is from a UID list. Either a single UID, a ``str`` containing
        UIDs delimited by commas or an iterable object of multiple UIDs may be provided. All of the
        UIDs must be from the database specified by `db`. For sequence databases (nuccore, nucest,
        nucgss, popset, protein), the UID list may be a mixed list of GI numbers and
        accession.version identifiers.
    query_key : int, optional
        Query key. Required when input is from the Entrez History server. This integer specifies
        which of the UID lists attached to the given Web Environment will be used as input to
        ESummary. Query keys are obtained from the output of previous ESearch, EPost or ELInk calls.
        The `query_key` parameter must be used in conjunction with `webenv`.
    webenv : str, optional
        Web Environment. Required when input is from the Entrez History server. This parameter
        specifies the Web Environment that contains the UID list to be provided as input to ESummary.
        Usually this WebEnv value is obtained from the output of a previous ESearch, EPost or ELink
        call. The `webenv` parameter must be used in conjunction with `query_key`.
    retstart : int, optional
        Sequential index of the first DocSum to be retrieved (default=0, corresponding to the first
        record of the entire set). This parameter can be used in conjunction with `retmax` to
        download an arbitrary subset of DocSums from the input set.
    retmax : int, optional
        Total number of DocSums from the input set to be retrieved, up to a maximum of 10,000. If
        the total set is larger than this maximum, the value of `retstart` can be iterated while
        holding `retmax` constant, thereby downloading the entire set in batches of size `retmax`.
    retmode : str, optional
        Retrieval type. Determines the format of the returned output. The default value is 'xml' for
        ESummary XML, but 'json' is also supported to return output in JSON format.
    version : float, optional
        Used to specify version 2.0 ESummary XML. The only supported value is '2.0'. When present,
        ESummary will return version 2.0 DocSum XML that is unique to each Entrez database and that
        often contains more data than the default DocSum XML.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``).
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.

    Examples
    --------
    Retrieve PubMed DocSums, version 2.0 XML parsed into an ``ElementTree`` object:

    >>> esummary_result = esummary("pubmed", id=[11850928, 11482001], version=2)

    Retrieve nucleotide DocSums:

    >>> esummary_result = esummary("nucleotide", id=[28864546, 28800981])

    Notes
    -----
    To find more details on how to use the ESummary function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESummary
    """
    _ensure_single_input_source(id, query_key, webenv)
    parameters = {
        "db": db,
        "id": id,
        "query_key": query_key,
        "webenv": webenv,
        "retstart": retstart,
        "retmax": retmax,
        "retmode": retmode,
        "version": version
    }
    return _make_request("esummary.fcgi",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def efetch(db,
           id=None,
           query_key=None,
           webenv=None,
           retmode=None,
           rettype=None,
           retstart=0,
           retmax=None,
           strand=None,
           seq_start=None,
           seq_stop=None,
           complexity=None,
           parse_mode="utf8",
           max_tries=2,
           interval=5):
    """Wrapper for Entrez's EFetch function.

    Functions:

    - Returns formatted data records for a list of input UIDs
    - Returns formatted data records for a set of UIDs stored on the Entrez History server

    Parameters
    ----------
    db : str
        Database from which to retrieve records. The value must be a valid Entrez database name.
        Currently EFetch does not support all Entrez databases. Please see Table 1 in Chapter 2 of
        the Entrez documentation for a list of available databases.
    id : str or int or list of (str or int), optional
        UID list. Required when input is from a UID list. Either a single UID, a ``str`` containing
        UIDs delimited by commas or an iterable object of multiple UIDs may be provided. All of the
        UIDs must be from the database specified by `db`. For sequence databases (nuccore, nucest,
        nucgss, popset, protein), the UID list may be a mixed list of GI numbers and
        accession.version identifiers. Special note for sequence databases: NCBI is no longer
        assigning GI numbers to a growing number of new sequence records. As such, these records are
        not indexed in Entrez, and so cannot be retrieved using ESearch or ESummary, and have no
        Entrez links accessible by ELink. EFetch can retrieve these records by including their
        accession.version identifier in the `id` parameter.
    query_key : int, optional
        Query key. Required when input is from the Entrez History server. This integer specifies
        which of the UID lists attached to the given Web Environment will be used as input to
        EFetch. Query keys are obtained from the output of previous ESearch, EPost or ELInk calls.
        The `query_key` parameter must be used in conjunction with `webenv`.
    webenv : str, optional
        Web Environment. Required when input is from the Entrez History server. This parameter
        specifies the Web Environment that contains the UID list to be provided as input to EFetch.
        Usually this WebEnv value is obtained from the output of a previous ESearch, EPost or ELink
        call. The `webenv` parameter must be used in conjunction with `query_key`.
    retmode : str, optional
        Retrieval mode. This parameter specifies the data format of the records returned, such as
        plain text, HMTL or XML. See Table 1 on EFetch's documentation for a full list of allowed
        values for each database.
    rettype : str, optional
        Retrieval type. This parameter specifies the record view returned, such as Abstract or
        MEDLINE from PubMed, or GenPept or FASTA from protein. Please see Table 1 on EFetch's
        documentation for a full list of allowed values for each database.
    retstart : int, optional
        Sequential index of the first record to be retrieved (default=0, corresponding to the first
        record of the entire set). This parameter can be used in conjunction with `retmax` to
        download an arbitrary subset of records from the input set.
    retmax : int, optional
        Total number of records from the input set to be retrieved, up to a maximum of 10,000.
        Optionally, for a large set the value of `retstart` can be iterated while holding `retmax`
        constant, thereby downloading the entire set in batches of size `retmax`.
    strand : int, optional
        Strand of DNA to retrieve. Applicable only for sequence databases. Available values are "1"
        for the plus strand and "2" for the minus strand.
    seq_start : int, optional
        First sequence base to retrieve. Applicable only for sequence databases. The value should be
        the integer coordinate of the first desired base, with "1" representing the first base of
        the sequence.
    seq_stop : int, optional
        Last sequence base to retrieve. Applicable only for sequence databases. The value should be
        the integer coordinate of the last desired base, with "1" representing the first base of the
        sequence.
    complexity : int, optional
        Data content to return. Applicable only for sequence databases. Many sequence records are
        part of a larger data structure or "blob", and the `complexity` parameter determines how
        much of that blob to return. For example, an mRNA may be stored together with its protein
        product. The available values are as follows: 0 - entire blob; 1 - bioseq; 2 - minimal
        bioseq-set; 3 - minimal nuc-prot; 4 - minimal pub-set.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``).
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.

    Examples
    --------
    Fetch PMIDs in XML and parse the result into an ``ElementTree`` object:

    >>> efetch_result = efetch("pubmed", id=[11748933, 11700088], retmode="xml", parse_mode="xml")

    Fetch the first 100 bases of the plus strand of GI 21614549 in FASTA format:

    >>> efetch_result = efetch("nuccore", id=21614549, strand=1, seq_start=1, seq_stop=100,
        rettype="fasta", retmode="text")

    Notes
    -----
    To find more details on how to use the EFetch function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
    """
    _ensure_single_input_source(id, query_key, webenv)
    parameters = {
        "db": db,
        "id": id,
        "query_key": query_key,
        "webenv": webenv,
        "retmode": retmode,
        "rettype": rettype,
        "retstart": retstart,
        "retmax": retmax,
        "strand": strand,
        "seq_start": seq_start,
        "seq_stop": seq_stop,
        "complexity": complexity
    }
    return _make_request("efetch.fcgi",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def elink(db,
          dbfrom,
          cmd="neighbor",
          id=None,
          query_key=None,
          webenv=None,
          retmode="xml",
          idtype=None,
          linkname=None,
          term=None,
          holding=None,
          datetype=None,
          reldate=None,
          mindate=None,
          maxdate=None,
          parse_mode="utf8",
          max_tries=2,
          interval=5):
    """Wrapper for Entrez's ELink function.

    Functions:

    - Returns UIDs linked to an input set of UIDs in either the same or a different Entrez
      database
    - Returns UIDs linked to other UIDs in the same Entrez database that match an Entrez query
    - Checks for the existence of Entrez links for a set of UIDs within the same database
    - Lists the available links for a UID
    - Lists LinkOut URLs and attributes for a set of UIDs
    - Lists hyperlinks to primary LinkOut providers for a set of UIDs
    - Creates hyperlinks to the primary LinkOut provider for a single UID

    Parameters
    ----------
    db : str
        Database from which to retrieve UIDs. The value must be a valid Entrez database name.
        This is the destination database for the link operation.
    dbfrom : str
        Database containing the input UIDs. The value must be a valid Entrez database name. This is
        the origin database of the link operation. If `db` and `dbfrom` are set to the same
        database value, then ELink will return computational neighbors within that database. Please
        see the full list of Entrez links in the ELink official documentation for available
        computational neighbors. Computational neighbors have linknames that begin with
        dbname_dbname (examples: protein_protein, pcassay_pcassay_activityneighbor).
    cmd : str, optional
        ELink command mode. The command mode specified which function ELink will perform. Some
        optional parameters only function for certain values of `cmd`. Please refer to ELink's
        official documentation for a full list of possible values for `cmd`.
    id : str or int or list of (str or int), optional
        UID list. Required when input is from a UID list. Either a single UID, a ``str`` containing
        UIDs delimited by commas or an iterable object of multiple UIDs may be provided. All of the
        UIDs must be from the database specified by `dbfrom`. If more than one `id` parameter is
        provided, ELink will perform a separate link operation for the set of UIDs specified by each
        `id` parameter. This effectively accomplishes "one-to-one" links and preserves the
        connection between the input and output UIDs. For sequence databases (nuccore, nucest,
        nucgss, popset, protein), the UID list may be a mixed list of GI numbers and
        accession.version identifiers.
    query_key : int, optional
        Query key. Required when input is from the Entrez History server. This integer specifies
        which of the UID lists attached to the given Web Environment will be used as input to
        ELink. Query keys are obtained from the output of previous ESearch, EPost or ELInk calls.
        The `query_key` parameter must be used in conjunction with `webenv`.
    webenv : str, optional
        Web Environment. Required when input is from the Entrez History server. This parameter
        specifies the Web Environment that contains the UID list to be provided as input to ELink.
        Usually this WebEnv value is obtained from the output of a previous ESearch, EPost or ELink
        call. The `webenv` parameter must be used in conjunction with `query_key`.
    retmode : str, optional
        Retrieval type. Determines the format of the returned output. The default value is 'xml' for
        ELink XML, but 'json' is also supported to return output in JSON format.
    idtype : str, optional
        Specifies the type of identifier to return for sequence databases (nuccore, nucest, nucgss,
        popset, protein). By default, ELink returns GI numbers in its output. If `idtype` is set
        to 'acc', ELink will return accession.version identifiers rather than GI numbers.
    linkname : str, optional
        Name of the Entrez link to retrieve. Every link in Entrez is given a name of the form
        `dbfrom_db_subset`. The values of subset vary depending on the values of `dbfrom` and `db`.
        Many `dbfrom`/`db` combinations have no subset values. See the list of Entrez links in
        ELink's official documentation for a listing of all available linknames. When `linkname` is
        used, only the links with that name will be retrieved. The `linkname` parameter only
        functions when `cmd` is set to 'neighbor' or 'neighbor_history'.
    term : str, optional
        Entrez query used to limit the output set of linked UIDs. The query in the `term` parameter
        will be applied after the link operation, and only those UIDs matching the query will be
        returned by ELink. The `term` parameter only functions when `db` and `dbfrom` are set to the
        same database value.
    holding : str, optional
        Name of LinkOut provider. Only URLs for the LinkOut provider specified by `holding` will be
        returned. The value provided to `holding` should be the abbreviation of the LinkOut
        provider's name found in the <NameAbbr> tag of the ELink XML output when `cmd` is set to
        'llinks' or 'llinkslib'. The `holding` parameter only functions when `cmd` is set to
        'llinks' or 'llinkslib'.
    datetype : str, optional
        Type of date used to limit a link operation. This parameter only functions when `cmd` is set
        to 'neighbor' or 'neighbor_history' and `dbfrom` is 'pubmed'. The allowed values vary
        between Entrez databases, but common values are 'mdat' (modification date), 'pdat'
        (publication date) and 'edat' (Entrez date). Generally an Entrez database will have only two
        allowed values for `datetype`.
    reldate : int, optional
        When `reldate` is set to an integer n, ELink returns only those items that have a date
        specified by `datetype` within the last n days.This parameter only functions when `cmd` is
        set to 'neighbor' or 'neighbor_history' and `dbfrom` is 'pubmed'.
    mindate : str, optional
        Minimum date of the range used to limit a search result by the date specified by `datetype`.
        The general date format is YYYY/MM/DD, and these variants are also allowed: YYYY, YYYY/MM.
        Musty be used along `maxdate`.
    maxdate : str, optional
        Maximum date of the range used to limit a search result by the date specified by `datetype`.
        The general date format is YYYY/MM/DD, and these variants are also allowed: YYYY, YYYY/MM.
        Musty be used along `mindate`.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``).
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.

    Examples
    --------
    Link from protein to gene:

    >>> elink_result = elink("gene", "protein", id=[15718680, 157427902])

    Find all links from gene to snp:

    >>> elink_result = elink("snp", "gene", id=93986)

    Notes
    -----
    To find more details on how to use the ELink function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ELink
    """
    if type(mindate) != type(maxdate):
        raise ValueError("'mindate' and 'maxdate' parameters must be either both defined or both"
                         "undefined")
    _ensure_single_input_source(id, query_key, webenv)
    parameters = {
        "db": db,
        "dbfrom": dbfrom,
        "cmd": cmd,
        "id": id,
        "query_key": query_key,
        "webenv": webenv,
        "retmode": retmode,
        "idtype": idtype,
        "linkname": linkname,
        "term": term,
        "holding": holding,
        "datetype": datetype,
        "reldate": reldate,
        "mindate": mindate,
        "maxdate": maxdate
    }
    return _make_request("elink.fcgi",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def egquery(term,
            parse_mode="utf8",
            max_tries=2,
            interval=5,
            **parameters):
    """Wrapper for Entrez's EGQuery function.

    Function: Provides the number of records retrieved in all Entrez databases by a single text
    query.

    Parameters
    ----------
    term : str
        Entrez text query. All special characters must be URL encoded. Spaces may be replaced by '+'
        signs. See the PubMed or Entrez help for information about search field descriptions and
        tags. Search fields and tags are database specific.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``).
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.
    **parameters : dict, optional
        Although no additional parameters for EGQuery are described in the official documentation,
        some common keywords such as `retmode` seem to also work for this function. Such parameters
        may be passed here as a convenience.

    Examples
    --------
    Search for the term 'asthma' in all databases and return the number of hits as XML:

    >>> egquery_result = egquery("asthma", retmode="xml")

    Notes
    -----
    To find more details on how to use the EGQuery function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EGQuery
    """
    parameters["term"] = term
    return _make_request("gquery",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def espell(db,
           term,
           parse_mode="utf8",
           max_tries=2,
           interval=5):
    """Wrapper for Entrez's ESpell function.

    Function: Provides spelling suggestions for terms within a single text query in a given database.

    Parameters
    ----------
    db : str
        Database to search. Value must be a valid Entrez database name.
    term : str
        Entrez text query. All special characters must be URL encoded. Spaces may be replaced by '+'
        signs. See the PubMed or Entrez help for information about search field descriptions and
        tags. Search fields and tags are database specific.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``).
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.

    Examples
    --------
    Spelling suggestions from PubMed for the term 'asthmaa+OR+alergies':

    >>> espell_result = espell("pubmed", "asthmaa+OR+alergies")

    Notes
    -----
    To find more details on how to use the ESpell function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESpell
    """
    parameters = {
        "db": db,
        "term": term
    }
    return _make_request("espell.fcgi",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def ecitmatch(bdata,
              db="pubmed",
              retmode="xml",
              parse_mode="utf8",
              max_tries=2,
              interval=5):
    """Wrapper for Entrez's ECitMatch function.

    Function: Retrieves PubMed IDs (PMIDs) that correspond to a set of input citation strings.

    Parameters
    ----------
    bdata : str
        Citation strings. Each input citation must be represented by a citation string in the
        following format: 'journal_title|year|volume|first_page|author_name|your_key|'. Multiple
        citation strings may be provided by separating the strings with a carriage return character
        (%0D). The your_key value is an arbitrary label provided by the user that may serve as a
        local identifier for the citation, and it will be included in the output. Be aware that all
        spaces must be replaced by '+' symbols and that citation strings should end with a final
        vertical bar '|'.
    db : str, optional
       Database to search. The only supported value is 'pubmed'.
    retmode : str, optional
        Retrieval type. The only supported value is 'xml'.
    parse_mode : str, optional
        Indicates how to parse the result from the Entrez API call. Currently supported values are
        ``None`` (returns the raw ``bytes`` object retrieved from the database) 'xml' (returns an
        ``ElementTree`` object), 'json' (returns a ``dict`` object) and text encodings understood
        by :py:func:`bytes.decode`, such as 'utf8' (returns a ``str``).
    max_tries : int, optional
        How many times to try to get a response from the database server before throwing an error.
    interval : float, optional
        How long to wait between tries to get a response from the database server.

    Examples
    -------
    Search for a citation in PubMed:

    >>> ecitmatch_result = ecitmatch("pubmed", bdata="proc+natl+acad+sci+u+s+a|1991|88|3248|mann+bj|Art1|%0Dscience|1987|235|182|palmenberg+ac|Art2|")

    Notes
    -----
    To find more details on how to use the ECitMatch function, please refer to its official
    documentation [1]_.

    References
    ----------
    .. [1] https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ECitMatch
    """
    parameters = {
        "db": db,
        "retmode": retmode,
        "bdata": bdata
    }
    return _make_request("ecitmatch.cgi",
                         parameters,
                         parse_mode=parse_mode,
                         max_tries=max_tries,
                         interval=interval)


def _build_request(cgi, parameters):
    """Formats parameters and builds a ``Request`` object."""
    # Formats the parameters into the expected by the Entrez API
    # Check if id is present in a format understood by the API and formats it otherwise
    id_ = parameters.get("id", None)
    if id_ is not None and not isinstance(id_, (str, int)):
        parameters["id"] = ",".join(str(uid) for uid in parameters["id"])
    # Check if parameter for version was passed as integer or float
    if parameters.get("version", None) is not None:
        parameters["version"] = str(float(parameters["version"]))
    # Add all provided identification parameters
    parameters["tool"] = tool
    parameters["email"] = email
    parameters["api_key"] = api_key
    # Remove undefined parameters to avoid polluting the url
    parameters = {k: v for k, v in parameters.items() if v is not None}

    # Request is made as POST so that there are no problems regarding character limits
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + cgi
    return Request(url, data=urlencode(parameters).encode("ascii"), method="POST")


def _parse_response(result, parse_mode):
    """Parses the result from an Entrez function into `parse_mode`."""
    if parse_mode.lower() == "json":
        return json.loads(result)
    if parse_mode.lower() == "xml":
        return ElementTree.fromstring(result)
    try:
        return result.decode(parse_mode)
    # parse_mode was not found as an encoding
    except LookupError as e:
        raise ValueError("only 'json', 'xml' and text encodings such as 'utf8' are currently supported "
                         "values for parse_mode") from e


def _make_request(cgi, parameters, parse_mode="utf8", max_tries=2, interval=5.0):
    """Generic function to make requests to Entrez's servers."""
    global _latest_request_time

    request = _build_request(cgi, parameters)

    # NCBI has a threshold for calls per second
    min_delay = 0.1 if api_key is not None else 1 / 3
    elapsed_time = time.time() - _latest_request_time
    if elapsed_time < min_delay:
        time.sleep(min_delay - elapsed_time)
    for i in range(max_tries):
        try:
            with urlopen(request) as response:
                result = response.read()
            break
        except HTTPError as e:
            # Raises exception anyway if the code is that of BAD REQUEST or the last try has been
            # reached
            if e.code == 400 or i == max_tries - 1:
                raise EntrezError(e.read()) from e
            time.sleep(interval)
    # Updates the global request time tracker
    _latest_request_time = time.time()

    # Parse results into a python object or returns the raw bytestring
    if parse_mode is not None:
        return _parse_response(result, parse_mode)
    return result


def _ensure_single_input_source(id, query_key, webenv):
    """Logic to ensure a valid form of input was given (either from id or both webenv and query_key)"""
    webenv_input = sum((webenv is None, query_key is None))
    if id is None:
        if webenv_input == 2:
            raise ValueError("a form of input must be provided from either setting the 'id' "
                             "parameter or setting both the 'webenv' and 'query_key' parameters")
        elif webenv_input == 1:
            raise ValueError("'webenv' and 'query_key' parameters are both required when retrieving "
                             "the input from a Web Environment")
    if webenv_input != 2:
        raise ValueError("the 'id' parameter can't be defined together 'webenv' or 'query_key', "
                         "these forms of input are mutually exclusive")
