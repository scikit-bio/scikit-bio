r"""
Database exploration objects and functions (:mod:`skbio.db.base`)
=================================================================

.. currentmodule:: skbio.db.base

This module provides a base class to work with querying RESTful interfaces
like NCBI's Entrez utilities.

Classes
-------

.. autosummary::
   :toctree: generated/

   URLGetter

Functions
---------

.. autosummary::
    :toctree: generated/

    expand_slice
    last_nondigit_index
    make_lists_of_accessions_of_set_size
    make_lists_of_expanded_slices_of_set_size

"""
from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from urllib import urlopen, urlretrieve, quote_plus


class URLGetter(object):
    """Base class to submit queries with EUtils

    This class provides the base functionality to work with the Entrez database
    as provided by the National Center of Biotechnology Information (NCBI).

    Attributes
    ----------
    defaults : dict
        Default values for the URL. Defaults to an empty dictionary.
    printed_fields : dict
        Values to be included in the URL. Defaults to an empty dictionary.
    base_url : str
        Base URL for the class (should end in `/`). Defaults to an empty
        string.
    key_value_delimiter : str
        Character that delimits keys & values. Defaults to `=`.
    field_delimiter : str
        Character that separates key-value pairs. Defaults to `&`.
    """
    defaults = {}
    printed_fields = {}
    base_url = ''
    key_value_delimiter = '='
    field_delimiter = '&'


    def __init__(self, **kwargs):
        self.__dict__.update(self.defaults)
        self.__dict__.update(kwargs)
        self._temp_args = {}

    def __str__(self):
        """Returns the formatted URL for this object

        Returns
        -------
        str
            Formatted URL that the class represents.

        Examples
        --------
        >>> from skbio.db.base import URLGetter
        >>> a = URLGetter(search='foo', db='bar', public='n')
        >>> a.base_url = "http://www.google.com/"
        >>> a.printed_fields['search'] = None
        >>> a.printed_fields['db'] = None
        >>> print(a)
        http://www.google.com/search=foo&db=bar

        ..shownumpydoc
        """
        to_get = self.__dict__.copy()
        to_get.update(self._temp_args)
        fields = []

        for key, value in to_get.items():
            if key in self.printed_fields:
                fields.append(quote_plus(key)+self.key_value_delimiter+
                              quote_plus(str(value)))

        return self.base_url + self.field_delimiter.join(fields)

    def open(self, **kwargs):
        """Return a stream to the URL constructed by this class

        Returns
        -------
        file-like object
            Returns a file-like object from which information can be read
        """
        self._temp_args = kwargs
        result = urlopen(str(self))
        self._temp_args = {}
        return result

    def read(self, **kwargs):
        """Reads the contents of the URL constructed by this class

        Returns
        -------
        object
            Read data from the URL, this can be a string, bytes or the binary
            data contained in the URL.
        """
        result = self.open(**kwargs)
        data = result.read()
        result.close()
        return data

    def retrieve(self, fp, **kwargs):
        """Reads and writes to a file the contents of the URL

        Parameters
        ----------
        fp : str
            File path where the contents will be written to.

        Notes
        -----
        This method produces no return value.
        """
        self._temp_args = kwargs
        urlretrieve(str(self), fp)
        self._temp_args = None

def expand_slice(s):
    """Takes a start and end accession, and generate the whole range.

    Both accessions must have the same non-numeric prefix.

    Parameters
    ----------
    s : list
        List of accession numbers.

    Notes
    -----
    Unlike standard slices, includes the last item in the range. In other
    words, `obj[AF1001:AF1010]` will include `AF1010`.

    Raises
    ------
    TypeError
       If the accessions are prefixes i. e. `AF100:`.

    Examples
    --------
    >>> from skbio.db.base import expand_slice
    >>> expand_slice(slice('AF1001', 'AF1003'))
    ['AF1001', 'AF1002', 'AF1003']
    """
    start, step, end = s.start, s.step, s.stop

    #find where the number is
    start_index = last_nondigit_index(start)
    end_index = last_nondigit_index(end)
    prefix = start[:start_index]

    if prefix != end[:end_index]:
        raise TypeError("Range start and end don't have same prefix")

    if not step:
        step = 1
    range_start = long(start[start_index:])
    range_end = long(end[end_index:])
    field_width = str(len(start) - start_index)
    format_string = '%'+field_width+'.'+field_width+'d'
    return [prefix + format_string % i \
        for i in range(range_start, range_end+1, step)]

def make_lists_of_expanded_slices_of_set_size(s, size_limit=200):
    """Returns a list of Accessions terms from `expand_slice`

    GenBank URLs are limited in size. This helps break up larger lists
    of Accessions (e.g. thousands) into GenBank friendly sizes for down
    stream fetching.

    Parameters
    ----------
    s : list
        Expanded list of accession strings
    size_limit : int
        max items each list should contain

    Returns
    -------
    list
        list of strings with accesions grouped in sets of `size_limit`.

    Examples
    --------
    >>> from skbio.db.base import make_lists_of_expanded_slices_of_set_size
    >>> a = make_lists_of_expanded_slices_of_set_size
    >>> a(slice('HM780503', 'HM780506'), size_limit=3)
    ['HM780503 HM780504 HM780505', 'HM780506']
    """
    full_list = expand_slice(s)
    ls = len(full_list)
    l = []

    # cast to int to be able to use range
    for i in range(int(ls/size_limit+1)):
        start = i * size_limit
        end = (i+1) * size_limit
        subset = full_list[start:end]
        l.append(' '.join(subset))
    return l

def make_lists_of_accessions_of_set_size(s,size_limit=200):
    """Returns list of search terms  containing accessions up to `size_limit`

    This is to help make friendly GenBank urls for fetching large lists 
    of accessions (1000s).

    Parameters
    ----------
    s : list
        List of accessions.
    size_limit : int
        Maximum number of items each list should contain.

    Returns
    -------
    list
        List of strings with the accessions.

    Examples
    --------
    >>> from skbio.db.base import make_lists_of_accessions_of_set_size
    >>> make_lists_of_accessions_of_set_size(['HM780503', 'HM780506',
    ...     'HM780660', 'HM780780'], size_limit=3)
    ['HM780503 HM780506 HM780660', 'HM780780']

    """
    ls = len(s)
    l = []
    for i in range(int(ls/size_limit+1)):
        start = i * size_limit
        end = (i+1) * size_limit
        subset = s[start:end]
        l.append(' '.join(subset))
    return l


def last_nondigit_index(s):
    """Returns the index of s such that `s[i:]` is numeric, or None

    Parameters
    ----------
    s : str
        From where the last non-digit index is retrieved.

    Returns
    -------
    int
        Index of the last non-digit character, will return None if there werer
        no trailing digits.

    Examples
    --------
    >>> from skbio.db.base import last_nondigit_index
    >>> s = 'foo1234'
    >>> idx = last_nondigit_index(s)
    >>> print(s[idx:])
    1234
    """
    for i in range(len(s)):
        if s[i:].isdigit():
            return i
    #if we get here, there weren't any trailing digits
    return None
