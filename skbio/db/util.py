r"""Retrieve information from web databases.
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


class UrlGetter(object):
    """Base class to submit queries with EUtils"""
    defaults = {}       #override in derived classes -- default values
    printed_fields = {}  #override in derived classes -- fields to print
    base_url = ''        #override in derived classes
    key_value_delimiter = '='
    field_delimiter = '&'

    def __init__(self, **kwargs):
        """Returns new instance with arbitrary kwargs."""
        self.__dict__.update(self.defaults)
        self.__dict__.update(kwargs)
        self._temp_args = {}

    def __str__(self):
        to_get = self.__dict__.copy()
        to_get.update(self._temp_args)
        fields = []

        for key, value in to_get.items():
            if key in self.printed_fields:
                fields.append(quote_plus(key)+self.key_value_delimiter+
                              quote_plus(str(value)))

        return self.base_url + self.field_delimiter.join(fields)

    def open(self, **kwargs):
        """Returns a stream handle to URL result, temporarily overriding kwargs.."""
        self._temp_args = kwargs
        result = urlopen(str(self))
        self._temp_args = {}
        return result

    def read(self, **kwargs):
        """Gets URL and reads into memory, temporarily overriding kwargs."""
        result = self.open(**kwargs)
        data = result.read()
        result.close()
        return data

    def retrieve(self, fname, **kwargs):
        """Gets URL and writes to file fname, temporarily overriding kwargs. 
        
        Note: produces no return value."""
        self._temp_args = kwargs
        urlretrieve(str(self), fname)
        self._temp_args = None
        
def expand_slice(s):
    """Takes a start and end accession, and gets the whole range.

    WARNING: Unlike standard slices, includes the last item in the range.
    In other words, obj[AF1001:AF1010] will include AF1010.
    
    Both accessions must have the same non-numeric prefix.
    """
    start, step, end = s.start, s.step, s.stop
    #find where the number is
    start_index = last_nondigit_index(start)
    end_index = last_nondigit_index(end)
    prefix = start[:start_index]
    if prefix != end[:end_index]:
        raise TypeError, "Range start and end don't have same prefix"
    if not step:
        step = 1
    range_start = long(start[start_index:])
    range_end = long(end[end_index:])
    field_width = str(len(start) - start_index)
    format_string = '%'+field_width+'.'+field_width+'d'
    return [prefix + format_string % i \
        for i in range(range_start, range_end+1, step)]

def make_lists_of_expanded_slices_of_set_size(s,size_limit=200):
    """Returns a list of Accessions terms from 'expand_slice'.
    GenBank URLs are limited in size. This helps break up larger lists
    of Accessions (e.g. thousands) into GenBank friendly sizes for down
    stream fetching.
       -s : slice of accessions
       -size_limit : max items each list should contain
    """
    full_list = expand_slice(s)
    ls = len(full_list)
    l = []
    for i in range(ls/size_limit+1):
        start = i * size_limit
        end = (i+1) * size_limit
        subset = full_list[start:end]
        l.append(' '.join(subset))
    return l

def make_lists_of_accessions_of_set_size(s,size_limit=200):
    """Returns list of search terms  that contain accessions up to the size
    'size_limit'
    This is to help make friendly GenBank urls for fetching large lists 
    of accessions (1000s).
        -s : list of accessions
        -size_limit : max items each list should contain
    """
    ls = len(s)
    l = []
    for i in range(ls/size_limit+1):
        start = i * size_limit
        end = (i+1) * size_limit
        subset = s[start:end]
        l.append(' '.join(subset))
    return l


def last_nondigit_index(s):
    """Returns the index of s such that s[i:] is numeric, or None."""
    for i in range(len(s)):
        if s[i:].isdigit():
            return i
    #if we get here, there weren't any trailing digits
    return None
