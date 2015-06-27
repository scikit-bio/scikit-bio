# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import hashlib
from os import remove, makedirs
from os.path import exists, isdir
from functools import partial
from warnings import warn
from types import FunctionType

from ._decorator import experimental, deprecated


class MiniRegistry(dict):
    def __call__(self, name):
        """Act as a decorator to register functions with self"""
        def decorator(func):
            self[name] = func
            return func
        return decorator

    def copy(self):
        """Useful for inheritance"""
        return self.__class__(super(MiniRegistry, self).copy())

    def formatted_listing(self):
        """Produce an RST list with descriptions."""
        if len(self) == 0:
            return "\tNone"
        else:
            return "\n".join(["\t%r\n\t  %s" %
                             (name, self[name].__doc__.split("\n")[0])
                              for name in sorted(self)])

    def interpolate(self, obj, name):
        """Inject the formatted listing in the second blank line of `name`."""
        # Py2/3 compatible way of calling getattr(obj, name).__func__
        f = getattr(obj, name).__get__(None, type(None))

        if hasattr(f, 'func_code'):
            f2 = FunctionType(f.func_code, f.func_globals, name=f.func_name,
                              argdefs=f.func_defaults, closure=f.func_closure)
        else:
            f2 = FunctionType(f.__code__, f.__globals__, name=f.__name__,
                              argdefs=f.__defaults__, closure=f.__closure__)
        # Conveniently the original docstring is on f2, not the new ones if
        # inheritence is happening. I have no idea why.
        t = f2.__doc__.split("\n\n")
        t.insert(2, self.formatted_listing())
        f2.__doc__ = "\n\n".join(t)

        setattr(obj, name, f2)


def chunk_str(s, n, char):
    """Insert `char` character every `n` characters in string `s`.

    Canonically pronounced "chunkster".

    """
    # Modified from http://stackoverflow.com/a/312464/3776794
    if n < 1:
        raise ValueError(
            "Cannot split string into chunks with n=%d. n must be >= 1." % n)
    return char.join((s[i:i+n] for i in range(0, len(s), n)))


@experimental(as_of="0.4.0")
def cardinal_to_ordinal(n):
    """Return ordinal string version of cardinal int `n`.

    Parameters
    ----------
    n : int
        Cardinal to convert to ordinal. Must be >= 0.

    Returns
    -------
    str
        Ordinal version of cardinal `n`.

    Raises
    ------
    ValueError
        If `n` is less than 0.

    Notes
    -----
    This function can be useful when writing human-readable error messages.

    Examples
    --------
    >>> from skbio.util import cardinal_to_ordinal
    >>> cardinal_to_ordinal(0)
    '0th'
    >>> cardinal_to_ordinal(1)
    '1st'
    >>> cardinal_to_ordinal(2)
    '2nd'
    >>> cardinal_to_ordinal(3)
    '3rd'

    """
    # Taken and modified from http://stackoverflow.com/a/20007730/3776794
    # Originally from http://codegolf.stackexchange.com/a/4712 by Gareth
    if n < 0:
        raise ValueError("Cannot convert negative integer %d to ordinal "
                         "string." % n)
    return "%d%s" % (n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


@experimental(as_of="0.4.0")
def is_casava_v180_or_later(header_line):
    """Check if the header looks like it is Illumina software post-casava v1.8

    Parameters
    ----------
    header_line : bytes
        A header line

    Returns
    -------
    bool
        ``True`` for if casava v1.8+, otherwise ``False``

    Examples
    --------
    >>> from skbio.util import is_casava_v180_or_later
    >>> print(is_casava_v180_or_later('@foo'))
    False
    >>> id_ = '@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0'
    >>> print(is_casava_v180_or_later(id_))
    True

    """
    if not header_line.startswith(b'@'):
        raise ValueError("Non-header line passed in!")
    fields = header_line.split(b':')

    return len(fields) == 10 and fields[7] in b'YN'


@experimental(as_of="0.4.0")
def safe_md5(open_file, block_size=2 ** 20):
    """Computes an md5 sum without loading the file into memory

    Parameters
    ----------
    open_file : file object
        open file handle to the archive to compute the checksum. It
        must be open as a binary file
    block_size : int, optional
        size of the block taken per iteration

    Returns
    -------
    md5 : md5 object from the hashlib module
        object with the loaded file

    Notes
    -----
    This method is based on the answers given in:
    http://stackoverflow.com/a/1131255/379593

    Examples
    --------
    >>> from StringIO import StringIO
    >>> from skbio.util import safe_md5
    >>> fd = StringIO("foo bar baz") # open file like object
    >>> x = safe_md5(fd)
    >>> x.hexdigest()
    'ab07acbb1e496801937adfa772424bf7'
    >>> fd.close()

    """
    md5 = hashlib.md5()
    data = True
    while data:
        data = open_file.read(block_size)
        if data:
            md5.update(data)
    return md5


@experimental(as_of="0.4.0")
def remove_files(list_of_filepaths, error_on_missing=True):
    """Remove list of filepaths, optionally raising an error if any are missing

    Parameters
    ----------
    list_of_filepaths : list of strings
        list with filepaths to remove
    error_on_missing : bool, optional
        whether or not the function should raise an ``OSError`` if a file is
        not found

    Raises
    ------
    OSError
        If a filepath in the list does not exist

    Examples
    --------
    >>> from tempfile import NamedTemporaryFile
    >>> from os.path import exists
    >>> from skbio.util import remove_files
    >>> h = NamedTemporaryFile(delete=False)
    >>> exists(h.name) # it exists
    True
    >>> remove_files([h.name])
    >>> exists(h.name) # and now it's gone
    False

    """
    missing = []
    for fp in list_of_filepaths:
        try:
            remove(fp)
        except OSError:
            missing.append(fp)

    if error_on_missing and missing:
        raise OSError("Some filepaths were not accessible: %s" %
                      '\t'.join(missing))


@experimental(as_of="0.4.0")
def create_dir(dir_name, fail_on_exist=False, handle_errors_externally=False):
    """Create a directory safely and fail meaningfully

    Parameters
    ----------
    dir_name: string
        name of directory to create

    fail_on_exist: bool, optional
        if true raise an error if ``dir_name`` already exists

    handle_errors_externally: bool, optional
        if True do not raise Errors, but return failure codes. This allows to
        handle errors locally and e.g. hint the user at a --force_overwrite
        options.

    Returns
    -------
    return_value : int
        These values are only returned if no error is raised:

        - ``0``:  directory was safely created
        - ``1``:  directory already existed
        - ``2``:  a file with the same name exists
        - ``3``:  any other unspecified ``OSError``

    Notes
    -----
    Depending  of how thorough we want to be we could add tests, e.g. for
    testing actual write permission in an existing dir.

    Examples
    --------
    >>> from skbio.util import create_dir
    >>> from os.path import exists, join
    >>> from tempfile import gettempdir
    >>> from os import rmdir
    >>> new_dir = join(gettempdir(), 'scikitbio')
    >>> create_dir(new_dir)
    0
    >>> exists(new_dir)
    True
    >>> rmdir(new_dir)

    """
    error_code_lookup = _get_create_dir_error_codes()
    # pre-instanciate function with
    ror = partial(_handle_error_codes, dir_name, handle_errors_externally)

    if exists(dir_name):
        if isdir(dir_name):
            # dir is there
            if fail_on_exist:
                return ror(error_code_lookup['DIR_EXISTS'])
            else:
                return error_code_lookup['DIR_EXISTS']
        else:
            # must be file with same name
            return ror(error_code_lookup['FILE_EXISTS'])
    else:
        # no dir there, try making it
        try:
            makedirs(dir_name)
        except OSError:
            return ror(error_code_lookup['OTHER_OS_ERROR'])

    return error_code_lookup['NO_ERROR']


@experimental(as_of="0.4.0")
def find_duplicates(iterable):
    """Find duplicate elements in an iterable.

    Parameters
    ----------
    iterable : iterable
        Iterable to be searched for duplicates (i.e., elements that are
        repeated).

    Returns
    -------
    set
        Repeated elements in `iterable`.

    """
    # modified from qiita.qiita_db.util.find_repeated
    # https://github.com/biocore/qiita
    # see licenses/qiita.txt
    seen, repeated = set(), set()
    for e in iterable:
        if e in seen:
            repeated.add(e)
        else:
            seen.add(e)
    return repeated

flatten_deprecation_reason = (
    "Solutions to this problem exist in the python standarnd library. "
    "Please refer to the following links for good alternatives:\n"
    "http://stackoverflow.com/a/952952/3639023\n"
    "http://stackoverflow.com/a/406199/3639023")


@deprecated(as_of="0.2.3-dev", until="0.4.1",
            reason=flatten_deprecation_reason)
def flatten(items):
    """Removes one level of nesting from items

    Parameters
    ----------
    items : iterable
        list of items to flatten one level

    Returns
    -------
    flattened_items : list
        list of flattened items, items can be any sequence, but flatten always
        returns a list.

    Examples
    --------
    >>> from skbio.util import flatten
    >>> h = [['a', 'b', 'c', 'd'], [1, 2, 3, 4, 5], ['x', 'y'], ['foo']]
    >>> print(flatten(h))
    ['a', 'b', 'c', 'd', 1, 2, 3, 4, 5, 'x', 'y', 'foo']

    """
    warn("skbio.util.flatten is deprecated. Please refer to the following "
         "links for solutions from the standard python library: "
         "http://stackoverflow.com/a/952952/3639023 "
         "http://stackoverflow.com/a/406199/3639023", DeprecationWarning)

    result = []
    for i in items:
        try:
            result.extend(i)
        except TypeError:
            result.append(i)
    return result


def _get_create_dir_error_codes():
    return {'NO_ERROR':      0,
            'DIR_EXISTS':    1,
            'FILE_EXISTS':   2,
            'OTHER_OS_ERROR': 3}


def _handle_error_codes(dir_name, suppress_errors=False,
                        error_code=None):
    """Wrapper function for error_handling.

    dir_name: name of directory that raised the error
    suppress_errors: if True raise Errors, otherwise return error_code
    error_code: the code for the error

    """
    error_code_lookup = _get_create_dir_error_codes()

    if error_code is None:
        error_code = error_code_lookup['NO_ERROR']

    error_strings = \
        {error_code_lookup['DIR_EXISTS']:
         "Directory already exists: %s" % dir_name,
         error_code_lookup['FILE_EXISTS']:
         "File with same name exists: %s" % dir_name,
         error_code_lookup['OTHER_OS_ERROR']:
         "Could not create output directory: %s. " % dir_name +
         "Check the permissions."}

    if error_code == error_code_lookup['NO_ERROR']:
        return error_code_lookup['NO_ERROR']
    if suppress_errors:
        return error_code
    else:
        raise OSError(error_strings[error_code])
