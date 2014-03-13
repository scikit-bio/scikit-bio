#!/usr/bin/env python
r"""
Miscellaneous Utilities (:mod:`bipy.util.misc`)
===============================================

.. currentmodule:: bipy.util.misc

This module provides miscellaneous useful utility classes and methods that do
not fit in any specific module.

Functions
---------

.. autosummary::
   :toctree: generated/

   safe_md5
   remove_files

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import hashlib
from os import remove, makedirs
from os.path import exists, isdir


def safe_md5(open_file, block_size=2**20):
    """Computes an md5 sum without loading the file into memory

    Parameters
    ----------
    open_file : file object
        open file handle to the archive to compute the checksum
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
    >>> from bipy.util.misc import safe_md5
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


def curry(f, *a, **kw):
    """curry(f,x)(y) = f(x,y) or =lambda y: f(x,y)
    
    modified from python cookbook"""
    def curried(*more_a, **more_kw):
        return f(*(a + more_a), **dict(kw, **more_kw))

    ## make docstring for curried funtion
    curry_params = []
    if a:
        curry_params.extend([e for e in a])
    if kw:
        curry_params.extend(['%s=%s' % (k, v) for k, v in kw.items()])
    #str it to prevent error in join()
    curry_params = map(str, curry_params)

    try:
        f_name = f.func_name
    except:  #e.g.  itertools.groupby failed .func_name 
        f_name = '?'

    curried.__doc__ = ' curry(%s,%s)\n'\
        '== curried from %s ==\n %s'\
        % (f_name, ', '.join(curry_params), f_name, f.__doc__)

    return curried


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
    >>> from bipy.util.misc import remove_files
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

def create_dir(dir_name, fail_on_exist=False, handle_errors_externally=False):
    """Create a dir safely and fail meaningful.

    dir_name: name of directory to create

    fail_on_exist: if true raise an error if dir already exists
    
    handle_errors_externally: if True do not raise Errors, but return
                   failure codes. This allows to handle errors locally and
                   e.g. hint the user at a --force_overwrite options.
                   
    returns values (if no Error raised):
    
         0:  dir could be safely made
         1:  directory already existed
         2:  a file with the same name exists          
         3:  any other unspecified OSError


    See qiime/denoiser.py for an example of how to use this mechanism.

    Note: Depending  of how thorough we want to be we could add tests,
          e.g. for testing actual write permission in an existing dir.
    """
    error_code_lookup = get_create_dir_error_codes()
    #pre-instanciate function with
    ror = curry(handle_error_codes, dir_name, handle_errors_externally)

    if exists(dir_name):
        if isdir(dir_name):
            #dir is there
            if fail_on_exist:
                return ror(error_code_lookup['DIR_EXISTS'])
            else:
                return error_code_lookup['DIR_EXISTS']
        else:
            #must be file with same name
            return ror(error_code_lookup['FILE_EXISTS'])
    else:
        #no dir there, try making it
        try:
            makedirs(dir_name)
        except OSError:
            return ror(error_code_lookup['OTHER_OS_ERROR'])
    
    return error_code_lookup['NO_ERROR']

#some error codes for creating a dir
def get_create_dir_error_codes():
    return {'NO_ERROR':      0,
            'DIR_EXISTS':    1,
            'FILE_EXISTS':   2,
            'OTHER_OS_ERROR':3}

def handle_error_codes(dir_name, supress_errors=False,
                       error_code=None):
    """Wrapper function for error_handling.

    dir_name: name of directory that raised the error
    suppress_errors: if True raise Errors, otherwise return error_code
    error_code: the code for the error
    """
    error_code_lookup = get_create_dir_error_codes()
    
    if error_code == None:
        error_code = error_code_lookup['NO_ERROR']
    
    error_strings = \
        {error_code_lookup['DIR_EXISTS'] :
          "Directory already exists: %s" % dir_name,
         error_code_lookup['FILE_EXISTS'] : 
          "File with same name exists: %s" % dir_name,
         error_code_lookup['OTHER_OS_ERROR']: 
          "Could not create output directory: %s. " % dir_name +
          "Check the permissions."}

    if error_code == error_code_lookup['NO_ERROR']:
        return error_code_lookup['NO_ERROR']
    if supress_errors:
        return error_code
    else:
        raise OSError, error_strings[error_code]
