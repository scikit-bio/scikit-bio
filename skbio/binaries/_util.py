# ----------------------------------------------------------------------------
# Copyright (c) 2025--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import ctypes

# ====================================================

# Internal, do not use directly
#
# Global values used for caching purposes
#

_skbb_first_try = True
# Set to invalid values
_skbb_dll = None
_skbb_version = 0

# ====================================================

# Load scikit-bio-binaries shared library object, if exists
# Returns None if it cannot be found
# Note: Should not be used directly, use get_dll instead
def _get_new_skbb_dll():
    import os
    try:
        dll = ctypes.CDLL("libskbb.so")
        if os.environ.get('SKBB_CPU_INFO', 'N') in ('Y','YES'):
            print("INFO (scikit-bio): Using shared library libskbb.so")
    except OSError:
        dll = None
    return dll


# Return scikit-bio-binaries shared library object, if exists
# Returns None if it cannot be found
# Note: Uses caching to minimize overheads
def get_dll():
    global _skbb_first_try
    global _skbb_dll
    if _skbb_first_try:
        _skbb_dll = _get_new_skbb_dll()
        _skbb_first_try = False
    return _skbb_dll

# ====================================================

# Is the scikit-bio-binaries shared library available
# Note that nother skbb function will work, if it does not
def available():
    return get_dll() is not None

# What API version does the scikit-bio-binaries shared library implement
# Returns 0, if the shared library is not available
# Note: Uses caching to minimize overheads
def get_api_version():
    global _skbb_version
    if _skbb_version==0:
        dll = get_dll()
        if dll is not None:
            _skbb_version = dll.skbb_get_api_version()
    return _skbb_version

# ====================================================

def set_random_seed(new_seed):
    dll = get_dll()
    if dll is not None:
        dll.skbb_set_random_seed(ctypes.c_uint(new_seed));
    # since there is no colateral impact
    # just do nothing if the shared library does not exist

