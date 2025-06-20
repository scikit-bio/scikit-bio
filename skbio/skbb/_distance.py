# ----------------------------------------------------------------------------
# Copyright (c) 2025--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import ctypes
import numpy as np
from numbers import Integral
from skbio.skbb._util import (skbb_get_api_version, get_skbb_dll)

# ====================================================

# Check if skbb_pcoa_fsvd is available
# same inputs as the full function, in case we support only a subset
# If it returns True, it is safe to call skbb_pcoa_fsvd
def skbb_permanova_available(
    distance_matrix,
    grouping,
    permutations,
    seed=None,
):
    if skbb_get_api_version()>=1: # minimum version that support permanova, includes check for library existence
        # check it is a positive number
        if not isinstance(permutations, Integral):
            return False
        elif permutations<1:
            return False
        if seed is not None: # None is OK
            # check it is a non-negative number
            if not isinstance(seed, Integral):
                return False
            elif (seed<0):
                return False
        return True
    else:
        return False


# perform permanova
# Note: Seed must be either a non-negative integer or None
def skbb_permanova(
    distance_matrix,
    grouping,
    permutations,
    seed=None,
):
    if skbb_get_api_version()>=1: # minimum version that support pcoa
        if permutations<1:
            raise ValueError("permutations must be a positive number")
        if (seed is not None):
            # just check it actually is a non-negative number
            if not isinstance(seed, Integral):
                raise TypeError("seed must be an integer")
            elif (seed>=0):
                int_seed = seed
            else:
                raise ValueError("seed must be a non-negative number")
        else:
            # the skbb API expects a negative number, when seed is invalid
            int_seed = -1
        if isinstance(distance_matrix,np.ndarray):
            # already a raw matrix, just use
            distance_matrix_data = distance_matrix
        else: 
            # we are assuming it is a DistanceMatrix object
            # get internal representations
            distance_matrix_data = distance_matrix.data
        distance_matrix_shape0 = distance_matrix_data.shape[0]
        if (len(distance_matrix_data.shape)!=2 or 
            distance_matrix_data.shape[1] != distance_matrix_shape0):
               raise TypeError("distance_matrix not square")
        # skbb expects uint32 array, so convert if needed (cheap)
        if grouping.dtype == np.dtype('uint32'):
            grouping_data = grouping
        else:
            grouping_data = grouping.astype(np.uint32)
        # ready to call the C functions
        dll = get_skbb_dll()
        i_mdim = ctypes.c_uint(distance_matrix_shape0)
        i_mat = distance_matrix_data.ctypes.data_as(ctypes.c_void_p)
        i_grp = grouping_data.ctypes.data_as(ctypes.c_void_p)
        i_n_perm = ctypes.c_uint(permutations)
        i_seed = ctypes.c_int(int_seed)
        if distance_matrix_data.dtype == np.dtype('float64'):
            o_fstat = ctypes.c_double()
            o_pvalue = ctypes.c_double()
            dll.skbb_permanova_fp64(
                    i_mdim, i_mat, i_grp, i_n_perm, i_seed,
                    ctypes.byref(o_fstat), ctypes.byref(o_pvalue))
        elif distance_matrix_data.dtype == np.dtype('float32'):
            o_fstat = ctypes.c_float()
            o_pvalue = ctypes.c_float()
            dll.skbb_permanova_fp32(
                    i_mdim, i_mat, i_grp, i_n_perm, i_seed,
                    ctypes.byref(o_fstat), ctypes.byref(o_pvalue))
        else:
            raise TypeError("distance_matrix type must be either float32 or float64")
        # if we got here, everything went well
        return o_fstat.value, o_pvalue.value
    else:
        raise ImportError("skbb_permanova not available")

