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
def skbb_pcoa_fsvd_available(
    distance_matrix,
    number_of_dimensions,
    inplace=False,
    seed=None,
):
    if skbb_get_api_version()>=1: # minimum version that support pcoa, includes check for library existence
        # check it is a positive number
        if not isinstance(number_of_dimensions, Integral):
            return False
        elif number_of_dimensions<1:
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


# perform PCoA using the FSVD method
# Note: Seed must be either a non-negative integer or None
def skbb_pcoa_fsvd(
    distance_matrix,
    number_of_dimensions,
    inplace=False,
    seed=None,
):
    if skbb_get_api_version()>=1: # minimum version that support pcoa
        if not isinstance(number_of_dimensions, Integral):
            raise ValueError("number_of_dimensions must be an integer value")
        if number_of_dimensions<1:
            raise ValueError("number_of_dimensions must be a positive number")
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
            # already a raw matrix, jsut use
            distance_matrix_data = distance_matrix
        else: 
            # we are assuming it is a DistanceMatrix object
            # get internal representations
            distance_matrix_data = distance_matrix.data
        distance_matrix_shape0 = distance_matrix_data.shape[0]
        if (len(distance_matrix_data.shape)!=2 or 
            distance_matrix_data.shape[1] != distance_matrix_shape0):
               raise TypeError("distance_matrix not square")
        if number_of_dimensions>distance_matrix_shape0:
            raise ValueError("number_of_dimensions cannot be larger that matrix size")
        # create output buffers
        eigenvalues = np.ndarray(
                shape=(number_of_dimensions,),
                dtype=distance_matrix_data.dtype )
        proportion_explained = np.ndarray(
                shape=(number_of_dimensions,),
                dtype=distance_matrix_data.dtype )
        samples = np.ndarray(
                shape=(distance_matrix_shape0, number_of_dimensions),
                dtype=distance_matrix_data.dtype ,
                order="C")
        # ready to call the C functions
        dll = get_skbb_dll()
        i_mdim = ctypes.c_uint(distance_matrix_shape0)
        i_mat = distance_matrix_data.ctypes.data_as(ctypes.c_void_p)
        i_n_eigh = ctypes.c_uint(number_of_dimensions)
        i_seed = ctypes.c_int(int_seed)
        o_ev = eigenvalues.ctypes.data_as(ctypes.c_void_p)
        o_sp = samples.ctypes.data_as(ctypes.c_void_p)
        o_pe = proportion_explained.ctypes.data_as(ctypes.c_void_p)
        if distance_matrix_data.dtype == np.dtype('float64'):
            if (inplace):
                dll.skbb_pcoa_fsvd_inplace_fp64(i_mdim, i_mat, i_n_eigh, i_seed, o_ev, o_sp, o_pe)
            else:
                dll.skbb_pcoa_fsvd_fp64(i_mdim, i_mat, i_n_eigh, i_seed, o_ev, o_sp, o_pe)
        elif distance_matrix_data.dtype == np.dtype('float32'):
            if (inplace):
                dll.skbb_pcoa_fsvd_inplace_fp32(i_mdim, i_mat, i_n_eigh, i_seed, o_ev, o_sp, o_pe)
            else:
                dll.skbb_pcoa_fsvd_fp32(i_mdim, i_mat, i_n_eigh, i_seed, o_ev, o_sp, o_pe)
        else:
            raise TypeError("distance_matrix type must be either float32 or float64")
        # if we got here, everything went well
        return eigenvalues, samples, proportion_explained
    else:
        raise ImportError("skbb_pcoa_fsvd not available")

