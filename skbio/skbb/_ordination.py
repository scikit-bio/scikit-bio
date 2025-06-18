# ----------------------------------------------------------------------------
# Copyright (c) 2025--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import ctypes
import numpy as np
from skbio.skbb._util import (skbb_get_api_version, get_skbb_dll, skbb_set_random_seed)

# ====================================================

# perform PCoA using the FSVD method
def skbb_pcoa_fsvd_inplace(
    distance_matrix,
    number_of_dimensions,
    seed=None,
):
    if skbb_get_api_version()>=1: # minimum version that support pcoa
        if (seed is not None):
            # TBD: Switch from global lo local
            skbb_set_random_seed(seed)
        # get internal representations
        distance_matrix_shape0 = distance_matrix.shape[0]
        distance_matrix_data = distance_matrix.data
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
        if distance_matrix_data.dtype == np.dtype('float64'):
            dll.skbb_pcoa_fsvd_inplace_fp64(
                ctypes.c_uint(distance_matrix_shape0),
                distance_matrix_data.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_uint(number_of_dimensions),
                eigenvalues.ctypes.data_as(ctypes.c_void_p),
                samples.ctypes.data_as(ctypes.c_void_p),
                proportion_explained.ctypes.data_as(ctypes.c_void_p) )
        elif distance_matrix_data.dtype == np.dtype('float32'):
            dll.skbb_pcoa_fsvd_inplace_fp64(
                ctypes.c_uint(distance_matrix_shape0),
                distance_matrix_data.ctypes.data_as(ctypes.c_void_p),
                ctypes.c_uint(number_of_dimensions),
                eigenvalues.ctypes.data_as(ctypes.c_void_p),
                samples.ctypes.data_as(ctypes.c_void_p),
                proportion_explained.ctypes.data_as(ctypes.c_void_p) )
        else:
            raise TypeError("distance_matrix must be either float32 or float64")
        # if we got here, everything went well
        return eigenvalues, samples, proportion_explained
    else:
        raise ImportError("skbb_pcoa_fsvd_inplace not avaialble")

