# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import struct


def _float_to_int(number: float) -> int:
    # convert float to a native-endian 8-byte sequence (a double):
    # (alignment doesn't matter because this isn't a struct, and since we are
    #  on the same hardware when we go to int, the endian-ness doesn't matter
    #  either)
    bytes_ = struct.pack("=d", number)
    (integer,) = struct.unpack("=Q", bytes_)
    return integer


def _int_to_float(number: int) -> float:
    bytes_ = struct.pack("=Q", number)
    (float_,) = struct.unpack("=d", bytes_)
    return float_


# 1954 is an homage to R's NA value which is a quiet NaN with a mantissa which
# appears to represent 1954, birth-year of Ross Ihaka (source: Hadley Wickham)
# http://www.markvanderloo.eu
#   /yaRb/2012/07/08/representation-of-numerical-nas-in-r-and-the-1954-enigma/
# It also happens to let us tolerate small negative values (such as -1 used
# in pd.Categorical.codes) without trashing the entire NaN.
# ('Q' from struct does fortunately catch this issue before it becomes a
#  larger problem)
_R_OFFSET = 1954
_DEFAULT_NAN_INT = _float_to_int(float("nan"))
# at this point, calling `bin(_DEFAULT_NAN_INT)` should produce a
# 64-bit positive quiet nan:
# 0 11111111111 1000000000000000000000000000000000000000000000000000
# https://www.csee.umbc.edu/courses/undergraduate/CMSC211/spring03
#     /burt/tech_help/IEEE-754references.html
# unless Python changes some implementation detail, which isn't a problem so
# long as XOR is used instead of AND


def make_nan_with_payload(payload: int, namespace: int = 255):
    """Construct a NaN with a namespace and payload.

    The payload must be in the range [-1953, 2141]
    The namespace must be in the range [0, 255]

    sign  exp                            mantissa
    v v---------v v----------------------------------------------------------v
     +qNaN "header" (includes 1 bit of the mantissa)  namespace     payload
    v-------------v                                   v-------v v------------v
    0 11111111111 10000000 00000000 00000000 00000000 0000 0000 0000 0000 0000

    The namespace + payload requires 20 bits of the mantissa, which will
    support both 32-bit floats and 64-bit doubles.

    The purpose is to allow enumerations to be identified and values preserved.
    Custom enumerations will have a namespace of 255 and will require user
    guidance. Other built-in schemes should organize themselves within an
    unsigned byte. The enumeration values are then stored in the payload which
    uses an offset of 1954 to distinguish between default +qNaN and an
    enumeration scheme. This also permits small negative values in the payload.

    """
    # To be safe, we will XOR our payload (instead of AND) so that we can take
    # the difference later, even if the default NaN changes to include a
    # mantissa payload for some reason
    nan_int = _DEFAULT_NAN_INT ^ (namespace << 12)
    nan_int = nan_int ^ (_R_OFFSET + payload)

    return _int_to_float(nan_int)


def get_payload_from_nan(nan: float):
    nan_int = _float_to_int(nan)
    namespaced_payload = nan_int ^ _DEFAULT_NAN_INT
    if namespaced_payload == 0:
        return (None, None)
    namespace = namespaced_payload >> 12
    payload = namespaced_payload - (namespace << 12)

    return (payload - _R_OFFSET, namespace)
