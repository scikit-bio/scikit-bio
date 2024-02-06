# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------


def _pprint_strs(
    strs,
    max_chars=80,
    delimiter=", ",
    suffix="...",
):
    """Pretty-print an iterable of strings, truncating if necessary."""
    # Adapted from http://stackoverflow.com/a/250373
    joined_str = delimiter.join(repr(s) for s in strs)

    if len(joined_str) > max_chars:
        truncated = joined_str[: max_chars + 1].split(delimiter)[0:-1]
        joined_str = delimiter.join(truncated)
        if joined_str:
            joined_str += delimiter
        joined_str += suffix

    return joined_str
