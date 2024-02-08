# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def find_duplicates(iterable):
    """Find duplicate values in an iterable.

    Parameters
    ----------
    iterable : iterable
        Iterable to search for duplicates.

    Returns
    -------
    set
        Values that are duplicated in `iterable`.

    Notes
    -----
    Values in `iterable` must be hashable.

    """
    # Modified from https://stackoverflow.com/a/9835819/3776794 to return
    # duplicates instead of remove duplicates from an iterable.
    seen = set()
    duplicates = set()
    for value in iterable:
        if value in seen:
            duplicates.add(value)
        else:
            seen.add(value)
    return duplicates
