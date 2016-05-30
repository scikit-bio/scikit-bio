# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import Mapping
from skbio.util._misc import merge_dicts


class Feature(Mapping):
    '''Store the metadata of a sequence feature.

    It is implemented as frozendict and can be used similarly
    as built-in ``dict``.

    Parameters
    ----------
    args : tuple
        Positional arguments that can be passed to ``dict``
    kwargs : dict
        Keyword arguments of feature name and feature value, which can
        be passed to ``dict``.
    '''
    def __init__(self, *args, **kwargs):
        self.__d = dict(*args, **kwargs)
        self._hash = None
        # make sure the values in the dict are also hashable/immutable
        for k in self.__d:
            hash(self.__d[k])

    def __len__(self):
        return len(self.__d)

    def __getitem__(self, key):
        return self.__d[key]

    def __iter__(self):
        return iter(self.__d)

    def __repr__(self):
        return ';'.join('{0}:{1}'.format(k, self[k]) for k in self)

    def __hash__(self):
        if self._hash is None:
            self._hash = hash(frozenset(self.items()))
        return self._hash

    def __lt__(self, other):
        return hash(self) < hash(other)

    def __gt__(self, other):
        return hash(self) > hash(other)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def update(self, *args, **kwargs):
        """
        Creates a new features object.

        Updates the existing attributes of the current Feature object
        and returns the modified Feature object.

        Parameters
        ----------
        args : tuple
            Positional arguments that can be passed to ``dict``
        kwargs : dict
            Keyword arguments of feature name and feature value, which can
            be passed to ``dict``.

        Returns
        -------
        skbio.sequence.Feature
        """
        __d = dict(*args, **kwargs)
        return Feature(**merge_dicts(self.__d, __d))
