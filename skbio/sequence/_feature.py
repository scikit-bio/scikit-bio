# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from collections import Mapping


class Feature(Mapping):
    '''Store the metadata of a sequence feature.

    Parameters
    ----------
    kwargs :
    '''
    def __init__(self, **kwargs):
        self.__d = dict(**kwargs)

    # def __eq__(self, other):
    #     if self.__class__ != other.__class__:
    #         return False
    #     for attr, value in other:
    #         if getattr(self, attr) != value:
    #             return False
    #     return True

    # def __ne__(self, other):
    #     return not (self == other)

    def __iter__(self):
        return iter(self.__d)

    def __len__(self):
        return len(self.__d)

    def __getitem__(self, key):
        return self.__d[key]

    def __repr__(self):
        return ';'.join('{0}:{1}'.format(k, self[k]) for k in self)

    def __hash__(self):
        return hash(frozenset(self.items()))
