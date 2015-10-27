# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import Mapping


class Feature(Mapping):
    '''Store the metadata of a sequence feature.

    Parameters
    ----------
    kwargs :
    '''
    def __init__(self, **kwargs):
        for attr in kwargs:
            setattr(self, attr, kwargs[attr])

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        for attr, value in other:
            if getattr(self, attr) != value:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __iter__(self):
        for attr, value in self.__dict__.iteritems():
            yield attr, value

    def __getitem__(self, attr):
        return getattr(self, attr)

    def __len__(self):
        return len(self.__dict__)

    def __repr__(self):
        return '__'.join('{0}:{1}'.format(attr, value)
                         for attr, value in self)
