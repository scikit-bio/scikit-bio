# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import with_metaclass

from abc import ABCMeta, abstractmethod


class SkbioObject(with_metaclass(ABCMeta, object)):
    """Abstract base class defining core API common to all scikit-bio objects.

    Public scikit-bio classes should subclass this class to ensure a common,
    core API is present. All abstract methods and properties defined here must
    be implemented in subclasses, otherwise they will not be instantiable.

    """
    @abstractmethod
    def __str__(self):
        pass
