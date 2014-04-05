#! /usr/bin/env python
from __future__ import print_function, absolute_import

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import namedtuple


class OrdinationResults(namedtuple('OrdinationResults',
                                   ('eigvals', 'species', 'site', 'biplot',
                                    'site_constraints', 'proportion_explained',
                                    'ids'))):
    __slots__ = ()  # To avoid creating a dict, as a namedtuple
                    # doesn't have it

    def __new__(cls, eigvals, species, site=None, biplot=None,
                site_constraints=None, proportion_explained=None, ids=None):
        return super(OrdinationResults, cls).__new__(cls, eigvals, species,
                                                     site, biplot,
                                                     site_constraints,
                                                     proportion_explained,
                                                     ids)


class Ordination(object):
    short_method_name = 'Overwrite in subclass!'
    long_method_name = 'Overwrite in subclass!'
