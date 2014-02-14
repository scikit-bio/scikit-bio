#! /usr/bin/env python
from __future__ import print_function, absolute_import
from collections import namedtuple


class OrdinationResults(namedtuple('OrdinationResults',
                                   ('eigvals', 'species', 'site', 'biplot',
                                    'site_constraints'))):
    __slots__ = ()  # To avoid creating a dict, as a namedtuple
                    # doesn't have it

    def __new__(cls, eigvals, species, site=None, biplot=None,
                site_constraints=None):
        return super(OrdinationResults, cls).__new__(cls, eigvals, species,
                                                     site, biplot,
                                                     site_constraints)


class Ordination(object):
    short_method_name = 'Overwrite in subclass!'
    long_method_name = 'Overwrite in subclass!'
