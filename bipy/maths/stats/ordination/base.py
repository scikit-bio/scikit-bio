#! /usr/bin/env python
from __future__ import print_function, absolute_import
from collections import namedtuple
try:  # py2 compatibility
    from itertools import izip as zip, imap as map
except ImportError:
    pass
from itertools import count
import numpy as np
import matplotlib.pyplot as plt

# Ordination results should maybe only have two compulsory fields:
# eigvals and objects (aka rows, ?species?). Then, optional fields:
# descriptors (aka columns, ?site?), biplot, site_constraints.


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

    def biplot(self, scaling=2, choices=[0, 1]):
        choices = list(choices)

        def scatter_with_names(ax, data, names, marker, color, label):
            """Helper to plot scatter points and their names."""
            xdata, ydata = data[:, choices].T
            ax.scatter(xdata, ydata, marker=marker, color=color, label=label)
            for namei, xi, yi in zip(names, xdata, ydata):
                ax.annotate(namei, (xi, yi), xytext=(0, 5),
                            textcoords='offset points', ha='center')

        def make_names(prefix):
            prefix += "{n}"
            return map(lambda i: prefix.format(n=i), count())

        # Extract scores
        scores = self.scores(scaling)
        eigvals = scores.eigvals
        species_scores = scores.species
        site_scores = scores.site
        biplot_scores = scores.biplot  # Can be None
        # Actual plotting
        fig, ax = plt.subplots()
        plt.axvline(color='k')
        plt.axhline(color='k')
        scatter_with_names(ax, site_scores, make_names("R"), marker='o',
                           color='b', label='Sites (rows)')
        scatter_with_names(ax, species_scores, make_names("C"),
                           marker='s', color='r', label='Species (columns)')
        if biplot_scores is not None:
            for (xi, yi) in biplot_scores[:, choices]:
                ax.arrow(0, 0, xi, yi, head_width=0.1, facecolor='none')

        for i, set_label in enumerate([ax.set_xlabel, ax.set_ylabel]):
            set_label('{ax_name}{ax_n} - {proportion:.3%}'.format(
                ax_name=self.__class__.short_method_name,
                ax_n=choices[i],
                proportion=(eigvals / np.sum(eigvals))[choices[i]]))
        ax.legend(loc='best').draggable()
        ax.set_title('{0} (scaling = {1})'.format(
            self.__class__.long_method_name, scaling))
        ax.set_aspect('equal')
        plt.show()
