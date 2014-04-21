#! /usr/bin/env python
from __future__ import print_function, absolute_import

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import namedtuple
from itertools import izip
from os.path import exists

import numpy as np

from skbio.core.exception import FileFormatError


class OrdinationResults(namedtuple('OrdinationResults',
                                   ('eigvals', 'species', 'site', 'biplot',
                                    'site_constraints', 'proportion_explained',
                                    'species_ids', 'site_ids'))):
    # To avoid creating a dict, as a namedtuple doesn't have it:
    __slots__ = ()

    def __new__(cls, eigvals, species=None, site=None, biplot=None,
                site_constraints=None, proportion_explained=None,
                species_ids=None, site_ids=None):
        return super(OrdinationResults, cls).__new__(cls, eigvals, species,
                                                     site, biplot,
                                                     site_constraints,
                                                     proportion_explained,
                                                     species_ids, site_ids)

    @classmethod
    def from_file(cls, ord_res_f):
        """Load ordination results from text file.

        Creates a `OrdinationResults` instance from serialized results
        stored as text.

        `ord_res_f` must be a file-like object containing text.

        The ord_res_f format should look like::

            Eigvals<tab>2
            0.096<tab>0.040

            Proportion explained<tab>2
            0.512<tab>0.488

            Species<tab>3<tab>2
            Species1<tab>0.408<tab>0.069
            Species2<tab>-0.115<tab>-0.299
            Species3<tab>-0.309<tab>0.187

            Site<tab>3<tab>2
            Site1<tab>-0.848<tab>0.882
            Site2<tab>-0.220<tab>-1.344
            Site3<tab>1.666<tab>0.470

            Biplot<tab>4<tab>3
            0.422<tab>-0.559<tab>-0.713
            0.988<tab>0.150<tab>-0.011
            -0.556<tab>0.817<tab>0.147
            -0.404<tab>-0.905<tab>-0.127

            Site constraints<tab>3<tab>2
            Site1<tab>-0.848<tab>0.882
            Site2<tab>-0.220<tab>-1.344
            Site3<tab>1.666<tab>0.470

        If a given result attribute is not present (e.g. Biplot), it should be
        still defined and declare its size as 0::

            Biplot<tab>0<tab>0

        Parameters
        ----------
        ord_res_f : iterable of str or str
            Iterable of strings (e.g., open file handle, file-like object, list
            of strings, etc.) or a file path (a string) containing the
            serialized ordination results.

        Returns
        -------
        OrdinationResults
            Instance of type `cls` containing the parsed contents of
            `ord_res_f`.

        """
        # Currently we support either a file or a filepath.
        # This will change once we have a centralized function that
        # takes care of this.
        # Adapted from skbio.core.distance.DissimilarityMatrix.from_file
        if isinstance(ord_res_f, str) and exists(ord_res_f):
            # Check if it's a valid path, if so read the contents
            fd = open(ord_res_f, 'U')
            orf = fd.readlines()
            fd.close()
        elif hasattr(ord_res_f, 'readlines'):
            orf = ord_res_f.readlines()
        else:
            orf = ord_res_f

        # Starting at line 0, we should find the eigvals
        eigvals, curr_line = cls._parse_eigvals(orf)
        # The next line should be an empty line
        curr_line = cls._check_empty_line(orf, curr_line)
        # Now we should find the proportion explained section
        prop_expl, curr_line = cls._parse_proportion_explained(orf, curr_line,
                                                               len(eigvals))
        if len(prop_expl) != len(eigvals):
            raise ValueError('There should be as many proportion explained '
                             'values as eigvals: %d != %d' % (len(prop_expl),
                                                              len(eigvals)))
        # The next line should be an empty line
        curr_line = cls._check_empty_line(orf, curr_line)
        # Next section should be the species section
        species, species_ids, curr_line = cls._parse_coords(orf, curr_line,
                                                            'Species')
        # The next line should be an empty line
        curr_line = cls._check_empty_line(orf, curr_line)
        # Next section should be the site section
        site, site_ids, curr_line = cls._parse_coords(orf, curr_line, 'Site')
        # The next line should be an empty line
        curr_line = cls._check_empty_line(orf, curr_line)
        # Next section should be the biplot section
        biplot, curr_line = cls._parse_biplot(orf, curr_line)
        # The next line should be an empty line
        curr_line = cls._check_empty_line(orf, curr_line)
        # Next section should be the site constraints section
        cons, cons_ids, curr_line = cls._parse_coords(orf, curr_line,
                                                      'Site constraints')
        # The last line should be an empty line
        curr_line = cls._check_empty_line(orf, curr_line)

        return cls(eigvals=eigvals, species=species, site=site, biplot=biplot,
                   site_constraints=cons, proportion_explained=prop_expl,
                   species_ids=species_ids, site_ids=site_ids)

    @staticmethod
    def _parse_eigvals(lines):
        curr_line = 0
        # The first line should contain the Eigvals header:
        # Eigvals<tab>NumEigvals
        header = lines[curr_line].strip().split('\t')
        if len(header) != 2 or header[0] != 'Eigvals':
            raise FileFormatError('Eigvals header not found')

        # Parse how many eigvals are we waiting for
        num_eigvals = int(header[1])
        if num_eigvals == 0:
            raise ValueError('At least one eigval should be present')

        # Parse the eigvals, present on the next line
        # Eigval_1<tab>Eigval_2<tab>Eigval_3<tab>...
        curr_line += 1
        eigvals = np.asarray(lines[curr_line], dtype=np.float64)
        if len(eigvals) != num_eigvals:
            raise ValueError('Expected %d eigvals, but found %d.' %
                             (num_eigvals, len(eigvals)))

        # Update the line cunter to left it after the eigvals section
        return eigvals, curr_line + 1

    @staticmethod
    def _check_empty_line(lines, curr_line):
        if lines[curr_line].strip():
            raise FileFormatError('Expected an empty line')
        # Update the line cunter to left it after the empty line
        return curr_line + 1

    @staticmethod
    def _parse_proportion_explained(lines, curr_line):
        # Parse the proportion explained header:
        # Proportion explained<tab>NumPropExpl
        header = lines[curr_line].strip().split('\t')
        if (len(header) != 2 or
                header[0] != 'Proportion explained'):
            raise FileFormatError('Proportion explained header not found')

        # Parse how many prop expl values are we waiting for
        num_prop_expl = int(header[1])
        if num_prop_expl == 0:
            # The ordination method didn't generate the prop explained vector
            # set it to None
            prop_expl = None
        else:
            # Parse the line with the proportion explained values
            curr_line += 1
            prop_expl = np.asarray(lines[curr_line], dtype=np.float64)
            if len(prop_expl) != num_prop_expl:
                raise ValueError('Expected %d proportion explained values, but'
                                 ' found %d.' % (num_prop_expl,
                                                 len(prop_expl)))
        # Update the line cunter to left it after the prop expl section
        return prop_expl, curr_line + 1

    @staticmethod
    def _parse_coords(lines, curr_line, header):
        # Parse the coords header
        header = lines[curr_line].strip().split('\t')
        if len(header) != 3 or header[0] != header:
            raise FileFormatError('%s header not found.' % header)

        # Parse the dimensions of the coord matrix
        rows = int(header[1])
        cols = int(header[2])

        if rows == 0 and cols == 0:
            # The ordination method didn't generate the coords for 'header'
            # Set the results to None
            coords = None
            ids = None
        elif (rows == 0 and cols != 0) or (rows != 0 and cols == 0):
            # Both dimensions should be 0 or none of them are zero
            raise ValueError('One dimension of %s is 0: %d x %d' %
                             (header, rows, cols))
        else:
            # Parse the coord lines
            coords = np.empty((rows, cols), dtype=np.float64)
            ids = []
            for i in range(rows):
                # Parse the next row of data
                curr_line += 1
                vals = lines[curr_line].strip().split('\t')
                # The +1 comes from the row header (which contains the row id)
                if len(vals) != cols + 1:
                    raise ValueError('Expected %d values, but found %d in row '
                                     '%d.' % (cols, len(vals) - 1, i))
                ids.append(vals[0])
                coords[i, :] = np.asarray(vals[1:], dtype=np.float64)
        # Update the line cunter to left it after the coords section
        return coords, ids, curr_line + 1

    @staticmethod
    def _parse_biplot(lines, curr_line):
        # Parse the biplot header
        header = lines[curr_line].strip().split('\t')
        if len(header) != 3 or header[0] != 'Biplot':
            raise FileFormatError('Biplot header not found.')

        # Parse the dimensions of the Biplot matrix
        rows = int(header[1])
        cols = int(header[2])

        if rows == 0 and cols == 0:
            # The ordination method didn't generate the biplot matrix
            # Set the results to None
            biplot = None
        elif (rows == 0 and cols != 0) or (rows != 0 and cols == 0):
            # Both dimensions should be 0 or none of them are zero
            raise ValueError('One dimension of %s is 0: %d x %d' %
                             (header, rows, cols))
        else:
            # Parse the biplot matrix
            biplot = np.empty((rows, cols), dtype=np.float64)
            for i in range(rows):
                # Parse the next row of data
                curr_line += 1
                vals = lines[curr_line].strip().split('\t')
                if len(vals) != cols:
                    raise ValueError('Expected %d values, but founf %d in row '
                                     '%d.' % (cols, len(vals), i))
                biplot[i, :] = np.asarray(vals, dtype=np.float64)
        # Update the line cunter to left it after the coords section
        return biplot, curr_line + 1

    def to_file(self, out_f):
        """Save the ordination results to file in text format.

        See Also
        --------
        from_file

        Parameters
        ----------
        out_f : file-like object
            File-like object to write serialized data to. Must have a ``write``
            method. It is the caller's responsibility to close `out_f` when
            done (if necessary).
        """

        # Write eigvals
        out_f.write("Eigvals\t%d\n" % self.eigvals.shape)
        out_f.write("%s\n\n" % '\t'.join(np.asarray(self.eigvals,
                    dtype=np.str)))

        # Write proportion explained
        if self.proportion_explained is None:
            out_f.write("Proportion explained\t0\n\n")
        else:
            out_f.write("Proportion explained\t%d\n" %
                        self.proportion_explained.shape)
            out_f.write("%s\n\n" % '\t'.join(np.asarray(
                        self.proportion_explained, dtype=np.str)))

        # Write species
        if self.species is None:
            out_f.write("Species\t0\t0\n\n")
        else:
            out_f.write("Species\t%d\t%d\n" % self.species.shape)
            for id_, vals in izip(self.species_ids, self.species):
                out_f.write("%s\t%s\n" % (id_, '\t'.join(np.asarray(vals,
                            dtype=np.str))))
            out_f.write("\n")

        # Write site
        if self.site is None:
            out_f.write("Site\t0\t0\n\n")
        else:
            out_f.write("Site\t%d\t%d\n" % self.site.shape)
            for id_, vals in izip(self.site_ids, self.site):
                out_f.write("%s\t%s\n" % (id_, '\t'.join(np.asarray(vals,
                            dtype=np.str))))
            out_f.write("\n")

        # Write biplot
        if self.biplot is None:
            out_f.write("Biplot\t0\t0\n\n")
        else:
            out_f.write("Biplot\t%d\t%d\n" % self.biplot.shape)
            for vals in self.biplot:
                out_f.write("%s\n" % '\t'.join(np.asarray(vals, dtype=np.str)))
            out_f.write("\n")

        # Write site-constraints
        if self.site_constraints is None:
            out_f.write("Site constraints\t0\t0\n")
        else:
            out_f.write("Site constraints\t%d\t%d\n" %
                        self.site_constraints.shape)
            for id_, vals in izip(self.site_ids, self.site_constraints):
                out_f.write("%s\t%s\n" % (id_, '\t'.join(np.asarray(vals,
                            dtype=np.str))))


class Ordination(object):
    short_method_name = 'Overwrite in subclass!'
    long_method_name = 'Overwrite in subclass!'
