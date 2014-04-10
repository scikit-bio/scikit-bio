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
import numpy as np


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


        The first line
        (header) must contain the IDs of each object. The subsequent lines
        must contain an ID followed by each distance (float) between the
        current object and all other objects, where the order of objects is
        determined by the header line.

        For example, a 2x2 distance matrix with IDs ``'a'`` and ``'b'`` might
        look like::

            <del>a<del>b
            a<del>0.0<del>1.0
            b<del>1.0<del>0.0

        where ``<del>`` is the delimiter between elements.

        Parameters
        ----------
        ord_res_f : iterable of str
            Iterable of strings (e.g., open file handle, file-like object, list
            of strings, etc.) containing the serialized ordination results.

        Returns
        -------
        OrdinationResults
            Instance of type `cls` containing the parsed contents of
            `ord_res_f`.

        Notes
        -----
        Whitespace-only lines can occur anywhere throughout the "file" and are
        ignored. Lines starting with ``#`` are treated as comments and ignored.
        These comments can only occur *before* the ID header.

        IDs will have any leading/trailing whitespace removed when they are
        parsed.

        """
        # We aren't using np.loadtxt because it uses *way* too much memory
        # (e.g, a 2GB matrix eats up 10GB, which then isn't freed after parsing
        # has finished). See:
        # http://mail.scipy.org/pipermail/numpy-tickets/2012-August/006749.html

        # Strategy:
        #     - find the header
        #     - initialize an empty ndarray
        #     - for each row of data in the input file:
        #         - populate the corresponding row in the ndarray with floats

        # We use iter() as we want to take a single pass over the iterable and
        # maintain our current position after finding the header (mainly
        # necessary for something like a list of strings).
        dm_f = iter(dm_f)
        ids = cls._parse_ids(dm_f, delimiter)
        num_ids = len(ids)
        data = np.empty((num_ids, num_ids), dtype='float')

        # curr_row_idx keeps track of the row index within the data matrix.
        # We're not using enumerate() because there may be
        # empty/whitespace-only lines throughout the data matrix. We want to
        # ignore those and only count the actual rows of data.
        curr_row_idx = 0
        for line in dm_f:
            line = line.strip()

            if not line:
                continue
            elif curr_row_idx >= num_ids:
                # We've hit a nonempty line after we already filled the data
                # matrix. Raise an error because we shouldn't ignore extra
                # data.
                raise DistanceMatrixFormatError(
                    "Encountered extra rows without corresponding IDs in the "
                    "header.")

            tokens = line.split(delimiter)

            # -1 because the first element contains the current ID.
            if len(tokens) - 1 != num_ids:
                raise DistanceMatrixFormatError(
                    "There are %d values in row number %d, which is not equal "
                    "to the number of IDs in the header (%d)."
                    % (len(tokens) - 1, curr_row_idx + 1, num_ids))

            curr_id = tokens[0].strip()
            expected_id = ids[curr_row_idx]
            if curr_id == expected_id:
                data[curr_row_idx, :] = np.asarray(tokens[1:], dtype='float')
            else:
                raise IDMismatchError(curr_id, expected_id)

            curr_row_idx += 1

        if curr_row_idx != num_ids:
            raise MissingDataError(curr_row_idx, num_ids)

        return cls(data, ids)

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
            for vals in self.site_constraints:
                out_f.write("%s\t%s\n" % (id_, '\t'.join(np.asarray(vals,
                            dtype=np.str))))


class Ordination(object):
    short_method_name = 'Overwrite in subclass!'
    long_method_name = 'Overwrite in subclass!'
