r"""
Ordination results format (:mod:`skbio.io.ordres`)
==================================================

.. currentmodule:: skbio.io.ordres

The ordination results file format (``ordres``) stores the results of an
ordination method in a human-readable, text-based format. The format supports
storing the results of various ordination methods available in scikit-bio,
including (but not necessarily limited to) PCoA, CA, RDA, and CCA.

Format Support
--------------
**Has Sniffer: Yes**

+----------+----------+------------------------------------------------------+
|**Reader**|**Writer**|                   **Object Class**                   |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.stats.ordination.OrdinationResults`       |
+----------+----------+------------------------------------------------------+

Format Specification
--------------------
The format is text-based, consisting of six attributes that describe the
ordination results:

- ``Eigvals``: 1-D
- ``Proportion explained``: 1-D
- ``Species``: 2-D
- ``Site``: 2-D
- ``Biplot``: 2-D
- ``Site constraints``: 2-D

The attributes in the file *must* be in this order.

Each attribute is defined in its own section of the file, where sections are
separated by a blank (or whitespace-only) line. Each attribute begins with a
header line, which contains the attribute's name (as listed above), followed by
a tab character, followed by one or more tab-separated dimensions (integers)
that describe the shape of the attribute's data.

The attribute's data follows its header line, and is stored in tab-separated
format. ``Species``, ``Site``, and ``Site constraints`` store species and site
IDs, respectively, as the first column, followed by the 2-D data array.

An example of this file format might look like::

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

If a given result attribute is not present (e.g. ``Biplot``), it should still
be defined and declare its dimensions as 0. For example::

    Biplot<tab>0<tab>0

All attributes are optional except for ``Eigvals``.

Examples
--------
Assume we have the following tab-delimited text file storing the
ordination results in ``ordres`` format::

    Eigvals\t2
    0.0961330159181\t0.0409418140138

    Proportion explained\t0

    Species\t3\t2
    Species1\t0.408869425742\t0.0695518116298
    Species2\t-0.1153860437\t-0.299767683538
    Species3\t-0.309967102571\t0.187391917117

    Site\t3\t2
    Site1\t-0.848956053187\t0.882764759014
    Site2\t-0.220458650578\t-1.34482000302
    Site3\t1.66697179591\t0.470324389808

    Biplot\t0\t0

    Site constraints\t0\t0

Load the ordination results from the file:

>>> from StringIO import StringIO
>>> from skbio.stats.ordination import OrdinationResults
>>> or_f = StringIO("Eigvals\t2\n"
...                 "0.0961330159181\t0.0409418140138\n"
...                 "\n"
...                 "Proportion explained\t0\n"
...                 "\n"
...                 "Species\t3\t2\n"
...                 "Species1\t0.408869425742\t0.0695518116298\n"
...                 "Species2\t-0.1153860437\t-0.299767683538\n"
...                 "Species3\t-0.309967102571\t0.187391917117\n"
...                 "\n"
...                 "Site\t3\t2\n"
...                 "Site1\t-0.848956053187\t0.882764759014\n"
...                 "Site2\t-0.220458650578\t-1.34482000302\n"
...                 "Site3\t1.66697179591\t0.470324389808\n"
...                 "\n"
...                 "Biplot\t0\t0\n"
...                 "\n"
...                 "Site constraints\t0\t0\n")
>>> ord_res = OrdinationResults.read(or_f)

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip

import numpy as np

from skbio.stats.ordination import OrdinationResults
from skbio.io import (register_reader, register_writer, register_sniffer,
                      OrdResFormatError)


@register_sniffer('ordres')
def _ordres_sniffer(fh):
    # Smells an ordres file if *all* of the following lines are present *from
    # the beginning* of the file:
    #   - eigvals header (minimally parsed)
    #   - another line (contents ignored)
    #   - a whitespace-only line
    #   - proportion explained header (minimally parsed)
    try:
        _parse_header(fh, 'Eigvals', 1)
        next_line = next(fh, None)

        if next_line is not None:
            _check_empty_line(fh)
            _parse_header(fh, 'Proportion explained', 1)
            return True, {}
    except OrdResFormatError:
        pass

    return False, {}


@register_reader('ordres', OrdinationResults)
def _ordres_to_ordination_results(fh):
    eigvals = _parse_vector_section(fh, 'Eigvals')
    if eigvals is None:
        raise OrdResFormatError("At least one eigval must be present.")
    _check_empty_line(fh)

    prop_expl = _parse_vector_section(fh, 'Proportion explained')
    _check_length_against_eigvals(prop_expl, eigvals,
                                  'proportion explained values')
    _check_empty_line(fh)

    species, species_ids = _parse_array_section(fh, 'Species')
    _check_length_against_eigvals(species, eigvals,
                                  'coordinates per species')
    _check_empty_line(fh)

    site, site_ids = _parse_array_section(fh, 'Site')
    _check_length_against_eigvals(site, eigvals,
                                  'coordinates per site')
    _check_empty_line(fh)

    # biplot does not have ids to parse (the other arrays do)
    biplot, _ = _parse_array_section(fh, 'Biplot', has_ids=False)
    _check_empty_line(fh)

    cons, cons_ids = _parse_array_section(fh, 'Site constraints')

    if cons_ids is not None and site_ids is not None:
        if cons_ids != site_ids:
            raise OrdResFormatError(
                "Site constraints ids and site ids must be equal: %s != %s" %
                (cons_ids, site_ids))

    return OrdinationResults(
        eigvals=eigvals, species=species, site=site, biplot=biplot,
        site_constraints=cons, proportion_explained=prop_expl,
        species_ids=species_ids, site_ids=site_ids)


def _parse_header(fh, header_id, num_dimensions):
    line = next(fh, None)
    if line is None:
        raise OrdResFormatError("Reached end of file while looking for %s "
                                "header." % header_id)

    header = line.strip().split('\t')
    # +1 for the header ID
    if len(header) != num_dimensions + 1 or header[0] != header_id:
        raise OrdResFormatError("%s header not found." % header_id)
    return header


def _check_empty_line(fh):
    """Check that the next line in `fh` is empty or whitespace-only."""
    line = next(fh, None)
    if line is None:
        raise OrdResFormatError("Reached end of file while looking for blank "
                                "line separating sections.")

    if line.strip():
        raise OrdResFormatError("Expected an empty line.")


def _check_length_against_eigvals(data, eigvals, label):
    if data is not None:
        num_vals = data.shape[-1]
        num_eigvals = eigvals.shape[-1]

        if num_vals != num_eigvals:
            raise OrdResFormatError(
                "There should be as many %s as eigvals: %d != %d" %
                (label, num_vals, num_eigvals))


def _parse_vector_section(fh, header_id):
    header = _parse_header(fh, header_id, 1)

    # Parse how many values we are waiting for
    num_vals = int(header[1])
    if num_vals == 0:
        # The ordination method didn't generate the vector, so set it to None
        vals = None
    else:
        # Parse the line with the vector values
        line = next(fh, None)
        if line is None:
            raise OrdResFormatError("Reached end of file while looking for "
                                    "line containing values for %s section."
                                    % header_id)
        vals = np.asarray(line.strip().split('\t'), dtype=np.float64)
        if len(vals) != num_vals:
            raise OrdResFormatError(
                "Expected %d values in %s section, but found %d." %
                (num_vals, header_id, len(vals)))
    return vals


def _parse_array_section(fh, header_id, has_ids=True):
    """Parse an array section of `fh` identified by `header_id`."""
    # Parse the array header
    header = _parse_header(fh, header_id, 2)

    # Parse the dimensions of the array
    rows = int(header[1])
    cols = int(header[2])

    ids = None
    if rows == 0 and cols == 0:
        # The ordination method didn't generate the array data for 'header', so
        # set it to None
        data = None
    elif rows == 0 or cols == 0:
        # Both dimensions should be 0 or none of them are zero
        raise OrdResFormatError("One dimension of %s is 0: %d x %d" %
                                (header_id, rows, cols))
    else:
        # Parse the data
        data = np.empty((rows, cols), dtype=np.float64)

        if has_ids:
            ids = []

        for i in range(rows):
            # Parse the next row of data
            line = next(fh, None)
            if line is None:
                raise OrdResFormatError(
                    "Reached end of file while looking for row %d in %s "
                    "section." % (i + 1, header_id))
            vals = line.strip().split('\t')

            if has_ids:
                ids.append(vals[0])
                vals = vals[1:]

            if len(vals) != cols:
                raise OrdResFormatError(
                    "Expected %d values, but found %d in row %d." %
                    (cols, len(vals), i + 1))
            data[i, :] = np.asarray(vals, dtype=np.float64)
    return data, ids


@register_writer('ordres', OrdinationResults)
def _ordination_results_to_ordres(obj, fh):
    _write_vector_section(fh, 'Eigvals', obj.eigvals)
    _write_vector_section(fh, 'Proportion explained', obj.proportion_explained)
    _write_array_section(fh, 'Species', obj.species, obj.species_ids)
    _write_array_section(fh, 'Site', obj.site, obj.site_ids)
    _write_array_section(fh, 'Biplot', obj.biplot)
    _write_array_section(fh, 'Site constraints', obj.site_constraints,
                         obj.site_ids, include_section_separator=False)


def _write_vector_section(fh, header_id, vector):
    if vector is None:
        shape = 0
    else:
        shape = vector.shape[0]
    fh.write("%s\t%d\n" % (header_id, shape))

    if vector is not None:
        fh.write(_format_vector(vector))
    fh.write("\n")


def _write_array_section(fh, header_id, data, ids=None,
                         include_section_separator=True):
    # write section header
    if data is None:
        shape = (0, 0)
    else:
        shape = data.shape
    fh.write("%s\t%d\t%d\n" % (header_id, shape[0], shape[1]))

    # write section data
    if data is not None:
        if ids is None:
            for vals in data:
                fh.write(_format_vector(vals))
        else:
            for id_, vals in zip(ids, data):
                fh.write(_format_vector(vals, id_))

    if include_section_separator:
        fh.write("\n")


def _format_vector(vector, id_=None):
    formatted_vector = '\t'.join(np.asarray(vector, dtype=np.str))

    if id_ is None:
        return "%s\n" % formatted_vector
    else:
        return "%s\t%s\n" % (id_, formatted_vector)
