r"""Ordination results format (:mod:`skbio.io.format.ordination`)
=============================================================

.. currentmodule:: skbio.io.format.ordination

The ordination results file format (``ordination``) stores the results of an
ordination method in a human-readable, text-based format. The format supports
storing the results of various ordination methods available in scikit-bio,
including (but not necessarily limited to) PCoA, CA, RDA, and CCA.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.stats.ordination.OrdinationResults`                |
+------+------+---------------------------------------------------------------+

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

    Eigvals<tab>4
    0.36<tab>0.18<tab>0.07<tab>0.08

    Proportion explained<tab>4
    0.46<tab>0.23<tab>0.10<tab>0.10

    Species<tab>9<tab>4
    Species0<tab>0.11<tab>0.28<tab>-0.20<tab>-0.00
    Species1<tab>0.14<tab>0.30<tab>0.39<tab>-0.14
    Species2<tab>-1.01<tab>0.09<tab>-0.19<tab>-0.10
    Species3<tab>-1.03<tab>0.10<tab>0.22<tab>0.22
    Species4<tab>1.05<tab>0.53<tab>-0.43<tab>0.22
    Species5<tab>0.99<tab>0.57<tab>0.67<tab>-0.38
    Species6<tab>0.25<tab>-0.17<tab>-0.20<tab>0.43
    Species7<tab>0.14<tab>-0.85<tab>-0.01<tab>0.05
    Species8<tab>0.41<tab>-0.70<tab>0.21<tab>-0.69

    Site<tab>10<tab>4
    Site0<tab>0.71<tab>-3.08<tab>0.21<tab>-1.24
    Site1<tab>0.58<tab>-3.00<tab>-0.94<tab>2.69
    Site2<tab>0.76<tab>-3.15<tab>2.13<tab>-3.11
    Site3<tab>1.11<tab>1.07<tab>-1.87<tab>0.66
    Site4<tab>-0.97<tab>-0.06<tab>-0.69<tab>-0.61
    Site5<tab>1.04<tab>0.45<tab>-0.63<tab>0.28
    Site6<tab>-0.95<tab>-0.08<tab>0.13<tab>-0.42
    Site7<tab>0.94<tab>-0.10<tab>0.52<tab>-0.00
    Site8<tab>-1.14<tab>0.49<tab>0.47<tab>1.17
    Site9<tab>1.03<tab>1.03<tab>2.74<tab>-1.28

    Biplot<tab>3<tab>3
    -0.16<tab>0.63<tab>0.76
    -0.99<tab>0.06<tab>-0.04
    0.18<tab>-0.97<tab>0.03

    Site constraints<tab>10<tab>4
    Site0<tab>0.69<tab>-3.08<tab>-0.32<tab>-1.24
    Site1<tab>0.66<tab>-3.06<tab>0.23<tab>2.69
    Site2<tab>0.63<tab>-3.04<tab>0.78<tab>-3.11
    Site3<tab>1.10<tab>0.50<tab>-1.55<tab>0.66
    Site4<tab>-0.97<tab>0.06<tab>-1.12<tab>-0.61
    Site5<tab>1.05<tab>0.53<tab>-0.43<tab>0.28
    Site6<tab>-1.02<tab>0.10<tab>-0.00<tab>-0.42
    Site7<tab>0.99<tab>0.57<tab>0.67<tab>-0.00
    Site8<tab>-1.08<tab>0.13<tab>1.11<tab>1.17
    Site9<tab>0.94<tab>0.61<tab>1.79<tab>-1.28


If a given result attribute is not present (e.g. ``Biplot``), it should still
be defined and declare its dimensions as 0. For example::

    Biplot<tab>0<tab>0

All attributes are optional except for ``Eigvals``.

Examples
--------
Assume we have the following tab-delimited text file storing the
ordination results in ``ordination`` format::

    Eigvals<tab>4
    0.36<tab>0.18<tab>0.07<tab>0.08

    Proportion explained<tab>4
    0.46<tab>0.23<tab>0.10<tab>0.10

    Species<tab>9<tab>4
    Species0<tab>0.11<tab>0.28<tab>-0.20<tab>-0.00
    Species1<tab>0.14<tab>0.30<tab>0.39<tab>-0.14
    Species2<tab>-1.01<tab>0.09<tab>-0.19<tab>-0.10
    Species3<tab>-1.03<tab>0.10<tab>0.22<tab>0.22
    Species4<tab>1.05<tab>0.53<tab>-0.43<tab>0.22
    Species5<tab>0.99<tab>0.57<tab>0.67<tab>-0.38
    Species6<tab>0.25<tab>-0.17<tab>-0.20<tab>0.43
    Species7<tab>0.14<tab>-0.85<tab>-0.01<tab>0.05
    Species8<tab>0.41<tab>-0.70<tab>0.21<tab>-0.69

    Site<tab>10<tab>4
    Site0<tab>0.71<tab>-3.08<tab>0.21<tab>-1.24
    Site1<tab>0.58<tab>-3.00<tab>-0.94<tab>2.69
    Site2<tab>0.76<tab>-3.15<tab>2.13<tab>-3.11
    Site3<tab>1.11<tab>1.07<tab>-1.87<tab>0.66
    Site4<tab>-0.97<tab>-0.06<tab>-0.69<tab>-0.61
    Site5<tab>1.04<tab>0.45<tab>-0.63<tab>0.28
    Site6<tab>-0.95<tab>-0.08<tab>0.13<tab>-0.42
    Site7<tab>0.94<tab>-0.10<tab>0.52<tab>-0.00
    Site8<tab>-1.14<tab>0.49<tab>0.47<tab>1.17
    Site9<tab>1.03<tab>1.03<tab>2.74<tab>-1.28

    Biplot<tab>0<tab>0

    Site constraints<tab>0<tab>0

Load the ordination results from the file:

>>> from io import StringIO
>>> from skbio import OrdinationResults
>>> or_f = StringIO(
...  "Eigvals\t4\n"
...  "0.36\t0.18\t0.07\t0.08\n"
...  "\n"
...  "Proportion explained\t4\n"
...  "0.46\t0.23\t0.10\t0.10\n"
...  "\n"
...  "Species\t9\t4\n"
...  "Species0\t0.11\t0.28\t-0.20\t-0.00\n"
...  "Species1\t0.14\t0.30\t0.39\t-0.14\n"
...  "Species2\t-1.01\t0.09\t-0.19\t-0.10\n"
...  "Species3\t-1.03\t0.10\t0.22\t0.22\n"
...  "Species4\t1.05\t0.53\t-0.43\t0.22\n"
...  "Species5\t0.99\t0.57\t0.67\t-0.38\n"
...  "Species6\t0.25\t-0.17\t-0.20\t0.43\n"
...  "Species7\t0.14\t-0.85\t-0.01\t0.05\n"
...  "Species8\t0.41\t-0.70\t0.21\t-0.69\n"
...  "\n"
...  "Site\t10\t4\n"
...  "Site0\t0.71\t-3.08\t0.21\t-1.24\n"
...  "Site1\t0.58\t-3.00\t-0.94\t2.69\n"
...  "Site2\t0.76\t-3.15\t2.13\t-3.11\n"
...  "Site3\t1.11\t1.07\t-1.87\t0.66\n"
...  "Site4\t-0.97\t-0.06\t-0.69\t-0.61\n"
...  "Site5\t1.04\t0.45\t-0.63\t0.28\n"
...  "Site6\t-0.95\t-0.08\t0.13\t-0.42\n"
...  "Site7\t0.94\t-0.10\t0.52\t-0.00\n"
...  "Site8\t-1.14\t0.49\t0.47\t1.17\n"
...  "Site9\t1.03\t1.03\t2.74\t-1.28\n"
...  "\n"
...  "Biplot\t0\t0\n"
...  "\n"
...  "Site constraints\t0\t0\n")
>>> ord_res = OrdinationResults.read(or_f)


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd

from skbio.stats.ordination import OrdinationResults
from skbio.io import create_format, OrdinationFormatError

ordination = create_format("ordination")


@ordination.sniffer()
def _ordination_sniffer(fh):
    # Smells an ordination file if *all* of the following lines are present
    # *from the beginning* of the file:
    #   - eigvals header (minimally parsed)
    #   - another line (contents ignored)
    #   - a whitespace-only line
    #   - proportion explained header (minimally parsed)
    try:
        _parse_header(fh, "Eigvals", 1)
        next_line = next(fh, None)

        if next_line is not None:
            _check_empty_line(fh)
            _parse_header(fh, "Proportion explained", 1)
            return True, {}
    except OrdinationFormatError:
        pass

    return False, {}


@ordination.reader(OrdinationResults)
def _ordination_to_ordination_results(fh):
    eigvals = _parse_vector_section(fh, "Eigvals")
    if eigvals is None:
        raise OrdinationFormatError("At least one eigval must be present.")
    _check_empty_line(fh)

    prop_expl = _parse_vector_section(fh, "Proportion explained")
    _check_length_against_eigvals(prop_expl, eigvals, "proportion explained values")
    _check_empty_line(fh)

    species = _parse_array_section(fh, "Species")
    _check_length_against_eigvals(species, eigvals, "coordinates per species")
    _check_empty_line(fh)

    site = _parse_array_section(fh, "Site")
    _check_length_against_eigvals(site, eigvals, "coordinates per site")
    _check_empty_line(fh)

    # biplot does not have ids to parse (the other arrays do)
    biplot = _parse_array_section(fh, "Biplot", has_ids=False)
    _check_empty_line(fh)

    cons = _parse_array_section(fh, "Site constraints")

    if cons is not None and site is not None:
        if not np.array_equal(cons.index, site.index):
            raise OrdinationFormatError(
                "Site constraints ids and site ids must be equal: %s != %s"
                % (cons.index, site.index)
            )

    return OrdinationResults(
        short_method_name="",
        long_method_name="",
        eigvals=eigvals,
        features=species,
        samples=site,
        biplot_scores=biplot,
        sample_constraints=cons,
        proportion_explained=prop_expl,
    )


def _parse_header(fh, header_id, num_dimensions):
    line = next(fh, None)
    if line is None:
        raise OrdinationFormatError(
            "Reached end of file while looking for %s header." % header_id
        )

    header = line.strip().split("\t")
    # +1 for the header ID
    if len(header) != num_dimensions + 1 or header[0] != header_id:
        raise OrdinationFormatError("%s header not found." % header_id)
    return header


def _check_empty_line(fh):
    """Check that the next line in `fh` is empty or whitespace-only."""
    line = next(fh, None)
    if line is None:
        raise OrdinationFormatError(
            "Reached end of file while looking for blank line separating " "sections."
        )

    if line.strip():
        raise OrdinationFormatError("Expected an empty line.")


def _check_length_against_eigvals(data, eigvals, label):
    if data is not None:
        num_vals = data.shape[-1]
        num_eigvals = eigvals.shape[-1]

        if num_vals != num_eigvals:
            raise OrdinationFormatError(
                "There should be as many %s as eigvals: %d != %d"
                % (label, num_vals, num_eigvals)
            )


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
            raise OrdinationFormatError(
                "Reached end of file while looking for line containing values "
                "for %s section." % header_id
            )
        vals = pd.Series(np.asarray(line.strip().split("\t"), dtype=np.float64))
        if len(vals) != num_vals:
            raise OrdinationFormatError(
                "Expected %d values in %s section, but found %d."
                % (num_vals, header_id, len(vals))
            )
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
        raise OrdinationFormatError(
            "One dimension of %s is 0: %d x %d" % (header_id, rows, cols)
        )
    else:
        # Parse the data
        data = np.empty((rows, cols), dtype=np.float64)

        if has_ids:
            ids = []

        for i in range(rows):
            # Parse the next row of data
            line = next(fh, None)
            if line is None:
                raise OrdinationFormatError(
                    "Reached end of file while looking for row %d in %s "
                    "section." % (i + 1, header_id)
                )
            vals = line.strip().split("\t")

            if has_ids:
                ids.append(vals[0])
                vals = vals[1:]

            if len(vals) != cols:
                raise OrdinationFormatError(
                    "Expected %d values, but found %d in row %d."
                    % (cols, len(vals), i + 1)
                )
            data[i, :] = np.asarray(vals, dtype=np.float64)
        data = pd.DataFrame(data, index=ids)

    return data


@ordination.writer(OrdinationResults)
def _ordination_results_to_ordination(obj, fh):
    _write_vector_section(fh, "Eigvals", obj.eigvals)
    _write_vector_section(fh, "Proportion explained", obj.proportion_explained)
    _write_array_section(fh, "Species", obj.features)
    _write_array_section(fh, "Site", obj.samples)
    _write_array_section(fh, "Biplot", obj.biplot_scores, has_ids=False)
    _write_array_section(
        fh, "Site constraints", obj.sample_constraints, include_section_separator=False
    )


def _write_vector_section(fh, header_id, vector):
    if vector is None:
        shape = 0
    else:
        shape = vector.shape[0]
    fh.write("%s\t%d\n" % (header_id, shape))

    if vector is not None:
        fh.write(_format_vector(vector.values))
    fh.write("\n")


def _write_array_section(
    fh, header_id, data, has_ids=True, include_section_separator=True
):
    # write section header
    if data is None:
        shape = (0, 0)
    else:
        shape = data.shape
    fh.write("%s\t%d\t%d\n" % (header_id, shape[0], shape[1]))

    # write section data
    if data is not None:
        if not has_ids:
            for vals in data.values:
                fh.write(_format_vector(vals))
        else:
            for id_, vals in zip(data.index, data.values):
                fh.write(_format_vector(vals, id_))

    if include_section_separator:
        fh.write("\n")


def _format_vector(vector, id_=None):
    formatted_vector = "\t".join(np.asarray(vector, dtype=str))

    if id_ is None:
        return "%s\n" % formatted_vector
    else:
        return "%s\t%s\n" % (id_, formatted_vector)
