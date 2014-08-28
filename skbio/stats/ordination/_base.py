from __future__ import absolute_import, division, print_function
from future.builtins import zip

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.io import FileFormatError
from skbio.io.util import open_file


class OrdinationResults(object):
    """Store ordination results

    Attributes
    ----------
    eigvals : 1-D numpy array
        The result eigenvalues
    species : 2-D numpy array
        The result coordinates for each species
    site : 2-D numpy array
        The results coordinates for each site
    biplot : 2-D numpy array
        The result biplot coordinates
    site_constraints : 2-D numpy array
        The result coordinates for each site constraint
    proportion_explained : 1-D numpy array
        The proportion explained by each eigenvector
    species_ids : list of str
        The species identifiers
    site_ids : list of str
        The site identifiers

    """

    def __init__(self, eigvals, species=None, site=None, biplot=None,
                 site_constraints=None, proportion_explained=None,
                 species_ids=None, site_ids=None):
        self.eigvals = eigvals
        self.species = species
        self.site = site
        self.biplot = biplot
        self.site_constraints = site_constraints
        self.proportion_explained = proportion_explained
        self.species_ids = species_ids
        self.site_ids = site_ids

    @classmethod
    def read(cls, fp, **kwargs):
        """Load ordination results from file.

        Creates an ``OrdinationResults`` instance from a supported file format.

        Supported file formats include:

        - ``ordres`` (:mod:`skbio.io.ordres`)

        Parameters
        ----------
        fp : filepath or filehandle
            File to read from.
        kwargs : dict, optional
            Keyword arguments passed to :mod:`skbio.io.read` and the file
            format reader.

        Returns
        -------
        OrdinationResults
            Instance of type `cls` containing the parsed contents of `fp`.

        See Also
        --------
        write
        skbio.io.ordres
        skbio.io.read

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
        return skbio.io.read(fp, into=cls, **kwargs)

    def write(self, fp, format='ordres', **kwargs):
        """Save ordination results to file.

        Supported file formats include:

        - ``ordres`` (:mod:`skbio.io.ordres`)

        Parameters
        ----------
        fp : filepath or filehandle
            File to write to.
        format : str, optional
            File format to write.
        kwargs : dict, optional
            Keyword arguments passed to :mod:`skbio.io.write` and the file
            format writer.

        See Also
        --------
        read
        skbio.io.ordres
        skbio.io.write

        """
        skbio.io.write(self, into=fp, format=format, **kwargs)

    @classmethod
    def from_file(cls, ord_res_f):
        r"""Load ordination results from text file.

        .. note:: Deprecated in scikit-bio 0.2.0-dev
           ``from_file`` will be removed in scikit-bio 0.3.0. It is replaced by
           ``read``, which is a more general method for deserializing
           ordination results. ``read`` supports multiple file formats,
           automatic file format detection, etc. by taking advantage of
           scikit-bio's I/O registry system. See :mod:`skbio.io` for more
           details.

        Creates an `OrdinationResults` instance from a ``ordres`` formatted
        file. See :mod:`skbio.io.ordres` for the format specification.

        Parameters
        ----------
        ord_res_f: filepath or filehandle
            File to read from.

        Returns
        -------
        OrdinationResults
            Instance of type `cls` containing the parsed contents of
            `ord_res_f`.

        Raises
        ------
        OrdResFormatError
            If the format of the file is not valid, or if the shapes of the
            different sections of the file are not consistent.

        See Also
        --------
        read

        """
        warnings.warn(
            "OrdinationResults.from_file is deprecated and will be removed in "
            "scikit-bio 0.3.0. Please update your code to use "
            "OrdinationResults.read.", UserWarning)
        return cls.read(ord_res_f, format='ordres')

    def to_file(self, out_f):
        """Save ordination results to file in text format.

        .. note:: Deprecated in scikit-bio 0.2.0-dev
           ``to_file`` will be removed in scikit-bio 0.3.0. It is replaced by
           ``write``, which is a more general method for serializing ordination
           results. ``write`` supports multiple file formats by taking
           advantage of scikit-bio's I/O registry system. See :mod:`skbio.io`
           for more details.

        Serializes ordination results as an ``ordres`` formatted file. See
        :mod:`skbio.io.ordres` for the format specification.

        Parameters
        ----------
        out_f : filepath or filehandle
            File to write to.

        See Also
        --------
        write

        """
        warnings.warn(
            "OrdinationResults.to_file is deprecated and will be removed in "
            "scikit-bio 0.3.0. Please update your code to use "
            "OrdinationResults.write.", UserWarning)
        self.write(out_f, format='ordres')


class Ordination(object):
    short_method_name = 'Overwrite in subclass!'
    long_method_name = 'Overwrite in subclass!'
