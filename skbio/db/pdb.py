#!/usr/bin/env python
"""Retrieves records by id from PDB, the Protein Data Bank."""
from cogent.db.util import UrlGetter

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

pdb_base='http://www.rcsb.org/pdb/files/'

class Pdb(UrlGetter):
    """Returns a pdb file."""
    BaseUrl = pdb_base
    Suffix='.pdb'
    Key=None

    def __str__(self):
        return self.BaseUrl + str(self.Key) + self.Suffix

    def __getitem__(self, item):
        """Returns handle to file containing specified PDB id."""
        orig_key = self.Key
        self.Key = item.lower()
        result = self.open()
        self.Key = orig_key
        return result
