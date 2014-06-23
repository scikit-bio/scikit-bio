#!/usr/bin/env python
"""Retrieves records by id from RFAM, the RNA families database."""
from cogent.db.util import UrlGetter

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight","Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

rfam_base='http://rfam.sanger.ac.uk/family/alignment/download/format?'
rfam_formats = dict.fromkeys('stockholm pfam fasta fastau')
rfam_types = dict.fromkeys(['seed','full'])

class Rfam(UrlGetter):
    """Returns a pdb file."""
    Defaults={'alnType':'seed','format':'stockholm','acc':None,'nseLabels':'1',
        'download':'0'}
    PrintedFields=dict.fromkeys('acc alnType nseLabels format download'.split())
    BaseUrl = rfam_base

    def __getitem__(self, item):
        """Returns handle to file containing aln of specified rfam id."""
        orig_acc = self.acc
        item = str(item).upper()
        if not item.startswith('RF'):
            item = 'RF'+item
        
        self.acc = item
        result = self.open()
        self.acc = orig_acc
        return result
