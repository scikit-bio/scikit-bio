#!/usr/bin/env python
"""Tests of data retrieval from PDB."""
from cogent.util.unit_test import TestCase, main
from cogent.db.rfam import Rfam

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight","Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class RfamTests(TestCase):
    """Tests of the Rfam class."""
    def test_simple_get(self):
        """Simple access of an item should work."""
        rec = Rfam()
        result = rec['rf00100'].read()
        assert result.startswith('# STOCKHOLM')
        assert 'AM773434.1/1-304' in result

if __name__ == '__main__':
    main()
