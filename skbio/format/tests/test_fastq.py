#!/usr/bin/env python

import numpy as np
from unittest import TestCase, main
from skbio.format.fastq import (format_fastq_record, _phred_to_ascii33,
                                _phred_to_ascii64)


class FASTQFormatTests(TestCase):
    def test_format_fastq_record(self):
        """Construt a FASTQ record"""
        exp = b"@abc\ndef\n+\nfgh\n"
        obs = format_fastq_record(b'abc', b'def',
                                  np.array([38, 39, 40], dtype=np.int8), 64)
        self.assertEqual(obs, exp)

    def test_phred_to_ascii33(self):
        """Write out terrible FASTQ quality scores"""
        exp = b'GHI'
        obs = _phred_to_ascii33(np.array([38, 39, 40], dtype=np.int8))
        self.assertEqual(obs, exp)

    def test_phred_to_ascii64(self):
        """Write out terrible FASTQ quality scores"""
        exp = b'fgh'
        obs = _phred_to_ascii64(np.array([38, 39, 40], dtype=np.int8))
        self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
