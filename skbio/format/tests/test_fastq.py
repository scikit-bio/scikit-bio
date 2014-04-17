#!/usr/bin/env python

import numpy as np
from unittest import TestCase, main
from skbio.format.fastq import format_fastq_record


class FASTQFormatTests(TestCase):
    def test_format_fastq_record(self):
        exp = "@abc\ndef\n+\nfgh\n"
        obs = format_fastq_record('abc', 'def',
                                  np.array([38, 39, 40], dtype=np.int8))
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
