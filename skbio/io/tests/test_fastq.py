# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip
from six import StringIO

from unittest import TestCase, main

from skbio import (BiologicalSequence, NucleotideSequence, DNA, RNA, Protein,
                   SequenceCollection, Alignment)
from skbio.sequence import BiologicalSequenceError
from skbio.io import FASTQFormatError
from skbio.io.fastq import (
    _fastq_sniffer, _fastq_to_generator, _fastq_to_biological_sequence,
    _fastq_to_nucleotide_sequence, _fastq_to_dna_sequence,
    _fastq_to_rna_sequence, _fastq_to_protein_sequence,
    _fastq_to_sequence_collection, _fastq_to_alignment, _generator_to_fastq,
    _biological_sequence_to_fastq, _nucleotide_sequence_to_fastq,
    _dna_sequence_to_fastq, _rna_sequence_to_fastq, _protein_sequence_to_fastq,
    _sequence_collection_to_fastq, _alignment_to_fastq)

from skbio.util import get_data_path

# Note: the example FASTQ files with file extension .fastq are taken from the
# following open-access publication's supplementary data:
#
# P.J.A. Cock, C.J. Fields, N. Goto, M.L. Heuer and P.M. Rice (2009). The
# Sanger FASTQ file format for sequences with quality scores, and the
# Solexa/Illumina FASTQ variants.
#
# See licenses/fastq-example-files-readme.txt for the original README that
# accompanied these files, which includes the terms of use.
#
# The example files have not been modified from their original form.


class ReaderTests(TestCase):
    def setUp(self):
        # empty file shouldn't yield sequences
        self.empty = ([], {'variant': 'sanger'}, [get_data_path('empty')])

        self.invalid_fps = [(get_data_path(e[0]), e[1], e[2]) for e in [
            ('error_diff_ids.fastq', FASTQFormatError,
             "header lines do not match: "
             "'SLXA-B3_649_FC8437_R1_1_1_850_123' != "
             "'SLXA-B3_649_FC8437_R1_1_1_850_124'"),

            ('error_double_qual.fastq', FASTQFormatError,
             "Extra quality.*'\+SLXA-B3_649_FC8437_R1_1_1_850_123'"),

            ('error_double_seq.fastq', FASTQFormatError,
             "FASTQ record that is missing a quality \(\+\) header line"),

            ('error_long_qual.fastq', FASTQFormatError, "Extra quality.*'Y'"),

            ('error_no_qual.fastq', FASTQFormatError,
             'blank line.*FASTQ'),

            ('error_qual_del.fastq', ValueError, '94.*[0, 93]'),

            ('error_qual_escape.fastq', ValueError, '-6.*[0, 93]'),

            ('error_qual_null.fastq', ValueError, '-33.*[0, 93]'),

            ('error_qual_space.fastq', ValueError, '-1.*[0, 93]'),

            ('error_qual_tab.fastq', FASTQFormatError,
             'blank line.*FASTQ'),

            ('error_qual_unit_sep.fastq', ValueError, '-2.*[0, 93]'),

            ('error_qual_vtab.fastq', ValueError, '-22.*[0, 93]'),

            ('error_short_qual.fastq', FASTQFormatError,
             "Extra quality.*'SLXA-B3_649_FC8437_R1_1_1_362_549'"),

            ('error_spaces.fastq', FASTQFormatError,
             "whitespace.*sequence data: 'GATGTGCAA TACCTTTGTA GAGGAA'"),

            ('error_tabs.fastq', FASTQFormatError,
             r"whitespace.*sequence data: 'GATGTGCAA\\tTACCTTTGTA\\tGAGGAA'"),

            ('error_trunc_at_seq.fastq', FASTQFormatError,
             'blank line.*FASTQ'),

            ('error_trunc_at_plus.fastq', FASTQFormatError,
             'blank line.*FASTQ'),

            ('error_trunc_at_qual.fastq', FASTQFormatError,
             "incomplete/truncated.*end of file.*missing quality scores"),

            ('error_trunc_in_title.fastq', FASTQFormatError,
             "incomplete/truncated.*end of file.*missing sequence data"),

            ('error_trunc_in_seq.fastq', FASTQFormatError,
             "incomplete/truncated.*end of file.*missing a quality.*header"),

            ('error_trunc_in_plus.fastq', FASTQFormatError,
             "header lines do not match: "
             "'SLXA-B3_649_FC8437_R1_1_1_183_714' != 'SLXA-B3_649_FC'"),

            ('error_trunc_in_qual.fastq', FASTQFormatError,
             'end of file.*different number.*24 != 25')
        ]]

    def test_fasta_to_generator_valid_files(self):
        test_cases = (self.empty,)

        for exp, kwargs, fps in test_cases:
            for fp in fps:
                obs = list(_fastq_to_generator(fp, **kwargs))
                self._assert_generator_results_equal(obs, exp)

    def test_fastq_to_generator_invalid_files(self):
        for fp, error_type, error_msg_regex in self.invalid_fps:
            with self.assertRaisesRegexp(error_type, error_msg_regex):
                # TODO test that errors are raised for all variants
                list(_fastq_to_generator(fp, variant='sanger'))

    def _assert_generator_results_equal(self, obs, exp):
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertTrue(o.equals(e))


if __name__ == '__main__':
    main()
