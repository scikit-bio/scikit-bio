# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import pandas as pd
import tempfile
import os


from skbio.util import get_data_path
from skbio.util._decorator import overrides

from skbio.metadata._metadata import SampleMetadata
from skbio.io.format.sample_metadata import (
    _sample_metadata_sniffer, _sample_metadata_read)

class SnifferTests(TestCase):
    def setUp(self):
        self.positive_fps = list(map(get_data_path, [
            'sample-metadata-comments-comment-char-id.tsv',
            'sample-metadata-comments.tsv',
            'sample-metadata-comments-mixed-case.tsv',
            'sample-metadata-complete-types-directive.tsv',
            'sample-metadata-empty-rows.tsv',
            'sample-metadata-leading-trailing-whitespace.tsv',
            'sample-metadata-leading-trailing-whitespace-split-id.tsv',
        ]))

        self.negative_fps = list(map(get_data_path, [
            'empty',
            'whitespace_only',
        ]))

    def test_positives(self):
        for fp in self.positive_fps:
            self.assertEqual(_sample_metadata_sniffer(fp), (True, {}))

    def test_negatives(self):
        for fp in self.negative_fps:
            self.assertEqual(_sample_metadata_sniffer(fp), (False, {}))


class TestSampleMetadataReader(TestCase):
    def test_reader(self):
        fp = get_data_path('sample-metadata-complete-types-directive.tsv')
        obs_md = SampleMetadata.read(fp)
        exp_index = pd.Index(['id1', 'id2', 'id3'], name='id')
        exp_df = pd.DataFrame({'col1': ['1', '2', '3'],
                               'col2': ['a', 'b', 'c'],
                               'col3': ['foo', 'bar', '42']},
                              index=exp_index)
        exp_md = SampleMetadata(exp_df)
        self.assertEqual(obs_md, exp_md)


class TestRoundtrip(TestCase):
    def setUp(self):
        self.temp_dir_obj = tempfile.TemporaryDirectory(
            prefix='sample-metadata-temp')
        self.temp_dir = self.temp_dir_obj.name

        self.filepath = os.path.join(self.temp_dir, 'metadata.tsv')

    def tearDown(self):
        self.temp_dir_obj.cleanup()

    def test_simple(self):
        fp = get_data_path('sample-metadata-comments-mixed-case.tsv')
        md1 = SampleMetadata.read(fp)
        md1.write(self.filepath)
        md2 = SampleMetadata.read(self.filepath)
        self.assertEqual(md1, md2)


if __name__ == '__main__':
    main()
