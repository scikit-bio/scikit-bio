import collections
import unittest
import warnings

import pandas as pd
import numpy as np

from skbio.metadata._metadata import (SampleMetadata, CategoricalMetadataColumn,
                                      NumericMetadataColumn)


class TestInvalidMetadataConstruction(unittest.TestCase):
    def test_non_dataframe(self):
        with self.assertRaisesRegex(
                TypeError, 'Metadata constructor.*DataFrame.*not.*Series'):
            SampleMetadata(pd.Series([1, 2, 3], name='col',
                               index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_no_ids(self):
        with self.assertRaisesRegex(ValueError, 'Metadata.*at least one ID'):
            SampleMetadata(pd.DataFrame({}, index=pd.Index([], name='id')))

        with self.assertRaisesRegex(ValueError, 'Metadata.*at least one ID'):
            SampleMetadata(pd.DataFrame({'column': []},
                                  index=pd.Index([], name='id')))

    def test_invalid_id_header(self):
        # default index name
        with self.assertRaisesRegex(ValueError, r'Index\.name.*None'):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3]}, index=pd.Index(['a', 'b', 'c'])))

        with self.assertRaisesRegex(ValueError, r'Index\.name.*my-id-header'):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', 'b', 'c'], name='my-id-header')))

    def test_non_str_id(self):
        with self.assertRaisesRegex(
                TypeError, 'non-string metadata ID.*type.*float.*nan'):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', np.nan, 'c'], name='id')))

    def test_non_str_column_name(self):
        with self.assertRaisesRegex(
                TypeError, 'non-string metadata column name.*type.*'
                           'float.*nan'):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3],
                 np.nan: [4, 5, 6]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_empty_id(self):
        with self.assertRaisesRegex(
                ValueError, 'empty metadata ID.*at least one character'):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3]}, index=pd.Index(['a', '', 'c'], name='id')))

    def test_empty_column_name(self):
        with self.assertRaisesRegex(
                ValueError, 'empty metadata column name.*'
                            'at least one character'):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3],
                 '': [4, 5, 6]}, index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_pound_sign_id(self):
        with self.assertRaisesRegex(
                ValueError, "metadata ID.*begins with a pound sign.*'#b'"):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', '#b', 'c'], name='id')))

    def test_id_conflicts_with_id_header(self):
        with self.assertRaisesRegex(
                ValueError, "metadata ID 'sample-id'.*conflicts.*reserved.*"
                            "ID header"):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', 'sample-id', 'c'], name='id')))

    def test_column_name_conflicts_with_id_header(self):
        with self.assertRaisesRegex(
                ValueError, "metadata column name 'featureid'.*conflicts.*"
                            "reserved.*ID header"):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3],
                 'featureid': [4, 5, 6]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_duplicate_ids(self):
        with self.assertRaisesRegex(ValueError, "Metadata IDs.*unique.*'a'"):
            SampleMetadata(pd.DataFrame(
                {'col': [1, 2, 3]},
                index=pd.Index(['a', 'b', 'a'], name='id')))

    def test_duplicate_column_names(self):
        data = [[1, 2, 3],
                [4, 5, 6],
                [7, 8, 9]]
        with self.assertRaisesRegex(ValueError,
                                    "Metadata column names.*unique.*'col1'"):
            SampleMetadata(pd.DataFrame(data, columns=['col1', 'col2', 'col1'],
                                  index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_unsupported_column_dtype(self):
        with self.assertRaisesRegex(
                TypeError, "Metadata column 'col2'.*unsupported.*dtype.*bool"):
            SampleMetadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 'col2': [True, False, True]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_categorical_column_unsupported_type(self):
        with self.assertRaisesRegex(
                TypeError, "CategoricalMetadataColumn.*strings or missing "
                           r"values.*42\.5.*float.*'col2'"):
            SampleMetadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 'col2': ['foo', 'bar', 42.5]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_categorical_column_empty_str(self):
        with self.assertRaisesRegex(
                ValueError, "CategoricalMetadataColumn.*empty strings.*"
                            "column 'col2'"):
            SampleMetadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 'col2': ['foo', '', 'bar']},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_numeric_column_infinity(self):
        with self.assertRaisesRegex(
                ValueError, "NumericMetadataColumn.*positive or negative "
                            "infinity.*column 'col2'"):
            SampleMetadata(pd.DataFrame(
                {'col1': ['foo', 'bar', 'baz'],
                 'col2': [42, float('+inf'), 4.3]},
                index=pd.Index(['a', 'b', 'c'], name='id')))

    def test_unknown_missing_scheme(self):
        with self.assertRaisesRegex(ValueError, "BAD:SCHEME"):
            SampleMetadata(pd.DataFrame(
                {'col1': [1, 2, 3],
                 'col2': ['foo', 'bar', 'bar']},
                index=pd.Index(['a', 'b', 'c'], name='id')),
                default_missing_scheme='BAD:SCHEME')

    def test_missing_q2_error(self):
        index = pd.Index(['None', 'nan', 'NA', 'foo'], name='id')
        df = pd.DataFrame(collections.OrderedDict([
            ('col1', [1.0, np.nan, np.nan, np.nan]),
            ('NA', [np.nan, np.nan, np.nan, np.nan]),
            ('col3', ['null', 'N/A', np.nan, 'NA']),
            ('col4', np.array([np.nan, np.nan, np.nan, np.nan],
                              dtype=object))]),
            index=index)

        with self.assertRaisesRegex(ValueError, 'col1.*no-missing'):
            SampleMetadata(df, default_missing_scheme='no-missing')


class TestMetadataConstructionAndProperties(unittest.TestCase):
    def assertEqualColumns(self, obs_columns, exp):
        obs = [(name, props.type) for name, props in obs_columns.items()]
        self.assertEqual(obs, exp)

    def test_minimal(self):
        md = SampleMetadata(pd.DataFrame({}, index=pd.Index(['a'], name='id')))

        self.assertEqual(md.id_count, 1)
        self.assertEqual(md.column_count, 0)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('a',))
        self.assertEqualColumns(md.columns, [])

    def test_single_id(self):
        index = pd.Index(['id1'], name='id')
        df = pd.DataFrame({'col1': [1.0], 'col2': ['a'], 'col3': ['foo']},
                          index=index)
        md = SampleMetadata(df)

        self.assertEqual(md.id_count, 1)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1',))
        self.assertEqualColumns(md.columns,
                                [('col1', 'numeric'), ('col2', 'categorical'),
                                 ('col3', 'categorical')])

    def test_no_columns(self):
        index = pd.Index(['id1', 'id2', 'foo'], name='id')
        df = pd.DataFrame({}, index=index)
        md = SampleMetadata(df)

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 0)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1', 'id2', 'foo'))
        self.assertEqualColumns(md.columns, [])

    def test_single_column(self):
        index = pd.Index(['id1', 'a', 'my-id'], name='id')
        df = pd.DataFrame({'column': ['foo', 'bar', 'baz']}, index=index)
        md = SampleMetadata(df)

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 1)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1', 'a', 'my-id'))
        self.assertEqualColumns(md.columns, [('column', 'categorical')])

    def test_retains_column_order(self):
        # Supply DataFrame constructor with explicit column ordering instead of
        # a dict.
        index = pd.Index(['id1', 'id2', 'id3'], name='id')
        columns = ['z', 'a', 'ch']
        data = [
            [1.0, 'a', 'foo'],
            [2.0, 'b', 'bar'],
            [3.0, 'c', '42']
        ]
        df = pd.DataFrame(data, index=index, columns=columns)
        md = SampleMetadata(df)

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1', 'id2', 'id3'))
        self.assertEqualColumns(md.columns,
                                [('z', 'numeric'), ('a', 'categorical'),
                                 ('ch', 'categorical')])

    def test_supported_id_headers(self):
        case_insensitive = {
            'id', 'sampleid', 'sample id', 'sample-id', 'featureid',
            'feature id', 'feature-id'
        }

        exact_match = {
            '#SampleID', '#Sample ID', '#OTUID', '#OTU ID', 'sample_name'
        }

        # Build a set of supported headers, including exact matches and headers
        # with different casing.
        headers = set()
        for header in case_insensitive:
            headers.add(header)
            headers.add(header.upper())
            headers.add(header.title())
        for header in exact_match:
            headers.add(header)

        count = 0
        for header in headers:
            index = pd.Index(['id1', 'id2'], name=header)
            df = pd.DataFrame({'column': ['foo', 'bar']}, index=index)
            md = SampleMetadata(df)

            self.assertEqual(md.id_header, header)
            count += 1

        # Since this test case is a little complicated, make sure that the
        # expected number of comparisons are happening.
        self.assertEqual(count, 26)

    def test_recommended_ids(self):
        index = pd.Index(['c6ca034a-223f-40b4-a0e0-45942912a5ea', 'My.ID'],
                         name='id')
        df = pd.DataFrame({'col1': ['foo', 'bar']}, index=index)
        md = SampleMetadata(df)

        self.assertEqual(md.id_count, 2)
        self.assertEqual(md.column_count, 1)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids,
                         ('c6ca034a-223f-40b4-a0e0-45942912a5ea', 'My.ID'))
        self.assertEqualColumns(md.columns, [('col1', 'categorical')])

    def test_non_standard_characters(self):
        index = pd.Index(['©id##1', '((id))2', "'id_3<>'", '"id#4"',
                          'i d\r\t\n5'], name='id')
        columns = ['↩c@l1™', 'col(#2)', "#col'3", '"<col_4>"', 'col\t  \r\n5']
        data = [
            ['ƒoo', '(foo)', '#f o #o', 'fo\ro', np.nan],
            ["''2''", 'b#r', 'ba\nr', np.nan, np.nan],
            ['b"ar', 'c\td', '4\r\n2', np.nan, np.nan],
            ['b__a_z', '<42>', '>42', np.nan, np.nan],
            ['baz', np.nan, '42']
        ]
        df = pd.DataFrame(data, index=index, columns=columns)
        md = SampleMetadata(df)

        self.assertEqual(md.id_count, 5)
        self.assertEqual(md.column_count, 5)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(
            md.ids, ('©id##1', '((id))2', "'id_3<>'", '"id#4"', 'i d\r\t\n5'))
        self.assertEqualColumns(md.columns, [('↩c@l1™', 'categorical'),
                                             ('col(#2)', 'categorical'),
                                             ("#col'3", 'categorical'),
                                             ('"<col_4>"', 'categorical'),
                                             ('col\t  \r\n5', 'numeric')])

    def test_missing_data(self):
        index = pd.Index(['None', 'nan', 'NA', 'foo'], name='id')
        df = pd.DataFrame(collections.OrderedDict([
            ('col1', [1.0, np.nan, np.nan, np.nan]),
            ('NA', [np.nan, np.nan, np.nan, np.nan]),
            ('col3', ['null', 'N/A', np.nan, 'NA']),
            ('col4', np.array([np.nan, np.nan, np.nan, np.nan],
                              dtype=object))]),
            index=index)
        md = SampleMetadata(df)

        self.assertEqual(md.id_count, 4)
        self.assertEqual(md.column_count, 4)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('None', 'nan', 'NA', 'foo'))
        self.assertEqualColumns(md.columns, [('col1', 'numeric'),
                                             ('NA', 'numeric'),
                                             ('col3', 'categorical'),
                                             ('col4', 'categorical')])

    def test_missing_data_insdc(self):
        index = pd.Index(['None', 'nan', 'NA', 'foo'], name='id')
        df = pd.DataFrame(collections.OrderedDict([
            ('col1', [1.0, np.nan, 'missing', np.nan]),
            # TODO: it is not currently possible to have an ENTIRELY numeric
            # column from missing terms, as the dtype of the series is object
            # and there is not way to indicate the dtype beyond that.
            # ('NA', [np.nan, np.nan, 'not applicable', np.nan]),
            ('col3', ['null', 'N/A', 'not collected', 'NA']),
            ('col4', np.array([np.nan, np.nan, 'restricted access', np.nan],
                              dtype=object))]),
            index=index)
        md = SampleMetadata(df, default_missing_scheme='INSDC:missing')

        self.assertEqual(md.id_count, 4)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('None', 'nan', 'NA', 'foo'))
        self.assertEqualColumns(md.columns, [('col1', 'numeric'),
                                             ('col3', 'categorical'),
                                             ('col4', 'categorical')])

        pd.testing.assert_frame_equal(md.to_dataframe(), pd.DataFrame(
            {'col1': [1.0, np.nan, np.nan, np.nan],
             'col3': ['null', 'N/A', np.nan, 'NA'],
             'col4': np.array([np.nan, np.nan, np.nan, np.nan], dtype=object)},
            index=index))

    def test_missing_data_insdc_column_missing(self):
        index = pd.Index(['None', 'nan', 'NA', 'foo'], name='id')
        df = pd.DataFrame(collections.OrderedDict([
            ('col1', [1.0, np.nan, 'missing', np.nan]),
            # TODO: it is not currently possible to have an ENTIRELY numeric
            # column from missing terms, as the dtype of the series is object
            # and there is not way to indicate the dtype beyond that.
            # ('NA', [np.nan, np.nan, 'not applicable', np.nan]),
            ('col3', ['null', 'N/A', 'not collected', 'NA']),
            ('col4', np.array([np.nan, np.nan, 'restricted access', np.nan],
                              dtype=object))]),
            index=index)
        md = SampleMetadata(df, column_missing_schemes={
                              'col1': 'INSDC:missing',
                              'col3': 'INSDC:missing',
                              'col4': 'INSDC:missing'
                          })

        self.assertEqual(md.id_count, 4)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('None', 'nan', 'NA', 'foo'))
        self.assertEqualColumns(md.columns, [('col1', 'numeric'),
                                             ('col3', 'categorical'),
                                             ('col4', 'categorical')])

        pd.testing.assert_frame_equal(md.to_dataframe(), pd.DataFrame(
            {'col1': [1.0, np.nan, np.nan, np.nan],
             'col3': ['null', 'N/A', np.nan, 'NA'],
             'col4': np.array([np.nan, np.nan, np.nan, np.nan], dtype=object)},
            index=index))

    def test_missing_data_default_override(self):
        index = pd.Index(['None', 'nan', 'NA', 'foo'], name='id')
        df = pd.DataFrame(collections.OrderedDict([
            ('col1', [1.0, np.nan, 'missing', np.nan]),
            # TODO: it is not currently possible to have an ENTIRELY numeric
            # column from missing terms, as the dtype of the series is object
            # and there is not way to indicate the dtype beyond that.
            # ('NA', [np.nan, np.nan, 'not applicable', np.nan]),
            ('col3', ['null', 'N/A', 'not collected', 'NA']),
            ('col4', np.array([np.nan, np.nan, 'restricted access', np.nan],
                              dtype=object))]),
            index=index)
        md = SampleMetadata(df, column_missing_schemes={
                              'col1': 'INSDC:missing',
                              'col3': 'INSDC:missing',
                              'col4': 'INSDC:missing'
                          }, default_missing_scheme='no-missing')

        self.assertEqual(md.id_count, 4)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('None', 'nan', 'NA', 'foo'))
        self.assertEqualColumns(md.columns, [('col1', 'numeric'),
                                             ('col3', 'categorical'),
                                             ('col4', 'categorical')])

        pd.testing.assert_frame_equal(md.to_dataframe(), pd.DataFrame(
            {'col1': [1.0, np.nan, np.nan, np.nan],
             'col3': ['null', 'N/A', np.nan, 'NA'],
             'col4': np.array([np.nan, np.nan, np.nan, np.nan], dtype=object)},
            index=index))

    def test_does_not_cast_ids_or_column_names(self):
        index = pd.Index(['0.000001', '0.004000', '0.000000'], dtype=object,
                         name='id')
        columns = ['42.0', '1000', '-4.2']
        data = [
            [2.0, 'b', 2.5],
            [1.0, 'b', 4.2],
            [3.0, 'c', -9.999]
        ]
        df = pd.DataFrame(data, index=index, columns=columns)
        md = SampleMetadata(df)

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 3)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('0.000001', '0.004000', '0.000000'))
        self.assertEqualColumns(md.columns, [('42.0', 'numeric'),
                                             ('1000', 'categorical'),
                                             ('-4.2', 'numeric')])

    def test_mixed_column_types(self):
        md = SampleMetadata(
            pd.DataFrame({'col0': [1.0, 2.0, 3.0],
                          'col1': ['a', 'b', 'c'],
                          'col2': ['foo', 'bar', '42'],
                          'col3': ['1.0', '2.5', '-4.002'],
                          'col4': [1, 2, 3],
                          'col5': [1, 2, 3.5],
                          'col6': [1e-4, -0.0002, np.nan],
                          'col7': ['cat', np.nan, 'dog'],
                          'col8': ['a', 'a', 'a'],
                          'col9': [0, 0, 0]},
                         index=pd.Index(['id1', 'id2', 'id3'], name='id')))

        self.assertEqual(md.id_count, 3)
        self.assertEqual(md.column_count, 10)
        self.assertEqual(md.id_header, 'id')
        self.assertEqual(md.ids, ('id1', 'id2', 'id3'))
        self.assertEqualColumns(md.columns, [('col0', 'numeric'),
                                             ('col1', 'categorical'),
                                             ('col2', 'categorical'),
                                             ('col3', 'categorical'),
                                             ('col4', 'numeric'),
                                             ('col5', 'numeric'),
                                             ('col6', 'numeric'),
                                             ('col7', 'categorical'),
                                             ('col8', 'categorical'),
                                             ('col9', 'numeric')])

    def test_case_insensitive_duplicate_ids(self):
        index = pd.Index(['a', 'b', 'A'], name='id')
        df = pd.DataFrame({'column': ['1', '2', '3']}, index=index)
        metadata = SampleMetadata(df)

        self.assertEqual(metadata.ids, ('a', 'b', 'A'))

    def test_case_insensitive_duplicate_column_names(self):
        index = pd.Index(['a', 'b', 'c'], name='id')
        df = pd.DataFrame({'column': ['1', '2', '3'],
                           'Column': ['4', '5', '6']}, index=index)
        metadata = SampleMetadata(df)

        self.assertEqual(set(metadata.columns), {'column', 'Column'})

    def test_categorical_column_leading_trailing_whitespace_value(self):
        md1 = SampleMetadata(pd.DataFrame(
            {'col1': [1, 2, 3],
             'col2': ['foo', ' bar ', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))
        md2 = SampleMetadata(pd.DataFrame(
            {'col1': [1, 2, 3],
             'col2': ['foo', 'bar', 'baz']},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(md1, md2)

    def test_leading_trailing_whitespace_id(self):
        md1 = SampleMetadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': [4, 5, 6]},
            index=pd.Index(['a', ' b ', 'c'], name='id')))
        md2 = SampleMetadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': [4, 5, 6]},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(md1, md2)

    def test_leading_trailing_whitespace_column_name(self):
        md1 = SampleMetadata(pd.DataFrame(
            {'col1': [1, 2, 3], ' col2 ': [4, 5, 6]},
            index=pd.Index(['a', 'b', 'c'], name='id')))
        md2 = SampleMetadata(pd.DataFrame(
            {'col1': [1, 2, 3], 'col2': [4, 5, 6]},
            index=pd.Index(['a', 'b', 'c'], name='id')))

        self.assertEqual(md1, md2)


class TestRepr(unittest.TestCase):
    def test_singular(self):
        md = SampleMetadata(pd.DataFrame({'col1': [42]},
                                   index=pd.Index(['a'], name='id')))

        obs = repr(md)

        self.assertIn('Metadata', obs)
        self.assertIn('1 ID x 1 column', obs)
        self.assertIn("col1: ColumnProperties(type='numeric',"
                      " missing_scheme='blank')", obs)

    def test_plural(self):
        md = SampleMetadata(pd.DataFrame({'col1': [42, 42], 'col2': ['foo', 'bar']},
                                   index=pd.Index(['a', 'b'], name='id')))

        obs = repr(md)

        self.assertIn('Metadata', obs)
        self.assertIn('2 IDs x 2 columns', obs)
        self.assertIn("col1: ColumnProperties(type='numeric',"
                      " missing_scheme='blank')", obs)
        self.assertIn("col2: ColumnProperties(type='categorical',"
                      " missing_scheme='blank')", obs)

    def test_column_name_padding(self):
        data = [[0, 42, 'foo']]
        index = pd.Index(['my-id'], name='id')
        columns = ['col1', 'longer-column-name', 'c']
        md = SampleMetadata(pd.DataFrame(data, index=index, columns=columns))

        obs = repr(md)

        self.assertIn('Metadata', obs)
        self.assertIn('1 ID x 3 columns', obs)
        self.assertIn(
            "col1:               ColumnProperties(type='numeric',"
            " missing_scheme='blank')", obs)
        self.assertIn(
            "longer-column-name: ColumnProperties(type='numeric',"
            " missing_scheme='blank')", obs)
        self.assertIn(
            "c:                  ColumnProperties(type='categorical',"
            " missing_scheme='blank')", obs)


class TestToDataframe(unittest.TestCase):
    def test_minimal(self):
        df = pd.DataFrame({}, index=pd.Index(['id1'], name='id'))
        md = SampleMetadata(df)

        obs = md.to_dataframe()

        pd.testing.assert_frame_equal(obs, df)

    def test_id_header_preserved(self):
        df = pd.DataFrame({'col1': [42, 2.5], 'col2': ['foo', 'bar']},
                          index=pd.Index(['id1', 'id2'], name='#SampleID'))
        md = SampleMetadata(df)

        obs = md.to_dataframe()

        pd.testing.assert_frame_equal(obs, df)
        self.assertEqual(obs.index.name, '#SampleID')

    def test_dataframe_copy(self):
        df = pd.DataFrame({'col1': [42, 2.5], 'col2': ['foo', 'bar']},
                          index=pd.Index(['id1', 'id2'], name='id'))
        md = SampleMetadata(df)

        obs = md.to_dataframe()

        pd.testing.assert_frame_equal(obs, df)
        self.assertIsNot(obs, df)

    def test_retains_column_order(self):
        index = pd.Index(['id1', 'id2'], name='id')
        columns = ['z', 'a', 'ch']
        data = [
            [1.0, 'a', 'foo'],
            [2.0, 'b', 'bar']
        ]
        df = pd.DataFrame(data, index=index, columns=columns)
        md = SampleMetadata(df)

        obs = md.to_dataframe()

        pd.testing.assert_frame_equal(obs, df)
        self.assertEqual(obs.columns.tolist(), ['z', 'a', 'ch'])

    def test_missing_data(self):
        # Different missing data representations should be normalized to np.nan
        index = pd.Index(['None', 'nan', 'NA', 'id1'], name='id')
        df = pd.DataFrame(collections.OrderedDict([
            ('col1', [42.5, np.nan, float('nan'), 3]),
            ('NA', [np.nan, 'foo', float('nan'), None]),
            ('col3', ['null', 'N/A', np.nan, 'NA']),
            ('col4', np.array([np.nan, np.nan, np.nan, np.nan],
                              dtype=object))]),
            index=index)
        md = SampleMetadata(df)

        obs = md.to_dataframe()

        exp = pd.DataFrame(collections.OrderedDict([
            ('col1', [42.5, np.nan, np.nan, 3.0]),
            ('NA', [np.nan, 'foo', np.nan, np.nan]),
            ('col3', ['null', 'N/A', np.nan, 'NA']),
            ('col4', np.array([np.nan, np.nan, np.nan, np.nan],
                              dtype=object))]),
            index=index)

        pd.testing.assert_frame_equal(obs, exp)
        self.assertEqual(obs.dtypes.to_dict(),
                         {'col1': np.float64, 'NA': object, 'col3': object,
                          'col4': object})
        self.assertTrue(np.isnan(obs['col1']['NA']))
        self.assertTrue(np.isnan(obs['NA']['NA']))
        self.assertTrue(np.isnan(obs['NA']['id1']))

    def test_dtype_int_normalized_to_dtype_float(self):
        index = pd.Index(['id1', 'id2', 'id3'], name='id')
        df = pd.DataFrame({'col1': [42, -43, 0],
                           'col2': [42.0, -43.0, 0.0],
                           'col3': [42, np.nan, 0]},
                          index=index)

        self.assertEqual(df.dtypes.to_dict(),
                         {'col1': np.int64, 'col2': np.float64,
                          'col3': np.float64})

        md = SampleMetadata(df)
        obs = md.to_dataframe()

        exp = pd.DataFrame({'col1': [42.0, -43.0, 0.0],
                            'col2': [42.0, -43.0, 0.0],
                            'col3': [42.0, np.nan, 0.0]},
                           index=index)

        pd.testing.assert_frame_equal(obs, exp)
        self.assertEqual(obs.dtypes.to_dict(),
                         {'col1': np.float64, 'col2': np.float64,
                          'col3': np.float64})

    def test_encode_missing_no_missing(self):
        df = pd.DataFrame({'col1': [42.0, 50.0],
                           'col2': ['foo', 'bar']},
                          index=pd.Index(['id1', 'id2'], name='id'))
        md = SampleMetadata(df, default_missing_scheme='INSDC:missing')

        obs = md.to_dataframe(encode_missing=True)

        pd.testing.assert_frame_equal(obs, df)
        self.assertIsNot(obs, df)

    def test_insdc_missing_encode_missing_true(self):
        df = pd.DataFrame({'col1': [42, 'missing'],
                           'col2': ['foo', 'not applicable']},
                          index=pd.Index(['id1', 'id2'], name='id'))
        md = SampleMetadata(df, default_missing_scheme='INSDC:missing')

        obs = md.to_dataframe(encode_missing=True)

        pd.testing.assert_frame_equal(obs, df)
        self.assertIsNot(obs, df)

    def test_insdc_missing_encode_missing_false(self):
        df = pd.DataFrame({'col1': [42, 'missing'],
                           'col2': ['foo', 'not applicable']},
                          index=pd.Index(['id1', 'id2'], name='id'))
        md = SampleMetadata(df, default_missing_scheme='INSDC:missing')

        obs = md.to_dataframe()

        exp = pd.DataFrame({'col1': [42, np.nan],
                            'col2': ['foo', np.nan]},
                           index=pd.Index(['id1', 'id2'], name='id'))

        pd.testing.assert_frame_equal(obs, exp)
        self.assertIsNot(obs, df)


class TestGetIDs(unittest.TestCase):
    def test_default(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = SampleMetadata(df)

        actual = metadata.get_ids()
        expected = {'S1', 'S2', 'S3'}
        self.assertEqual(actual, expected)

    def test_incomplete_where(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='sampleid'))
        metadata = SampleMetadata(df)

        where = "Subject='subject-1' AND SampleType="
        with self.assertRaises(ValueError):
            metadata.get_ids(where)

        where = "Subject="
        with self.assertRaises(ValueError):
            metadata.get_ids(where)

    def test_invalid_where(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='sampleid'))
        metadata = SampleMetadata(df)

        where = "not-a-column-name='subject-1'"
        with self.assertRaises(ValueError):
            metadata.get_ids(where)

    def test_empty_result(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = SampleMetadata(df)

        where = "Subject='subject-3'"
        actual = metadata.get_ids(where)
        expected = set()
        self.assertEqual(actual, expected)

    def test_simple_expression(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = SampleMetadata(df)

        where = "Subject='subject-1'"
        actual = metadata.get_ids(where)
        expected = {'S1', 'S2'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-2'"
        actual = metadata.get_ids(where)
        expected = {'S3'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-3'"
        actual = metadata.get_ids(where)
        expected = set()
        self.assertEqual(actual, expected)

        where = "SampleType='gut'"
        actual = metadata.get_ids(where)
        expected = {'S1', 'S3'}
        self.assertEqual(actual, expected)

        where = "SampleType='tongue'"
        actual = metadata.get_ids(where)
        expected = {'S2'}
        self.assertEqual(actual, expected)

    def test_more_complex_expressions(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = SampleMetadata(df)

        where = "Subject='subject-1' OR Subject='subject-2'"
        actual = metadata.get_ids(where)
        expected = {'S1', 'S2', 'S3'}
        self.assertEqual(actual, expected)

        where = "Subject='subject-1' AND Subject='subject-2'"
        actual = metadata.get_ids(where)
        expected = set()
        self.assertEqual(actual, expected)

        where = "Subject='subject-1' AND SampleType='gut'"
        actual = metadata.get_ids(where)
        expected = {'S1'}
        self.assertEqual(actual, expected)

    def test_query_by_id(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = SampleMetadata(df)

        actual = metadata.get_ids(where="id='S2' OR id='S1'")
        expected = {'S1', 'S2'}
        self.assertEqual(actual, expected)

    def test_query_by_alternate_id_header(self):
        metadata = SampleMetadata(pd.DataFrame(
            {}, index=pd.Index(['id1', 'id2', 'id3'], name='#OTU ID')))

        obs = metadata.get_ids(where="\"#OTU ID\" IN ('id2', 'id3')")

        exp = {'id2', 'id3'}
        self.assertEqual(obs, exp)

    def test_no_columns(self):
        metadata = SampleMetadata(
            pd.DataFrame({}, index=pd.Index(['a', 'b', 'my-id'], name='id')))

        obs = metadata.get_ids()

        exp = {'a', 'b', 'my-id'}
        self.assertEqual(obs, exp)

    def test_query_mixed_column_types(self):
        df = pd.DataFrame({'Name': ['Foo', 'Bar', 'Baz', 'Baaz'],
                           # numbers that would sort incorrectly as strings
                           'Age': [9, 10, 11, 101],
                           'Age_Str': ['9', '10', '11', '101'],
                           'Weight': [80.5, 85.3, np.nan, 120.0]},
                          index=pd.Index(['S1', 'S2', 'S3', 'S4'], name='id'))
        metadata = SampleMetadata(df)

        # string pattern matching
        obs = metadata.get_ids(where="Name LIKE 'Ba_'")
        exp = {'S2', 'S3'}
        self.assertEqual(obs, exp)

        # string comparison
        obs = metadata.get_ids(where="Age_Str >= 11")
        exp = {'S1', 'S3'}
        self.assertEqual(obs, exp)

        # numeric comparison
        obs = metadata.get_ids(where="Age >= 11")
        exp = {'S3', 'S4'}
        self.assertEqual(obs, exp)

        # numeric comparison with missing data
        obs = metadata.get_ids(where="Weight < 100")
        exp = {'S1', 'S2'}
        self.assertEqual(obs, exp)

    def test_column_with_space_in_name(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'Sample Type': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = SampleMetadata(df)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            metadata.get_ids()
            # The list of captured warnings should be empty
            self.assertFalse(w)


if __name__ == '__main__':
    unittest.main()
