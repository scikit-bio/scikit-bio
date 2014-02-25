#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from bipy.util.unit_test import TestCase, main
from bipy.app.parameters import (FlagParameter, ValuedParameter,
                                 MixedParameter, Parameters,
                                 ParameterError, FilePath)


class FlagParameterTests(TestCase):

    """ Tests of the FlagParameter class """

    def setUp(self):
        """Setup some variables for the tests to use """
        self.p_modify_prefix = [FlagParameter(Name='d', Prefix='-'),
                                FlagParameter(Name='d', Prefix='--'),
                                FlagParameter(Name='d', Prefix='')]

        self.p_modify_name = [FlagParameter(Name='d', Prefix='-'),
                              FlagParameter(Name='D', Prefix='-'),
                              FlagParameter(Name=4, Prefix='-'),
                              FlagParameter(Name='abcdef', Prefix='-')]

        self.p_On = [FlagParameter(Name='d', Prefix='-', Value=True),
                     FlagParameter(Name='d', Prefix='-', Value=5),
                     FlagParameter(Name='d', Prefix='-', Value=[1]),
                     FlagParameter(Name='d', Prefix='-', Value='F')]

        self.p_Off = [FlagParameter(Name='d', Prefix='-', Value=False),
                      FlagParameter(Name='d', Prefix='-', Value=None),
                      FlagParameter(Name='d', Prefix='-', Value=[]),
                      FlagParameter(Name='d', Prefix='-', Value=0),
                      FlagParameter(Name='d', Prefix='-', Value='')]

        self.ID_tests = [FlagParameter(Name='d', Prefix='-'),
                         FlagParameter(Name='d', Prefix=''),
                         FlagParameter(Name='', Prefix='-'),
                         FlagParameter(Name=4, Prefix='-'),
                         FlagParameter(Name=None, Prefix='-'),
                         FlagParameter(Name=4, Prefix=None),
                         FlagParameter(Name='abcdef', Prefix='-')]

    def test_init(self):
        """FlagParameter: init functions as expected """
        param = FlagParameter(Name='a', Prefix='-', Value=42)
        self.assertEqual(param.Name, 'a')
        self.assertEqual(param.Prefix, '-')
        self.assertEqual(param.Value, 42)
        self.assertEqual(param.Delimiter, None)
        self.assertEqual(param.Quote, None)
        self.assertEqual(param.Id, '-a')

    def test_init_defaults(self):
        """FlagParameter: init functions as expected with default values"""
        p = FlagParameter(Name='a', Prefix='-')
        self.assertEqual(p.Name, 'a')
        self.assertEqual(p.Prefix, '-')
        self.assertEqual(p.Value, False)
        self.assertEqual(p.Delimiter, None)
        self.assertEqual(p.Quote, None)
        self.assertEqual(p.Id, '-a')

    def test_get_id(self):
        """FlagParameter: _get_id functions as expected """

        expected_results = ['-d', 'd', '-', '-4', '-', '4', '-abcdef']

        for param, exp in zip(self.ID_tests, expected_results):
            self.assertEqual(param._get_id(), exp)

    def test_eq(self):
        """FlagParameter: eq functions as expected """
        p1 = FlagParameter(Name='a', Prefix='-', Value=True)
        p2 = FlagParameter(Name='a', Prefix='-', Value=True)
        p3 = FlagParameter(Name='a', Prefix='-')
        p4 = FlagParameter(Name='i', Prefix='-', Value=True)
        p5 = FlagParameter(Name='a', Prefix='--', Value=True)

        assert p1 == p2
        assert not p1 == p3
        assert not p1 == p4
        assert not p1 == p5
        assert not p3 == p4
        assert not p3 == p5
        assert not p4 == p5

    def test_ne(self):
        """FlagParameter: ne functions as expected """
        p1 = FlagParameter(Name='a', Prefix='-', Value=True)
        p2 = FlagParameter(Name='a', Prefix='-', Value=True)
        p3 = FlagParameter(Name='a', Prefix='-')
        p4 = FlagParameter(Name='i', Prefix='-', Value=True)
        p5 = FlagParameter(Name='a', Prefix='--', Value=True)

        assert not p1 != p2
        assert p1 != p3
        assert p1 != p4
        assert p1 != p5
        assert p3 != p4
        assert p3 != p5
        assert p4 != p5

    def test_isOn_True(self):
        """FlagParameter: isOn functions as expected with True Values """
        for param in self.p_On:
            assert param.isOn()

    def test_isOn_False(self):
        """FlagParameter: isOn functions as expected with False Values """
        for param in self.p_Off:
            assert not param.isOn()

    def test_isOff_True(self):
        """FlagParameter: isOff functions as expected with True values """
        for param in self.p_Off:
            assert param.isOff()

    def test_isOff_False(self):
        """FlagParameter: isOff functions as expected with False values """
        for param in self.p_On:
            assert not param.isOff()

    def test_on(self):
        """FlagParameter: on functions as expected """
        for param in self.p_On + self.p_Off:
            param.on()
            assert param.isOn()

    def test_off(self):
        """FlagParameter: off functions as expected """
        for param in self.p_On + self.p_Off:
            param.off()
            assert param.isOff()

    def test_str_modify_prefix(self):
        """FlagParameter: str functions as expected with different prefixes """

        expected_results = ['-d', '--d', 'd']

        for param, exp in zip(self.p_modify_prefix, expected_results):
            param.on()
            self.assertEqual(str(param), exp)

    def test_str_modify_name(self):
        """FlagParameter: str functions as expected with different names """

        expected_results = ['-d', '-D', '-4', '-abcdef']

        for param, exp in zip(self.p_modify_name, expected_results):
            param.on()
            self.assertEqual(str(param), exp)


class ValuedParameterTests(TestCase):

    """ Tests of the ValuedParameter class """
    constructor = ValuedParameter
    s = 'Valued'

    def setUp(self):
        """Setup some variables for the tests to use """
        self.p_modify_prefix = [self.constructor(Name='d', Prefix='-'),
                                self.constructor(Name='d', Prefix='--'),
                                self.constructor(Name='d', Prefix='')]

        self.p_modify_name = [self.constructor(Name='d', Prefix='-'),
                              self.constructor(Name='D', Prefix='-'),
                              self.constructor(Name=4, Prefix='-'),
                              self.constructor(Name='abcdef', Prefix='-')]

        self.p_On = [self.constructor(Name='d', Prefix='-', Value=True),
                     self.constructor(Name='d', Prefix='-', Value=5),
                     self.constructor(Name='d', Prefix='-', Value=[1]),
                     self.constructor(Name='d', Prefix='-', Value=False),
                     self.constructor(Name='d', Prefix='-', Value='F')]

        self.p_Off = [self.constructor(Name='d', Prefix='-', Value=None)]

        self.p_full = [self.constructor(Name='a', Prefix='-',
                                        Value=42, Delimiter=' ', Quote='\'')]

        self.p_default = [self.constructor(Name='a', Prefix='-')]

        self.p_modified_prefix = [self.constructor(Name='d', Prefix='-'),
                                  self.constructor(Name='d', Prefix='--'),
                                  self.constructor(Name='d', Prefix='')]

        self.p_modified_name = [self.constructor(Name='d', Prefix='-'),
                                self.constructor(Name='D', Prefix='-'),
                                self.constructor(Name=4, Prefix='-'),
                                self.constructor(Name='abcdef', Prefix='-')]

        self.p_modified_delimiter =\
            [self.constructor(Name='d', Prefix='-', Value=42),
             self.constructor(Name='d', Prefix='-', Value=42, Delimiter=''),
             self.constructor(Name='d', Prefix='-', Value=42, Delimiter=' '),
             self.constructor(Name='d', Prefix='-', Value=42, Delimiter=9),
             self.constructor(Name='d', Prefix='-', Value=42, Delimiter='=')]

        self.p_modified_value =\
            [self.constructor(Name='d', Prefix='-', Value=42, Delimiter=' '),
             self.constructor(
                 Name='d',
                 Prefix='-',
                 Value='pbl',
                 Delimiter=' '),
                self.constructor(
                    Name='d',
                    Prefix='-',
                    Value='2-2',
                    Delimiter=' '),
                self.constructor(Name='d', Prefix='-', Value='evo/t.txt',
                                 Delimiter=' '),
                self.constructor(Name='d', Prefix='-', Value='\'',
                                 Delimiter=' ')]

        self.p_modified_quote =\
            [self.constructor(Name='d', Prefix='-', Value=42, Quote=''),
             self.constructor(Name='d', Prefix='-', Value=42),
             self.constructor(Name='d', Prefix='-', Value=42, Quote=' '),
             self.constructor(Name='d', Prefix='-', Value=42, Quote='\''),
             self.constructor(Name='d', Prefix='-', Value=42, Quote='\"'),
             self.constructor(Name='d', Prefix='-', Value=42, Quote='x')]

        self.ID_tests = [self.constructor(Name='d', Prefix='-'),
                         self.constructor(Name='d', Prefix=''),
                         self.constructor(Name='', Prefix='-'),
                         self.constructor(Name=4, Prefix='-'),
                         self.constructor(Name=None, Prefix='-'),
                         self.constructor(Name=4, Prefix=None),
                         self.constructor(Name='abcdef', Prefix='-')]

        self.p_modified_is_path =\
            [self.constructor(Name='d', Prefix='-', Delimiter=' ',
                              Value='test.txt', IsPath=True),
             self.constructor(Name='d', Prefix='-', Delimiter=' ',
                              Value='test.txt', IsPath=False),
             self.constructor(Name='d', Prefix='-', Delimiter=' ',
                              Value='test.txt', Quote='"', IsPath=True)]

    def test_init(self):
        """Parameter: init functions as expected """

        for param in self.p_full:
            self.assertEqual(param.Name, 'a')
            self.assertEqual(param.Prefix, '-')
            self.assertEqual(param.Value, 42)
            self.assertEqual(param.Delimiter, ' ')
            self.assertEqual(param.Quote, '\'')
            self.assertEqual(param.Id, '-a')

    def test_init_defaults(self):
        """Parameter: init functions as expected with default values"""
        for p in self.p_default:
            self.assertEqual(p.Name, 'a')
            self.assertEqual(p.Prefix, '-')
            self.assertEqual(p.Value, None)
            self.assertEqual(p.Delimiter, None)
            self.assertEqual(p.Quote, None)
            self.assertEqual(p.Id, '-a')

    def test_get_id(self):
        """Parameter: _get_id functions as expected """

        expected_results = ['-d', 'd', '-', '-4', '-', '4', '-abcdef']

        for param, exp in zip(self.ID_tests, expected_results):
            self.assertEqual(param._get_id(), exp)

    def test_eq(self):
        """Parameter: eq functions as expected """
        p1 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=')
        p2 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=')
        p3 = self.constructor(Name='dsf', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=')
        p4 = self.constructor(Name='a', Prefix='--', Value=42, Quote='\'',
                              Delimiter='=')
        p5 = self.constructor(Name='a', Prefix='-', Value=942, Quote='\'',
                              Delimiter='=')
        p6 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\"',
                              Delimiter='=')
        p7 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='!!!')
        p8 = self.constructor(Name='wwwww', Prefix='-------')
        p9 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=', IsPath=True)

        assert p1 == p2
        assert not p1 == p3
        assert not p1 == p4
        assert not p1 == p5
        assert not p1 == p6
        assert not p1 == p7
        assert not p1 == p8
        assert not p1 == p9
        # test default setting
        p5.Value = 42
        assert not p1 == p5

    def test_ne(self):
        """Parameter: ne functions as expected """
        p1 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=')
        p2 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=')
        p3 = self.constructor(Name='dsf', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=')
        p4 = self.constructor(Name='a', Prefix='--', Value=42, Quote='\'',
                              Delimiter='=')
        p5 = self.constructor(Name='a', Prefix='-', Value=942, Quote='\'',
                              Delimiter='=')
        p6 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\"',
                              Delimiter='=')
        p7 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='!!!')
        p8 = self.constructor(Name='wwwww', Prefix='-------')
        p9 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=', IsPath=True)

        assert not p1 != p2
        assert p1 != p3
        assert p1 != p4
        assert p1 != p5
        assert p1 != p6
        assert p1 != p7
        assert p1 != p8
        assert p1 != p9
        # test default setting
        p5.Value = 42
        assert p1 != p5

    def test_get_default(self):
        """Parameter: default behaves as expected """
        p1 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=')
        self.assertEqual(p1._get_default(), 42)
        p1.Value = 43
        self.assertEqual(p1._get_default(), 42)

    def test_get_default_w_IsPath(self):
        """Parameter: default is a FilePath object when IsPath is set """
        p = self.constructor(
            Name='a', Prefix='-', Value='test.txt', Quote='\'',
            Delimiter='=', IsPath=True)
        self.assertEqual(p._get_default(), 'test.txt')
        self.assertEqual(p.Default, 'test.txt')
        p.Value = 'test2.txt'
        self.assertEqual(p._get_default(), 'test.txt')
        self.assertEqual(p.Default, 'test.txt')
        assert isinstance(p._get_default(), FilePath)
        assert isinstance(p.Default, FilePath)

    def test_reset(self):
        """Parameter: reset correctly set Value to _default """
        p1 = self.constructor(Name='a', Prefix='-', Value=42, Quote='\'',
                              Delimiter='=')
        p1.Value = 43
        self.assertNotEqual(p1.Default, p1.Value)
        p1.reset()
        self.assertEqual(p1.Default, p1.Value)

    def test_isOn_True(self):
        """Parameter: isOn functions as expected with True Values """
        for param in self.p_On:
            assert param.isOn()

    def test_isOn_False(self):
        """Parameter: isOn functions as expected with False Values """
        for param in self.p_Off:
            assert not param.isOn()

    def test_isOff_True(self):
        """Parameter: isOff functions as expected with True values """
        for param in self.p_Off:
            assert param.isOff()

    def test_isOff_False(self):
        """Parameter: isOff functions as expected with False values """
        for param in self.p_On:
            assert not param.isOff()

    def test_on(self):
        """Parameter: on functions as expected """
        for param in self.p_On + self.p_Off:
            param.on('a')
            assert param.isOn()
        p = self.p_On[0]
        self.assertRaises(ParameterError, p.on, None)

    def test_off(self):
        """Parameter: off functions as expected """
        for param in self.p_On + self.p_Off:
            param.off()
            assert param.isOff()

    def test_str_off(self):
        """Parameter: str() prints empty string when off """
        for p in self.p_Off:
            self.assertEqual(str(p), '')

    def test_str_modify_prefix(self):
        """Parameter: str functions as expected with different prefixes """

        expected_results = ['-d', '--d', 'd']

        for param, exp in zip(self.p_modified_prefix, expected_results):
            param.on('')
            self.assertEqual(str(param), exp)

    def test_str_modify_name(self):
        """Parameter: str functions as expected with different names """

        expected_results = ['-d', '-D', '-4', '-abcdef']

        for param, exp in zip(self.p_modified_name, expected_results):
            param.on('')
            self.assertEqual(str(param), exp)

    def test_str_modify_delimiter(self):
        """Parameter: str functions as expected with different delimiter """

        expected_results = ['-d42', '-d42', '-d 42', '-d942', '-d=42']

        for param, exp in zip(self.p_modified_delimiter, expected_results):
            self.assertEqual(str(param), exp)

    def test_str_modify_values(self):
        """Parameter: str functions as expected with different values """

        expected_results = ['-d 42',
                            '-d pbl', '-d 2-2', '-d evo/t.txt', '-d \'']

        for param, exp in zip(self.p_modified_value, expected_results):
            self.assertEqual(str(param), exp)

    def test_str_modify_quotes(self):
        """Parameter: str functions as expected with different quotes """

        expected_results = ['-d42', '-d42', '-d 42 ', '-d\'42\'',
                            '-d\"42\"', '-dx42x']

        for param, exp in zip(self.p_modified_quote, expected_results):
            self.assertEqual(str(param), exp)

    def test_str_modify_is_path(self):
        """Parameter: str functions as expected with different IsPath """

        expected_results = ['-d "test.txt"', '-d test.txt', '-d "test.txt"']

        for param, exp in zip(self.p_modified_is_path, expected_results):
            self.assertEqual(str(param), exp)

    def test_str_full(self):
        """Parameter: str functions as expected with all values non-default """
        for p in self.p_full:
            self.assertEqual(str(p), '-a \'42\'')


class MixedParameterTests(ValuedParameterTests):

    """ Tests of the MixedParameter class """

    constructor = MixedParameter

    def setUp(self):
        """Setup some variables for the tests to use """
        super(MixedParameterTests, self).setUp()
        self.p_On = [self.constructor(Name='d', Prefix='-', Value=True),
                     self.constructor(Name='d', Prefix='-', Value=5),
                     self.constructor(Name='d', Prefix='-', Value=[1]),
                     self.constructor(Name='d', Prefix='-', Value=None),
                     self.constructor(Name='d', Prefix='-', Value='F')]

        self.p_Off = [self.constructor(Name='d', Prefix='-', Value=False)]

        # This is different from the superclass variable b/c we need to make
        # sure that specifying IsPath with Value=None functions as expected
        self.p_modified_is_path =\
            [self.constructor(Name='d', Prefix='-', Delimiter=' ',
                              Value='test.txt', IsPath=True),
             self.constructor(Name='d', Prefix='-', Delimiter=' ',
                              Value='test.txt', Quote='"', IsPath=True),
             self.constructor(Name='d', Prefix='-', Delimiter=' ',
                              Value='test.txt', IsPath=False),
             self.constructor(Name='d', Prefix='-', Delimiter=' ',
                              Value=None, IsPath=True),
             self.constructor(Name='d', Prefix='-', Delimiter=' ',
                              Value=None, IsPath=False)]

    def test_on(self):
        """Parameter: on functions as expected """
        for param in self.p_On + self.p_Off:
            param.on('a')
            assert param.isOn()
        p = self.p_On[0]
        self.assertRaises(ParameterError, p.on, False)

    def test_init_defaults(self):
        """MixedParameter: init functions as expected with default values"""
        for p in self.p_default:
            self.assertEqual(p.Name, 'a')
            self.assertEqual(p.Prefix, '-')
            self.assertEqual(p.Value, False)
            self.assertEqual(p.Delimiter, None)
            self.assertEqual(p.Quote, None)
            self.assertEqual(p.Id, '-a')
            self.assertEqual(p.IsPath, False)

    def test_str_all_modes(self):
        """MixedParameter: str() functions in various modes """
        p = MixedParameter(Prefix='-', Name='d', Delimiter='=', Quote=']')
        self.assertEqual(str(p), '')
        p.on()
        self.assertEqual(str(p), '-d')
        p.on('a')
        self.assertEqual(str(p), '-d=]a]')

    def test_str_modify_is_path(self):
        """MixedParameter: str functions as expected with different IsPath """

        # This is different from the superclass test b/c we need to make
        # sure that specifying IsPath with Value=None functions as expected
        expected_results = ['-d "test.txt"', '-d "test.txt"',
                            '-d test.txt', '-d', '-d']

        for param, exp in zip(self.p_modified_is_path, expected_results):
            self.assertEqual(str(param), exp)


class ParametersTests(TestCase):

    """Tests of the Parameters class"""

    def setUp(self):
        self.fp = FlagParameter(Prefix='-', Name='d')
        self.vp = ValuedParameter(Name='p', Prefix='-', Value=[1])
        self.mp = MixedParameter(Prefix='--', Name='k', Delimiter=' ')
        self.all_params = {self.fp.Id: self.fp, self.vp.Id: self.vp,
                           self.mp.Id: self.mp}
        self.p1 = Parameters()
        self.p2 = Parameters(self.all_params)
        self._synonyms = {'Pino': '-p', 'K': 'k'}
        self.p3 = Parameters(self.all_params, self._synonyms)

    def test_init(self):
        """Parameters: init functions as expected"""
        self.assertEqualItems(self.p1, {})
        self.assertEqualItems(self.p2, self.all_params)
        self.assertEqualItems(self.p3, self.all_params)

    def test_lookup(self):
        """Parameters: test ability to lookup """
        self.assertEqual(self.p2['-p'], self.vp)
        self.assertEqual(self.p3['Pino'], self.vp)

    def test_immutability(self):
        """Parameters: attempt to modify object raises error """
        try:
            self.p2['-p'] = 42
        except TypeError:
            pass
        else:
            raise AttributeError("Parameters shouldn't support assignment.")

        try:
            del self.p2['-p']
        except TypeError:
            pass
        else:
            raise AttributeError("Parameters shouldn't support deletion.")

    def test_all_off(self):
        """Parameters: all_off() should turn all parameters off"""
        p = self.p2
        # turn everything on
        for v in p.values():
            try:
                v.on(3)
            except TypeError:
                v.on()
            self.assertTrue(v.isOn())

        # turn everything off
        p.all_off()
        for v in p.values():
            self.assertTrue(v.isOff())


class FilePathTests(TestCase):

    """ Tests of the FilePath class """

    def setUp(self):
        """ Initialize variables to be used by tests """
        self.filename = 'filename.txt'
        self.relative_dir_path = 'a/relative/path/'
        self.relative_dir_path_no_trailing_slash = 'a/relative/path'
        self.relative_file_path = 'a/relative/filepath.txt'
        self.absolute_dir_path = '/absolute/path/'
        self.absolute_file_path = '/absolute/filepath.txt'
        self.all_paths = [self.filename, self.relative_dir_path,
                          self.relative_file_path, self.absolute_dir_path,
                          self.absolute_file_path]

    def test_init(self):
        """FilePath: initialization returns w/o error """
        for p in self.all_paths:
            self.assertEqual(FilePath(p), p)
        self.assertEqual(FilePath(''), '')

    def test_str(self):
        """FilePath: str wraps path in quotes """
        # Do one explicit test (for sanity), then automatically run
        # through the examples
        self.assertEqual(str(FilePath(self.filename)), '"filename.txt"')
        for p in self.all_paths:
            self.assertEqual(str(FilePath(p)), '"' + p + '"')

    def test_str_path_is_None(self):
        """FilePath: str return empty string when path is None """
        self.assertEqual(str(FilePath(None)), '')

    def test_add(self):
        """FilePath: add (or joining of paths) functions as expected """
        actual = FilePath(self.relative_dir_path) + FilePath(self.filename)
        expected = FilePath('a/relative/path/filename.txt')
        self.assertEqual(actual, expected)
        # result is a FilePath
        assert isinstance(actual, FilePath)
        # appending a string to a FilePath results in a FilePath
        actual = FilePath(self.relative_dir_path) + 'filename.txt'
        expected = FilePath('a/relative/path/filename.txt')
        self.assertEqual(actual, expected)
        # result is a FilePath
        assert isinstance(actual, FilePath)

    def test_FilePath_identity_preserved(self):
        """FilePath: trivial actions on FilePaths yeild original FilePath
        """
        p = FilePath(self.filename)
        # Creating FilePath from FilePath results in FilePath
        # equal to original
        self.assertEqual(FilePath(p), p)
        for p in self.all_paths:
            self.assertEqual(FilePath(p), p)
        # Appending an empty FilePath to a FilePath results in FilePath
        # equal to original
        self.assertEqual(p + FilePath(''), p)


if __name__ == '__main__':
    main()
