#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from tempfile import gettempdir
from os import remove, system, mkdir, rmdir, removedirs, getcwd, walk

from bipy.app.parameters import *
from bipy.util.unit_test import TestCase, main
from bipy.app.util import (Application, CommandLineApplication,
                           CommandLineAppResult, ResultPath, ApplicationError, ParameterIterBase,
                           ParameterCombinations, cmdline_generator, ApplicationNotFoundError,
                           get_tmp_filename, guess_input_handler, app_path)


class ParameterCombinationsTests(TestCase):

    def setUp(self):
        """Setup for ParameterCombinations tests"""
        self.mock_app = ParameterCombinationsApp
        self.params = {'-flag1': True,
                       '--value1': range(0, 5),
                       '-delim': range(0, 2),
                       '-mix1': [None] + range(0, 3)}
        self.always_on = ['--value1']
        self.param_iter = ParameterCombinations(self.mock_app, self.params,
                                                self.always_on)

    def test_init_generator(self):
        """Tests generator capabilities"""
        all_params = list(self.param_iter)
        self.assertEqual(len(all_params), 150)
        params = {'-flag1': True,
                  '--value1': 1,
                  '-delim': ['choice1', 'choice2']}
        always_on = ['-flag1', '-delim']
        param_iter = ParameterCombinations(self.mock_app, params, always_on)

        exp = [self.mock_app._parameters.copy(),
               self.mock_app._parameters.copy(),
               self.mock_app._parameters.copy(),
               self.mock_app._parameters.copy()]

        # default is on in all these cases
        exp[0]['-flag1'].on()
        exp[0]['--value1'].on(1)
        exp[0]['-delim'].on('choice1')

        exp[1]['-flag1'].on()
        exp[1]['--value1'].on(1)
        exp[1]['-delim'].on('choice2')

        exp[2]['-flag1'].on()
        exp[2]['--value1'].off()
        exp[2]['-delim'].on('choice1')

        exp[3]['-flag1'].on()
        exp[3]['--value1'].off()
        exp[3]['-delim'].on('choice2')

        obs = list(param_iter)
        self.assertEqual(obs, exp)

    def test_reset(self):
        """Resets the iterator"""
        first = list(self.param_iter)
        self.assertRaises(StopIteration, self.param_iter.next)
        self.param_iter.reset()
        second = list(self.param_iter)
        self.assertEqual(first, second)


class ParameterIterBaseTests(TestCase):

    def setUp(self):
        """Setup for ParameterIterBase tests"""
        self.mock_app = ParameterCombinationsApp
        self.params = {'-flag1': True,
                       '--value1': range(0, 5),
                       '-delim': range(0, 2),
                       '-mix1': [None] + range(0, 3)}
        self.always_on = ['--value1']
        self.param_base = ParameterIterBase(self.mock_app, self.params,
                                            self.always_on)

    def test_init(self):
        """Test constructor"""
        exp_params = {'-flag1': [True, False],
                      '--value1': range(0, 5),
                      '-delim': range(0, 2) + [False],
                      '-mix1': [None, 0, 1, 2] + [False]}
        exp_keys = exp_params.keys()
        exp_values = exp_params.values()

        self.assertEqual(sorted(self.param_base._keys), sorted(exp_keys))
        self.assertEqual(sorted(self.param_base._values), sorted(exp_values))

        self.params['asdasda'] = 5
        self.assertRaises(ValueError, ParameterIterBase, self.mock_app,
                          self.params, self.always_on)

        self.params.pop('asdasda')
        self.always_on.append('asdasd')
        self.assertRaises(ValueError, ParameterIterBase, self.mock_app,
                          self.params, self.always_on)

    def test_make_app_params(self):
        """Returns app parameters with expected values set"""
        values = [0, 0, True, None]
        exp = self.mock_app._parameters.copy()
        exp['-flag1'].on()
        exp['--value1'].on(0)
        exp['-delim'].on(0)
        exp['-mix1'].on(None)
        obs = self.param_base._make_app_params(values)
        self.assertEqual(obs, exp)

        state = [4, False, False, False]
        exp = self.mock_app._parameters.copy()
        exp['-flag1'].off()
        exp['--value1'].on(4)
        exp['-delim'].off()
        exp['-mix1'].off()
        obs = self.param_base._make_app_params(values)
        self.assertEqual(obs, exp)


class CommandLineGeneratorTests(TestCase):

    def setUp(self):
        self.abs_path_to_bin = '/bin/path'
        self.abs_path_to_cmd = '/cmd/path'
        self.abs_path_to_input = '/input/path'
        self.abs_path_to_output = '/output/path'
        self.abs_path_to_stdout = '/stdout/path'
        self.abs_path_to_stderr = '/stderr/path'
        self.app = ParameterCombinationsApp

        params = {'-flag1': True,
                  '-delim': ['choice1', 'choice2']}
        always_on = ['-delim']
        self.mock_app = ParameterCombinationsApp
        self.param_iter = ParameterCombinations(
            self.mock_app,
            params,
            always_on)

    def test_cmdline_generator_easy(self):
        """Returns parameter combinations commandlines"""
        cmdgen = cmdline_generator(self.param_iter,
                                   PathToBin=self.abs_path_to_bin,
                                   PathToCmd=self.abs_path_to_cmd,
                                   PathsToInputs=self.abs_path_to_input,
                                   PathToOutput=self.abs_path_to_output,
                                   PathToStdout=self.abs_path_to_stdout,
                                   PathToStderr=self.abs_path_to_stderr,
                                   UniqueOutputs=False,
                                   InputParam='-input',
                                   OutputParam='-output')
        bin = self.abs_path_to_bin
        cmd = self.abs_path_to_cmd
        inputfile = self.abs_path_to_input
        outputfile = self.abs_path_to_output
        stdout = self.abs_path_to_stdout
        stderr = self.abs_path_to_stderr

        exp = [' '.join([bin, cmd, '-default=42', '-delimaaachoice1', '-flag1',
                        '-input="%s"' % inputfile, '-output="%s"' % outputfile,
                         '> "%s"' % stdout, '2> "%s"' % stderr])]
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                             '-input="%s"' % inputfile, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                            '-flag1', '-input="%s"' % inputfile,
                             '-output="%s"' % outputfile, '> "%s"' % stdout,
                             '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                             '-input="%s"' % inputfile, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))

        cmdlines = list(cmdgen)
        self.assertEqual(cmdlines, exp)

    def test_cmdline_generator_hard(self):
        """Returns parameter combinations commandlines. Test stdin/stdout"""
        cmdgen = cmdline_generator(self.param_iter,
                                   PathToBin=self.abs_path_to_bin,
                                   PathToCmd=self.abs_path_to_cmd,
                                   PathsToInputs=self.abs_path_to_input,
                                   PathToOutput=self.abs_path_to_output,
                                   PathToStdout=self.abs_path_to_stdout,
                                   PathToStderr=self.abs_path_to_stderr,
                                   UniqueOutputs=True,
                                   InputParam=None,
                                   OutputParam=None)
        bin = self.abs_path_to_bin
        cmd = self.abs_path_to_cmd
        inputfile = self.abs_path_to_input
        outputfile = self.abs_path_to_output
        stdout = self.abs_path_to_stdout
        stderr = self.abs_path_to_stderr

        # the extra '' is intentionally added. When stdout is used for actual
        # output, the stdout_ param gets set to '' which results in an extra
        # space being generated on the cmdline. this should be benign
        # across operating systems
        exp = [' '.join([bin, cmd, '-default=42', '-delimaaachoice1', '-flag1',
                        '< "%s"' % inputfile, '> "%s"0' % outputfile, '',
                         '2> "%s"' % stderr])]
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                            '< "%s"' % inputfile, '> "%s"1' % outputfile, '',
                             '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                            '-flag1', '< "%s"' % inputfile,
                             '> "%s"2' % outputfile, '', '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                            '< "%s"' % inputfile, '> "%s"3' % outputfile, '',
                             '2> "%s"' % stderr]))

        cmdlines = list(cmdgen)
        self.assertEqual(cmdlines, exp)

    def test_cmdline_generator_stdout_stderr_off(self):
        """Returns cmdlines with stdout and stderr disabled"""
        cmdgen = cmdline_generator(self.param_iter,
                                   PathToBin=self.abs_path_to_bin,
                                   PathToCmd=self.abs_path_to_cmd,
                                   PathsToInputs=self.abs_path_to_input,
                                   PathToOutput=self.abs_path_to_output,
                                   PathToStdout=None,
                                   PathToStderr=None,
                                   UniqueOutputs=False,
                                   InputParam='-input',
                                   OutputParam='-output')
        bin = self.abs_path_to_bin
        cmd = self.abs_path_to_cmd
        inputfile = self.abs_path_to_input
        outputfile = self.abs_path_to_output
        stdout = self.abs_path_to_stdout
        stderr = self.abs_path_to_stderr

        exp = [' '.join([bin, cmd, '-default=42', '-delimaaachoice1', '-flag1',
                        '-input="%s"' % inputfile, '-output="%s"' % outputfile,
                         '', ''])]
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                             '-input="%s"' % inputfile, '-output="%s"' % outputfile,
                             '', '']))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                            '-flag1', '-input="%s"' % inputfile,
                             '-output="%s"' % outputfile, '', '']))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                             '-input="%s"' % inputfile, '-output="%s"' % outputfile,
                             '', '']))

        cmdlines = list(cmdgen)
        self.assertEqual(cmdlines, exp)

    def test_cmdline_generator_multiple_inputs(self):
        """Tests the cmdline_generator for multiple input support"""
        paths_to_inputs = ['/some/dir/a', '/some/dir/b']
        cmdgen = cmdline_generator(self.param_iter,
                                   PathToBin=self.abs_path_to_bin,
                                   PathToCmd=self.abs_path_to_cmd,
                                   PathsToInputs=paths_to_inputs,
                                   PathToOutput=self.abs_path_to_output,
                                   PathToStdout=self.abs_path_to_stdout,
                                   PathToStderr=self.abs_path_to_stderr,
                                   UniqueOutputs=False,
                                   InputParam='-input',
                                   OutputParam='-output')
        bin = self.abs_path_to_bin
        cmd = self.abs_path_to_cmd
        inputfile1 = paths_to_inputs[0]
        inputfile2 = paths_to_inputs[1]
        outputfile = self.abs_path_to_output
        stdout = self.abs_path_to_stdout
        stderr = self.abs_path_to_stderr

        exp = [' '.join([bin, cmd, '-default=42', '-delimaaachoice1', '-flag1',
                        '-input="%s"' % inputfile1, '-output="%s"' % outputfile,
                         '> "%s"' % stdout, '2> "%s"' % stderr])]
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                   '-flag1', '-input="%s"' % inputfile2,
                             '-output="%s"' % outputfile, '> "%s"' % stdout,
                             '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                             '-input="%s"' % inputfile1, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                             '-input="%s"' % inputfile2, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                            '-flag1', '-input="%s"' % inputfile1,
                             '-output="%s"' % outputfile, '> "%s"' % stdout,
                             '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                            '-flag1', '-input="%s"' % inputfile2,
                             '-output="%s"' % outputfile, '> "%s"' % stdout,
                             '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                             '-input="%s"' % inputfile1, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                             '-input="%s"' % inputfile2, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))

        cmdlines = list(cmdgen)
        self.assertEqual(cmdlines, exp)

    def test_cmdline_generator_multiple_input_stdin(self):
        """Tests cmdline_generator for multiple inputs over stdin"""
        paths_to_inputs = ['/some/dir/a', '/some/dir/b']
        cmdgen = cmdline_generator(self.param_iter,
                                   PathToBin=self.abs_path_to_bin,
                                   PathToCmd=self.abs_path_to_cmd,
                                   PathsToInputs=paths_to_inputs,
                                   PathToOutput=self.abs_path_to_output,
                                   PathToStdout=self.abs_path_to_stdout,
                                   PathToStderr=self.abs_path_to_stderr,
                                   UniqueOutputs=False,
                                   InputParam=None,
                                   OutputParam='-output')
        bin = self.abs_path_to_bin
        cmd = self.abs_path_to_cmd
        inputfile1 = paths_to_inputs[0]
        inputfile2 = paths_to_inputs[1]
        outputfile = self.abs_path_to_output
        stdout = self.abs_path_to_stdout
        stderr = self.abs_path_to_stderr

        exp = [' '.join([bin, cmd, '-default=42', '-delimaaachoice1', '-flag1',
                        '< "%s"' % inputfile1, '-output="%s"' % outputfile,
                         '> "%s"' % stdout, '2> "%s"' % stderr])]
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                   '-flag1', '< "%s"' % inputfile2,
                             '-output="%s"' % outputfile, '> "%s"' % stdout,
                             '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                             '< "%s"' % inputfile1, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice1',
                             '< "%s"' % inputfile2, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                            '-flag1', '< "%s"' % inputfile1,
                             '-output="%s"' % outputfile, '> "%s"' % stdout,
                             '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                            '-flag1', '< "%s"' % inputfile2,
                             '-output="%s"' % outputfile, '> "%s"' % stdout,
                             '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                             '< "%s"' % inputfile1, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))
        exp.append(' '.join([bin, cmd, '-default=42', '-delimaaachoice2',
                             '< "%s"' % inputfile2, '-output="%s"' % outputfile,
                             '> "%s"' % stdout, '2> "%s"' % stderr]))

        cmdlines = list(cmdgen)
        self.assertEqual(cmdlines, exp)


class CommandLineApplicationTests(TestCase):

    """Tests for the CommandLineApplication class"""

    def setUp(self):
        """setUp for all CommandLineApplication tests"""

        f = open('/tmp/CLAppTester.py', 'w')
        f.write(script)
        f.close()
        system('chmod 777 /tmp/CLAppTester.py')

        # create a copy of the script with a space in the name
        f = open('/tmp/CLApp Tester.py', 'w')
        f.write(script)
        f.close()
        system('chmod 777 "/tmp/CLApp Tester.py"')

        self.app_no_params = CLAppTester()
        self.app_no_params_no_stderr = CLAppTester(SuppressStderr=True)
        self.app_params = CLAppTester({'-F': 'p_file.txt'})
        self.app_params_space_in_command =\
            CLAppTester_space_in_command({'-F': 'p_file.txt'})
        self.app_params_no_stderr = CLAppTester({'-F': 'p_file.txt'},
                                                SuppressStderr=True)
        self.app_params_no_stdout = CLAppTester({'-F': 'p_file.txt'},
                                                SuppressStdout=True)
        self.app_params_input_as_file = CLAppTester({'-F': 'p_file.txt'},
                                                    InputHandler='_input_as_lines')
        self.app_params_WorkingDir = CLAppTester({'-F': 'p_file.txt'},
                                                 WorkingDir='/tmp/test')
        self.app_params_WorkingDir_w_space = CLAppTester({'-F': 'p_file.txt'},
                                                         WorkingDir='/tmp/test space')
        self.app_params_TmpDir = CLAppTester({'-F': 'p_file.txt'},
                                             TmpDir='/tmp/tmp2')
        self.app_params_TmpDir_w_space = CLAppTester({'-F': 'p_file.txt'},
                                                     TmpDir='/tmp/tmp space')
        self.data = 42

    def test_base_command(self):
        """CLAppTester: BaseCommand correctly composed """
        # No parameters on
        app = CLAppTester()
        self.assertEqual(app.BaseCommand, 'cd "/tmp/"; /tmp/CLAppTester.py')
        # ValuedParameter on/off
        app.Parameters['-F'].on('junk.txt')
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; /tmp/CLAppTester.py -F "junk.txt"')
        app.Parameters['-F'].off()
        self.assertEqual(app.BaseCommand, 'cd "/tmp/"; /tmp/CLAppTester.py')
        # ValuedParameter accessed by synonym turned on/off
        app.Parameters['File'].on('junk.txt')
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; /tmp/CLAppTester.py -F "junk.txt"')
        app.Parameters['File'].off()
        self.assertEqual(app.BaseCommand, 'cd "/tmp/"; /tmp/CLAppTester.py')
        # Try multiple parameters, must check for a few different options
        # because parameters are printed in arbitrary order
        app.Parameters['-F'].on('junk.txt')
        app.Parameters['--duh'].on()
        self.failUnless(app.BaseCommand ==
                        'cd "/tmp/"; /tmp/CLAppTester.py -F "junk.txt" --duh'
                        or app.BaseCommand ==
                        'cd "/tmp/"; /tmp/CLAppTester.py --duh -F "junk.txt"')
        # Space in _command
        app = CLAppTester_space_in_command()
        self.assertEqual(app.BaseCommand, 'cd "/tmp/"; "/tmp/CLApp Tester.py"')

    def test_getHelp(self):
        """CLAppTester: getHelp() functions as expected """
        app = CLAppTester()
        self.assertEqual(app.getHelp(), 'Duh')

    def test_handle_app_result_build_failure(self):
        """_handle_app_result_build_failure called when CommandLineAppResult() fails
        """
        app = CLAppTester_bad_fixed_file()
        self.assertRaises(ApplicationError, app)

        app = CLAppTester_bad_fixed_file_w_handler()
        self.assertEqual(app(), "Called self._handle_app_result_build_failure")

    def test_error_on_missing_executable(self):
        """CLAppTester: Useful error message on executable not found
        """

        # fake command via self._command
        class Blah(CLAppTester):
            _command = 'fake_command_jasdlkfsadlkfskladfkladf'

        self.assertRaises(ApplicationNotFoundError, Blah)

        # real command but bad path via self._command
        class Blah(CLAppTester):
            _command = '/not/a/real/path/ls'

        self.assertRaises(ApplicationNotFoundError, Blah)

        # alt _error_on_missing_command function works as expected
        class Blah(CLAppTester):
            _command = 'ls'

            def _error_on_missing_application(self, data):
                raise ApplicationNotFoundError

        self.assertRaises(ApplicationNotFoundError, Blah)

        class Blah(CLAppTester):
            _command = 'fake_app_asfasdasdasdasdasd'

            def _error_on_missing_application(self, data):
                pass

        # no error raised
        Blah()

    def test_no_p_no_d(self):
        """CLAppTester: parameters turned off, no data"""
        app = self.app_no_params
        # test_init
        assert app.Parameters['-F'].isOff()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # test_command
        self.assertEqual(app.BaseCommand, 'cd "/tmp/"; /tmp/CLAppTester.py')
        # test_result
        result = app()
        self.assertEqual(result['StdOut'].read(), 'out\n')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')
        self.assertEqual(result['parameterized_file'], None)
        result.cleanUp()

    def test_no_p_data_as_str(self):
        """CLAppTester: parameters turned off, data as string"""
        app = self.app_no_params
        # test_init
        assert app.Parameters['-F'].isOff()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # test_command
        self.assertEqual(app.BaseCommand, 'cd "/tmp/"; /tmp/CLAppTester.py')
        # test_result
        result = app(self.data)
        self.assertEqual(result['StdOut'].read(), 'out 43\n')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')
        self.assertEqual(result['parameterized_file'], None)
        result.cleanUp()

    def test_p_data_as_str_suppress_stderr(self):
        """CLAppTester: parameters turned on, data as string, suppress stderr"""
        app = self.app_params_no_stderr
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert app.SuppressStderr
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; /tmp/CLAppTester.py -F "p_file.txt"')
        # test_result
        result = app(self.data)
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'], None)
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')
        self.assertEqual(result['parameterized_file'].read(),
                         'out 43 p_file.txt')
        result.cleanUp()

    def test_p_data_as_str_suppress_stdout(self):
        """CLAppTester: parameters turned on, data as string, suppress stdout"""
        app = self.app_params_no_stdout
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert app.SuppressStdout
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; /tmp/CLAppTester.py -F "p_file.txt"')
        # test_result
        result = app(self.data)
        self.assertEqual(result['StdOut'], None)
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')
        self.assertEqual(result['parameterized_file'].read(),
                         'out 43 p_file.txt')
        result.cleanUp()

    def test_p_no_data(self):
        """CLAppTester: parameters turned on, no data"""
        app = self.app_params
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; /tmp/CLAppTester.py -F "p_file.txt"')
        # test_result
        result = app()
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')
        self.assertEqual(result['parameterized_file'].read(),
                         'out p_file.txt')
        result.cleanUp()

    def test_p_space_in_command(self):
        """CLAppTester: parameters turned on, no data, space in command"""
        app = self.app_params_space_in_command
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; "/tmp/CLApp Tester.py" -F "p_file.txt"')
        # test_result
        result = app()
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')
        self.assertEqual(result['parameterized_file'].read(),
                         'out p_file.txt')
        result.cleanUp()

    def test_p_data_as_str(self):
        """CLAppTester: parameters turned on, data as str"""
        app = self.app_params
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; /tmp/CLAppTester.py -F "p_file.txt"')
        # test_result
        result = app(self.data)
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')
        self.assertEqual(result['parameterized_file'].read(),
                         'out 43 p_file.txt')
        result.cleanUp()

    def test_p_data_as_file(self):
        """CLAppTester: parameters turned on, data as file"""
        app = self.app_params_input_as_file
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_lines')
        assert not app.SuppressStderr
        # test_command
        # we don't test the command in this case, because we don't know what
        # the name of the input file is.
        # test_result
        result = app([self.data])
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')
        self.assertEqual(result['parameterized_file'].read(),
                         'out 43 p_file.txt')
        result.cleanUp()

    def test_WorkingDir(self):
        """CLAppTester: WorkingDir functions as expected """
        system('cp /tmp/CLAppTester.py /tmp/test/CLAppTester.py')
        app = self.app_params_WorkingDir
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # WorkingDir is what we expect
        self.assertEqual(app.WorkingDir, '/tmp/test/')
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/test/"; /tmp/CLAppTester.py -F "p_file.txt"')
        # test_result
        result = app()
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')

        # Make sure that the parameterized file is in the correct place
        self.assertEqual(result['parameterized_file'].name,
                         '/tmp/test/p_file.txt')
        self.assertEqual(result['parameterized_file'].read(),
                         'out p_file.txt')
        result.cleanUp()

    def test_WorkingDir_w_space(self):
        """CLAppTester: WorkingDir w/ space in path functions as expected """
        system('cp /tmp/CLAppTester.py "/tmp/test space/CLAppTester.py"')
        app = self.app_params_WorkingDir_w_space
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # WorkingDir is what we expect
        self.assertEqual(app.WorkingDir, '/tmp/test space/')
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/test space/"; /tmp/CLAppTester.py -F "p_file.txt"')
        # test_result
        result = app()
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')

        # Make sure that the parameterized file is in the correct place
        self.assertEqual(result['parameterized_file'].name,
                         '/tmp/test space/p_file.txt')
        self.assertEqual(result['parameterized_file'].read(),
                         'out p_file.txt')
        result.cleanUp()

    def test_TmpDir(self):
        """CLAppTester: Alternative TmpDir functions as expected"""
        app = self.app_params_TmpDir
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # TmpDir is what we expect
        self.assertEqual(app.TmpDir, '/tmp/tmp2')
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; /tmp/CLAppTester.py -F "p_file.txt"')
        # test_result
        result = app()
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')

        # Make sure that the parameterized file is in the correct place
        self.assertEqual(result['parameterized_file'].name,
                         '/tmp/p_file.txt')
        self.assertEqual(result['parameterized_file'].read(),
                         'out p_file.txt')
        result.cleanUp()

    def test_TmpDir_w_space(self):
        """CLAppTester: TmpDir functions as expected w space in name"""
        app = self.app_params_TmpDir_w_space
        # test_init
        assert app.Parameters['-F'].isOn()
        self.assertEqual(app.InputHandler, '_input_as_string')
        assert not app.SuppressStderr
        # TmpDir is what we expect
        self.assertEqual(app.TmpDir, '/tmp/tmp space')
        # test_command
        self.assertEqual(app.BaseCommand,
                         'cd "/tmp/"; /tmp/CLAppTester.py -F "p_file.txt"')
        # test_result
        result = app()
        self.assertEqual(result['StdOut'].read(), '')
        self.assertEqual(result['StdErr'].read(), 'I am stderr\n')
        self.assertEqual(result['ExitStatus'], 0)
        self.assertEqual(result['fixed_file'].read(), 'I am fixed file')
        self.assertEqual(result['base_dep_1'].read(), 'base dependent 1')
        self.assertEqual(result['base_dep_2'].read(), 'base dependent 2')

        # Make sure that the parameterized file is in the correct place
        self.assertEqual(result['parameterized_file'].name,
                         '/tmp/p_file.txt')
        self.assertEqual(result['parameterized_file'].read(),
                         'out p_file.txt')
        result.cleanUp()

    def test_input_as_string(self):
        """CLAppTester: _input_as_string functions as expected """
        self.assertEqual(self.app_no_params._input_as_string('abcd'), 'abcd')
        self.assertEqual(self.app_no_params._input_as_string(42), '42')
        self.assertEqual(self.app_no_params._input_as_string(None), 'None')
        self.assertEqual(self.app_no_params._input_as_string([1]), '[1]')
        self.assertEqual(self.app_no_params._input_as_string({'a': 1}),
                         "{'a': 1}")

    def test_input_as_lines_from_string(self):
        """CLAppTester: _input_as_lines functions as expected w/ data as str
        """
        filename = self.app_no_params._input_as_lines('abcd')
        self.assertEqual(filename[0], '/')
        f = open(filename)
        self.assertEqual(f.readline(), 'a\n')
        self.assertEqual(f.readline(), 'b\n')
        self.assertEqual(f.readline(), 'c\n')
        self.assertEqual(f.readline(), 'd')
        f.close()
        remove(filename)

    def test_input_as_lines_from_list(self):
        """CLAppTester: _input_as_lines functions as expected w/ data as list
        """
        filename = self.app_no_params._input_as_lines(['line 1', None, 3])
        self.assertEqual(filename[0], '/')
        f = open(filename)
        self.assertEqual(f.readline(), 'line 1\n')
        self.assertEqual(f.readline(), 'None\n')
        self.assertEqual(f.readline(), '3')
        f.close()
        remove(filename)

    def test_input_as_lines_from_list_w_newlines(self):
        """CLAppTester: _input_as_lines functions w/ data as list w/ newlines
        """
        filename = self.app_no_params._input_as_lines(['line 1\n', None, 3])
        self.assertEqual(filename[0], '/')
        f = open(filename)
        self.assertEqual(f.readline(), 'line 1\n')
        self.assertEqual(f.readline(), 'None\n')
        self.assertEqual(f.readline(), '3')
        f.close()
        remove(filename)

    def test_input_as_multiline_string(self):
        """CLAppTester: _input_as_multiline_string functions as expected
        """
        filename = self.app_no_params._input_as_multiline_string(
            'line 1\nNone\n3')
        self.assertEqual(filename[0], '/')
        f = open(filename)
        self.assertEqual(f.readline(), 'line 1\n')
        self.assertEqual(f.readline(), 'None\n')
        self.assertEqual(f.readline(), '3')
        f.close()
        remove(filename)

    def test_input_as_lines_from_list_single_entry(self):
        """CLAppTester: _input_as_lines functions as expected w/ 1 element list
        """
        filename = self.app_no_params._input_as_lines(['line 1'])
        self.assertEqual(filename[0], '/')
        f = open(filename)
        self.assertEqual(f.readline(), 'line 1')
        f.close()
        remove(filename)

    def test_input_as_multiline_string_single_line(self):
        """CLAppTester: _input_as_multiline_string functions w/ single line
        """
        # functions as expected with single line string
        filename = self.app_no_params._input_as_multiline_string(
            'line 1')
        self.assertEqual(filename[0], '/')
        f = open(filename)
        self.assertEqual(f.readline(), 'line 1')
        f.close()
        remove(filename)

    def test_getTmpFilename_non_default(self):
        """TmpFilename handles alt tmp_dir, prefix and suffix properly"""
        app = CLAppTester()
        obs = app.getTmpFilename(include_class_id=False)
        self.assertTrue(obs.startswith('/tmp/tmp'))
        self.assertTrue(obs.endswith('.txt'))

        obs = app.getTmpFilename(tmp_dir="/tmp/blah", prefix="app_ctl_test",
                                 suffix='.test', include_class_id=False)
        self.assertTrue(obs.startswith('/tmp/blah/app_ctl_test'))
        self.assertTrue(obs.endswith('.test'))

    def test_getTmpFilename_defaults_to_no_class_id(self):
        """CLAppTester: getTmpFilename doesn't include class id by default
        """
        # I want to explicitly test for this so people don't forget to
        # set the default to False if they change it for testing purposes
        app = CLAppTester()
        self.assertFalse(app.getTmpFilename().
                         startswith('/tmp/tmpCLAppTester'))
        self.assertTrue(app.getTmpFilename(include_class_id=True).
                        startswith('/tmp/tmpCLAppTester'))

    def test_input_as_path(self):
        """CLAppTester: _input_as_path casts data to FilePath"""
        actual = self.app_no_params._input_as_path('test.pdb')
        self.assertEqual(actual, 'test.pdb')
        self.assertEqual(str(actual), '"test.pdb"')
        actual = self.app_no_params._input_as_path('te st.pdb')
        self.assertEqual(actual, 'te st.pdb')
        self.assertEqual(str(actual), '"te st.pdb"')
        actual = self.app_no_params._input_as_path('./test.pdb')
        self.assertEqual(actual, './test.pdb')
        self.assertEqual(str(actual), '"./test.pdb"')
        actual = self.app_no_params._input_as_path('/this/is/a/test.pdb')
        self.assertEqual(actual, '/this/is/a/test.pdb')
        self.assertEqual(str(actual), '"/this/is/a/test.pdb"')
        actual = self.app_no_params._input_as_path('/this/i s/a/test.pdb')
        self.assertEqual(actual, '/this/i s/a/test.pdb')
        self.assertEqual(str(actual), '"/this/i s/a/test.pdb"')

    def test_input_as_paths(self):
        """CLAppTester: _input_as_paths casts each input to FilePath """
        input = ['test.pdb']
        actual = self.app_no_params._input_as_paths(input)
        expected = '"test.pdb"'
        self.assertEqual(actual, expected)

        input = ['test1.pdb', 'test2.pdb']
        actual = self.app_no_params._input_as_paths(input)
        expected = '"test1.pdb" "test2.pdb"'
        self.assertEqual(actual, expected)

        input = ['/path/to/test1.pdb', 'test2.pdb']
        actual = self.app_no_params._input_as_paths(input)
        expected = '"/path/to/test1.pdb" "test2.pdb"'
        self.assertEqual(actual, expected)

        input = ['test1.pdb', '/path/to/test2.pdb']
        actual = self.app_no_params._input_as_paths(input)
        expected = '"test1.pdb" "/path/to/test2.pdb"'
        self.assertEqual(actual, expected)

        input = ['/path/to/test1.pdb', '/path/to/test2.pdb']
        actual = self.app_no_params._input_as_paths(input)
        expected = '"/path/to/test1.pdb" "/path/to/test2.pdb"'
        self.assertEqual(actual, expected)

        input = ['/pa th/to/test1.pdb', '/path/to/te st2.pdb']
        actual = self.app_no_params._input_as_paths(input)
        expected = '"/pa th/to/test1.pdb" "/path/to/te st2.pdb"'
        self.assertEqual(actual, expected)

    def test_absolute(self):
        """CLAppTester: _absolute converts relative paths to absolute paths
        """
        absolute = self.app_no_params._absolute
        self.assertEqual(absolute('/tmp/test.pdb'), '/tmp/test.pdb')
        self.assertEqual(absolute('test.pdb'), '/tmp/test.pdb')

    def test_working_dir_setting(self):
        """CLAppTester: WorkingDir is set correctly """
        app = CLAppTester_no_working_dir()
        self.assertEqual(app.WorkingDir, getcwd() + '/')

    def test_error_raised_on_command_None(self):
        """CLAppTester: An error is raises when _command == None """
        app = CLAppTester()
        app._command = None
        self.assertRaises(ApplicationError, app._get_base_command)

    def test_rejected_exit_status(self):
        """CLAppTester_reject_exit_status results in useful error """
        app = CLAppTester_reject_exit_status()
        self.assertRaises(ApplicationError, app)

    def test_getTmpFilename(self):
        """TmpFilename should return filename of correct length"""
        app = CLAppTester()
        obs = app.getTmpFilename(include_class_id=True)
        # leaving the strings in this statement so it's clear where the expected
        # length comes from
        self.assertEqual(len(obs), len(app.TmpDir) + len('/') + app.TmpNameLen
                         + len('tmp') + len('CLAppTester') + len('.txt'))
        assert obs.startswith(app.TmpDir)
        chars = set(obs[18:])
        assert len(chars) > 1

        obs = app.getTmpFilename(include_class_id=False)
        # leaving the strings in this statement so it's clear where the expected
        # length comes from
        self.assertEqual(len(obs), len(app.TmpDir) + len('/') + app.TmpNameLen
                         + len('tmp') + len('.txt'))
        assert obs.startswith(app.TmpDir)

    def test_getTmpFilename_prefix_suffix_result_constructor(self):
        """TmpFilename: result has correct prefix, suffix, type"""
        app = CLAppTester()
        obs = app.getTmpFilename(prefix='blah', include_class_id=False)
        self.assertTrue(obs.startswith('/tmp/blah'))
        obs = app.getTmpFilename(suffix='.blah', include_class_id=False)
        self.assertTrue(obs.endswith('.blah'))
        # Prefix defaults to not include the class name
        obs = app.getTmpFilename(include_class_id=False)
        self.assertFalse(obs.startswith('/tmp/tmpCLAppTester'))
        self.assertTrue(obs.endswith('.txt'))
        # including class id functions correctly
        obs = app.getTmpFilename(include_class_id=True)
        self.assertTrue(obs.startswith('/tmp/tmpCLAppTester'))
        self.assertTrue(obs.endswith('.txt'))

        # result as FilePath
        obs = app.getTmpFilename(result_constructor=FilePath)
        self.assertEqual(type(obs), FilePath)

        # result as str (must check that result is a str and is not a FilePath
        # since a FilePath is a str)
        obs = app.getTmpFilename(result_constructor=str)
        self.assertEqual(type(obs), str)
        self.assertNotEqual(type(obs), FilePath)


class ConvenienceFunctionTests(TestCase):

    """
    """

    def setUp(self):
        """
        """
        self.tmp_dir = gettempdir()
        self.tmp_name_len = 20

    def test_guess_input_handler(self):
        """guess_input_handler should correctly identify input"""
        gih = guess_input_handler
        self.assertEqual(gih('abc.txt'), '_input_as_string')
        self.assertEqual(gih('>ab\nTCAG'), '_input_as_multiline_string')
        self.assertEqual(gih(['ACC', 'TGA'], True), '_input_as_seqs')
        self.assertEqual(gih(['>a', 'ACC', '>b', 'TGA']), '_input_as_lines')
        self.assertEqual(gih([('a', 'ACC'), ('b', 'TGA')]),
                         '_input_as_seq_id_seq_pairs')
        self.assertEqual(gih([]), '_input_as_lines')

    def test_get_tmp_filename(self):
        """get_tmp_filename should return filename of correct length

            Adapted from the CommandLineApplication tests of the member function

        """
        obs = get_tmp_filename()
        # leaving the strings in this statement so it's clear where the expected
        # length comes from
        self.assertEqual(len(obs), len(self.tmp_dir) + len('/') + self.tmp_name_len
                         + len('tmp') + len('.txt'))
        self.assertTrue(obs.startswith(self.tmp_dir))

        # different results on different calls
        self.assertNotEqual(get_tmp_filename(), get_tmp_filename())

        obs = get_tmp_filename()
        # leaving the strings in this statement so it's clear where the expected
        # length comes from
        self.assertEqual(len(obs), len(self.tmp_dir) + len('/') + self.tmp_name_len
                         + len('tmp') + len('.txt'))
        assert obs.startswith(self.tmp_dir)

    def test_get_tmp_filename_prefix_suffix_constructor(self):
        """get_tmp_filename: result has correct prefix, suffix, type

            Adapted from the CommandLineApplication tests of the member function
        """
        obs = get_tmp_filename(prefix='blah')
        self.assertTrue(obs.startswith('%s/blah' % self.tmp_dir))
        obs = get_tmp_filename(suffix='.blah')
        self.assertTrue(obs.endswith('.blah'))

        # result as FilePath
        obs = get_tmp_filename(result_constructor=FilePath)
        self.assertEqual(type(obs), FilePath)

        # result as str (must check that result is a str and is not a FilePath
        # since a FilePath is a str)
        obs = get_tmp_filename(result_constructor=str)
        self.assertEqual(type(obs), str)
        self.assertNotEqual(type(obs), FilePath)

    def test_app_path(self):
        """app_path should return correct paths"""
        self.assertEqual(app_path('ls'), '/bin/ls')
        self.assertEqual(app_path('lsxxyyx'), False)


class RemoveTests(TestCase):

    def test_remove(self):
        """This will remove the test script. Not actually a test!"""

        for dir, n, fnames in walk('/tmp/test/'):
            for f in fnames:
                try:
                    remove(dir + f)
                except OSError as e:
                    pass

        remove('/tmp/CLAppTester.py')
        remove('/tmp/test space/CLAppTester.py')
        remove('/tmp/CLApp Tester.py')
        rmdir('/tmp/tmp space')
        rmdir('/tmp/test')
        rmdir('/tmp/test space')
        rmdir('/tmp/tmp2')
        rmdir('/tmp/blah')


#=====================END OF TESTS===================================

script = """#!/usr/bin/env python
#This is a test script intended to test the CommandLineApplication
#class and CommandLineAppResult class

from sys import argv, stderr,stdin
from os import isatty

out_file_name = None

input_arg = None

# parse input
try:
    if argv[1] == '-F':
        out_file_name = argv[2]
except IndexError:
    pass
try:
    if out_file_name:
        input_arg = argv[3]
    else:
        input_arg = argv[1]
except IndexError:
    pass
# Create the output string
out = 'out'
# get input
try:
    f = open(str(input_arg))
    data = int(f.readline().strip())
except IOError:
    try:
        data = int(input_arg)
    except TypeError:
        data = None

if data:
    data = str(data + 1)
    out = ' '.join([out,data])

# Write base dependent output files
base = 'BASE'
f = open('/tmp/' + base + '.1','w')
f.writelines(['base dependent 1'])
f.close()
f = open('/tmp/' + base + '.2','w')
f.writelines(['base dependent 2'])
f.close()

# If output to file, open the file and write output to it
if out_file_name:
    filename = argv[2]
    f = open(''.join([out_file_name]),'w')
    out = ' '.join([out,out_file_name])
    f.writelines(out)
    f.close()
else:
    print out

#generate some stderr
print >> stderr, 'I am stderr'

# Write the fixed file
f = open('/tmp/fixed.txt','w')
f.writelines(['I am fixed file'])
f.close()

"""


class CLAppTester(CommandLineApplication):
    _parameters = {
        '-F': ValuedParameter(Prefix='-', Name='F', Delimiter=' ',
                              Value=None, Quote="\""),
        '--duh': FlagParameter(Prefix='--', Name='duh')}
    _command = '/tmp/CLAppTester.py'
    _synonyms = {'File': '-F', 'file': '-F'}
    _working_dir = '/tmp'

    def _get_result_paths(self, data):

        if self.Parameters['-F'].isOn():
            param_path = ''.join(
                [self.WorkingDir, self.Parameters['-F'].Value])
        else:
            param_path = None

        result = {}
        result['fixed_file'] = ResultPath(Path='/tmp/fixed.txt')
        result['parameterized_file'] = ResultPath(Path=param_path,
                                                  IsWritten=self.Parameters['-F'].isOn())
        result['base_dep_1'] = ResultPath(Path=self._build_name(suffix='.1'))
        result['base_dep_2'] = ResultPath(Path=self._build_name(suffix='.2'))
        return result

    def _build_name(self, suffix):
        return '/tmp/BASE' + suffix

    def getHelp(self):
        return """Duh"""


class CLAppTester_no_working_dir(CLAppTester):
    _working_dir = None


class CLAppTester_reject_exit_status(CLAppTester):

    def _accept_exit_status(self, exit_status):
        return False


class CLAppTester_bad_fixed_file(CLAppTester):

    def _get_result_paths(self, data):

        if self.Parameters['-F'].isOn():
            param_path = ''.join(
                [self.WorkingDir, self.Parameters['-F'].Value])
        else:
            param_path = None

        result = {}
        result['fixed_file'] = ResultPath(Path='/tmp/fixed.txt')
        result['fixed_file_bad'] = ResultPath(Path='/tmp/i_dont_exist.txt')
        result['parameterized_file'] = ResultPath(Path=param_path,
                                                  IsWritten=self.Parameters['-F'].isOn())
        result['base_dep_1'] = ResultPath(Path=self._build_name(suffix='.1'))
        result['base_dep_2'] = ResultPath(Path=self._build_name(suffix='.2'))
        return result


class CLAppTester_bad_fixed_file_w_handler(CLAppTester_bad_fixed_file):

    def _handle_app_result_build_failure(
            self,
            out,
            err,
            exit_status,
            result_paths):
        return "Called self._handle_app_result_build_failure"


class CLAppTester_space_in_command(CLAppTester):
    _command = '"/tmp/CLApp Tester.py"'


class ParameterCombinationsApp(CommandLineApplication):

    """ParameterCombinations mock application to wrap"""
    _command = 'testcmd'
    _parameters = {'-flag1': FlagParameter(Prefix='-', Name='flag1'),
                   '-flag2': FlagParameter(Prefix='-', Name='flag2'),
                   '--value1': ValuedParameter(Prefix='--', Name='value1'),
                   '-value2': ValuedParameter(Prefix='-', Name='value2'),
                   '-mix1': MixedParameter(Prefix='-', Name='mix1'),
                   '-mix2': MixedParameter(Prefix='-', Name='mix2'),
                   '-delim': ValuedParameter(Prefix='-', Name='delim',
                                             Delimiter='aaa'),
                   '-default': ValuedParameter(Prefix='-', Name='default',
                                               Value=42, Delimiter='='),
                   '-input': ValuedParameter(Prefix='-', Name='input',
                                             Delimiter='='),
                   '-output': ValuedParameter(Prefix='-', Name='output',
                                              Delimiter='=')}
if __name__ == '__main__':

    main()
