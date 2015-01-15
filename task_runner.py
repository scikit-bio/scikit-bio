#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import sys
import argparse
import subprocess
import yaml
from urllib import urlretrieve

TEST_DIRNAME = ".test_sandbox"
CONDA_DIRNAME = "conda"
#CONDA_SERVER = "http://repo.continuum.io/miniconda/"
CONDA_SERVER = "http://localhost:8000/"
CONDA_FILENAME = "Miniconda3-3.7.3-Linux-x86_64.sh"

SKBIO_DIR = os.getcwd()
TEST_DIR = os.path.join(SKBIO_DIR, TEST_DIRNAME)
CONDA_DIR = os.path.join(TEST_DIR, CONDA_DIRNAME)

class Environment(object):
    def __init__(self, env_vars):
        self._vars = env_vars
        self._prev_var_vals = {}

    def set_up(self):
        for var_name in self._vars:
            self._prev_var_vals[var_name] = os.environ.get(var_name)
            os.environ[var_name] = self._vars[var_name]

    def tear_down(self):
        for var_name in self._prev_var_vals:
            if self._prev_var_vals[var_name] is not None:
                os.environ[var_name] = self._prev_var_vals[var_name]
            else:
                del os.environ[var_name]

    def __str__(self):
        return '\n'.join([ str((key,val)) for key,val in self._vars.items() ])


class Environment_Manager(object):
    def __init__(self):
        with open(".travis.yml", "r") as stream:
            travis_config = yaml.load(stream)
        self._envs = []
        for line in travis_config["env"]:
            env_vars = {}
            for var in line.split():
                var_name, var_val = self._read_var(var)
                env_vars[var_name] = var_val
            env = Environment(env_vars)
            self._envs.append(env)

    def _read_var(self, line):
        eq_idx = line.index('=')
        var = line[:eq_idx]
        val = line[eq_idx+1:].replace('\'', '')
        return var, val

    def get_environments(self):
        return self._envs
        
class Runner(object):
    def __init__(self):
        self._conda_path = os.path.join(CONDA_DIR, 'bin', 'conda')

    def _run_shell_command(self, command):
        print(command)
        subprocess.check_call(command, stdout=sys.stdout)

    def _install_conda(self):
        if not os.path.exists(TEST_DIR):
            os.mkdir(TEST_DIR)
        conda_installer_path = os.path.join(TEST_DIR, CONDA_FILENAME)
        if not os.path.exists(self._conda_path):
            if not os.path.exists(conda_installer_path):
                conda_url = CONDA_SERVER + CONDA_FILENAME
                urlretrieve(conda_url, conda_installer_path)
            install_cmd = [
                'bash', conda_installer_path, '-b',
                '-p', CONDA_DIR
            ]
            self._run_shell_command(install_cmd)

    def _setup_environment(self):
        python_version = os.environ.get('PYTHON_VERSION') 
        numpy_version = os.environ.get('NUMPY_VERSION') 
        mpl_version = os.environ.get('MATPLOTLIB_VERSION')
        env_name = 'py_'+python_version
        if len(numpy_version) > 0:
            env_name += '_np_'+numpy_version
        if len(mpl_version) > 0:
            env_name += '_mpl_'+mpl_version
        self._env_name = env_name
            
        #activate_script = os.path.join(CONDA_DIR, 'bin', 'activate')
        #self._run_shell_command(['bash', '-c', 'source', activate_script, env_name])

        self._pip_path = os.path.join(CONDA_DIR, 'envs', env_name, 'bin', 'pip')
        self._nosetests_path = os.path.join(CONDA_DIR, 'envs', env_name, 'bin', 'nosetests')
        env_dir = os.path.join(CONDA_DIR, 'envs', env_name)
        if not os.path.exists(env_dir):
            shell_command = [
                self._conda_path, 'create', '--yes',
                '-n', env_name,
                'python=%s' %  python_version,
                'pip',
                'numpy%s' % numpy_version,
                'scipy',
                'matplotlib%s' % mpl_version,
                'pandas',
                'nose',
                'pep8',
                'Sphinx=1.2.2',
                'IPython',
            ]
            self._run_shell_command(shell_command)

            use_cython = os.environ.get('USE_CYTHON')
            if use_cython:
                shell_command = [
                    self._conda_path, 'install', '--yes',
                    '-n', env_name,
                    'cython',
                ]
                self._run_shell_command(shell_command)


            shell_command = [
                self._pip_path, 'install',
                'sphinx-bootstrap-theme',
                'future',
                'six',
                'coveralls',
                'natsort',
                'pyflakes',
                'flake8',
                'python-dateutil',
            ]
            self._run_shell_command(shell_command)

    def _install(self):
        skbio_env_path = os.path.join(TEST_DIR, 'skbio-%s' % self._env_name)
        if not os.path.exists(skbio_env_path):
            os.mkdir(skbio_env_path)
            copy_command = [
                'rsync', '-av', '--exclude=.*',
                os.path.normpath(SKBIO_DIR) + os.sep, skbio_env_path,
            ]
            self._run_shell_command(copy_command)
            install_command = [
                self._pip_path, 'install', '-e', skbio_env_path, '--no-deps',
            ]
            self._run_shell_command(install_command)

    def _tests(self, target):
        prev_val = os.environ.get('PYTHONWARNINGS')
        os.environ['PYTHONWARNINGS'] = 'ignore'
        test_command = [
            self._nosetests_path,
            target,
            '--with-coverage',
            '-I DONOTIGNOREANYTHING',
        ]
        
        with_doctest = os.environ.get('WITH_DOCTEST')
        if with_doctest == 'True':
            test_command.append('--with-doctest')

        self._run_shell_command(test_command)
        self._run_shell_command(['pep8', target, 'setup.py', 'checklist.py'])
        self._run_shell_command(['flake8', target, 'setup.py', 'checklist.py'])
        self._run_shell_command(['./checklist.py'])

        if prev_val is not None:
            os.environ['PYTHONWARNINGS'] = prev_val
        else:
            del os.environ['PYTHONWARNINGS']

        #self._run_shell_command(['bash', '-c', 'source', 'deactivate'])

    def _docs(self):
        os.chdir(os.path.join(SKBIO_DIR, 'doc'))
        self._run_shell_command(['make', 'clean'])
        self._run_shell_command(['make', 'html'])
        os.chdir(SKBIO_DIR)

    def _coverage(self):
        self._run_shell_command(['coveralls'])

    def run(self, target):
        self._install_conda()
        self._setup_environment()
        self._install()
        self._tests(target)
        self._coverage()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simple task runner.')
    parser.add_argument('target')
    args = parser.parse_args()
    target = args.target

    is_travis = os.environ.get('TRAVIS') == 'TRUE'

    # If running under travis, it will take care of setting environment
    # variables. Otherwise, we need to parse .travis.yml and set the
    # variables manually
    if is_travis:
        runner = Runner()
        runner.run()
    else:
        envs = Environment_Manager().get_environments()
        for env in envs:
            print
            print(env)
            print
            env.set_up()
            runner = Runner()
            runner.run(target)
            env.tear_down()

