#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__version__ = '0.0.0-dev'

from setuptools import find_packages, setup
from distutils.command.build_py import build_py

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

long_description = """The BiPy Project"""

setup(name='bipy',
      cmdclass={'build_py': build_py},
      version=__version__,
      license='BSD',
      description='bipy',
      long_description=long_description,
      author="The BiPy Developers",
      author_email="gregcaporaso@gmail.com",
      maintainer="The BiPy Developers",
      maintainer_email="gregcaporaso@gmail.com",
      url='https://github.com/biocore/bipy',
      packages=find_packages(),
      install_requires=['numpy >= 1.5.1', 'matplotlib >= 1.1.0',
                        'scipy >=0.13.0'],
      extras_require={'test': ["nose >= 0.10.1", "pep8"],
                      'doc': "Sphinx >= 0.3"},
      classifiers=classifiers
      )
