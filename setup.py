#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The BiPy project"
__credits__ = ["Greg Caporaso", "Rob Knight", "Daniel McDonald",
               "Jai Ram Rideout", "Antiono Gonzalez", "Yoshiki Vazquez Baeza",
               "Emily TerAvest"]
__license__ = "BSD"
__version__ = '0.0.0-dev'
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from setuptools import setup
import sys

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

setup(name='BiPy',
      cmdclass={'build_py':build_py},
      version=__version__,
      license=__license__,
      description='BiPy',
      long_description=long_description,
      author=__maintainer__,
      author_email=__email__,
      maintainer=__maintainer__,
      maintainer_email=__email__,
      url='https://github.com/gregcaporaso/bipy', # will soon be replaced
      packages=['bipy'],
      install_requires=['numpy >= 1.5.1, <=1.7.1', 'matplotlib >= 1.1.0',
                        'scipy'],
      extras_require={'test':["nose >= 0.10.1", "tox >= 1.6.1"],
                      'doc':"Sphinx >= 0.3"
                     },
      classifiers=classifiers
      )
