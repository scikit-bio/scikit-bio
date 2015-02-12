#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import platform
from setuptools import find_packages, setup
from setuptools.extension import Extension

import numpy as np

__version__ = "0.2.2-dev"

classes = """
    Development Status :: 1 - Planning
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python
    Programming Language :: Python :: 2
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3.3
    Programming Language :: Python :: 3.4
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Data structures, algorithms and educational '
               'resources for bioinformatics.')

with open('README.rst') as f:
    long_description = f.read()

# Dealing with Cython
USE_CYTHON = os.environ.get('USE_CYTHON', False)
ext = '.pyx' if USE_CYTHON else '.c'

# There's a bug in some versions of Python 3.4 that propagates
# -Werror=declaration-after-statement to extensions, instead of just affecting
# the compilation of the interpreter. See http://bugs.python.org/issue21121 for
# details. This acts as a workaround until the next Python 3 release -- thanks
# Wolfgang Maier (wolma) for the workaround!
ssw_extra_compile_args = ['-Wno-error=declaration-after-statement']

# Users with i686 architectures have reported that adding this flag allows
# SSW to be compiled. See https://github.com/biocore/scikit-bio/issues/409 and
# http://stackoverflow.com/q/26211814/3776794 for details.
if platform.machine() == 'i686':
    ssw_extra_compile_args.append('-msse2')

extensions = [
    Extension("skbio.stats.__subsample",
              ["skbio/stats/__subsample" + ext]),
    Extension("skbio.alignment._ssw_wrapper",
              ["skbio/alignment/_ssw_wrapper" + ext,
               "skbio/alignment/_lib/ssw.c"],
              extra_compile_args=ssw_extra_compile_args)
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup(name='scikit-bio',
      version=__version__,
      license='BSD',
      description=description,
      long_description=long_description,
      author="scikit-bio development team",
      author_email="gregcaporaso@gmail.com",
      maintainer="scikit-bio development team",
      maintainer_email="gregcaporaso@gmail.com",
      url='http://scikit-bio.org',
      test_suite='nose.collector',
      packages=find_packages(),
      ext_modules=extensions,
      include_dirs=[np.get_include()],
      install_requires=['numpy >= 1.7', 'matplotlib >= 1.1.0',
                        'scipy >= 0.13.0', 'pandas', 'future', 'six',
                        'natsort', 'IPython'],
      extras_require={'test': ["nose >= 0.10.1", "pep8", "flake8",
                               "python-dateutil"],
                      'doc': ["Sphinx == 1.2.2", "sphinx-bootstrap-theme"]},
      classifiers=classifiers,
      package_data={
          'skbio.io.tests': ['data/*'],
          'skbio.stats.tests': ['data/*'],
          'skbio.stats.distance.tests': ['data/*'],
          'skbio.stats.ordination.tests': ['data/*'],
          'skbio.parse.sequences.tests': ['data/*'],
          }
      )
