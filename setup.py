#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import platform
import re
import ast
import sys
import sysconfig
import subprocess

from setuptools import find_packages, setup
from setuptools.extension import Extension

import numpy as np


if sys.version_info.major != 3:
    sys.exit("scikit-bio can only be used with Python 3. You are currently "
             "running Python %d." % sys.version_info.major)

clang = False
icc = False
try:
    if os.environ['CC'] == "clang":
        clang = True
except KeyError:
    pass

if not clang:
    try:
        if subprocess.check_output(
                ["gcc", "--version"],
                universal_newlines=True).find("clang") != -1:
            clang = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

try:
    if os.environ['CC'] == "icc":
        icc = True
except KeyError:
    pass

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('skbio/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: 3.9
    Programming Language :: Python :: 3.10
    Programming Language :: Python :: 3.11
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
USE_CYTHON = os.environ.get('USE_CYTHON')
if USE_CYTHON is None or USE_CYTHON.lower() in {'false', 'no'}:
    USE_CYTHON = False
else:
    USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.c'

ssw_extra_compile_args = ['-I.']

if platform.system() != 'Windows':
    if icc or sysconfig.get_config_vars()['CC'] == 'icc':
        ssw_extra_compile_args.extend(['-qopenmp-simd',
                                       '-DSIMDE_ENABLE_OPENMP'])
    elif not (clang or sysconfig.get_config_vars()['CC'] == 'clang'):
        ssw_extra_compile_args.extend(['-fopenmp-simd',
                                       '-DSIMDE_ENABLE_OPENMP'])
elif platform.system() == 'Windows':
    ssw_extra_compile_args.extend(['-openmp:experimental'])

# Users with i686 architectures have reported that adding this flag allows
# SSW to be compiled. See https://github.com/biocore/scikit-bio/issues/409 and
# http://stackoverflow.com/q/26211814/3776794 for details.
if platform.machine() == 'i686':
    ssw_extra_compile_args.append('-msse2')

extensions = [
    Extension("skbio.metadata._intersection",
              ["skbio/metadata/_intersection" + ext]),
    Extension("skbio.stats.__subsample",
              ["skbio/stats/__subsample" + ext],
              include_dirs=[np.get_include()]),
    Extension("skbio.alignment._ssw_wrapper",
              ["skbio/alignment/_ssw_wrapper" + ext,
               "skbio/alignment/_lib/ssw.c"],
              extra_compile_args=ssw_extra_compile_args,
              include_dirs=[np.get_include()]),
    Extension("skbio.diversity._phylogenetic",
              ["skbio/diversity/_phylogenetic" + ext],
              include_dirs=[np.get_include()]),
    Extension("skbio.stats.ordination._cutils",
              ["skbio/stats/ordination/_cutils" + ext],
              extra_compile_args=ssw_extra_compile_args),
    Extension("skbio.stats.distance._cutils",
              ["skbio/stats/distance/_cutils" + ext],
              extra_compile_args=ssw_extra_compile_args),
]

if USE_CYTHON:
    from Cython.Build import cythonize
    # Always recompile the pyx files to C if USE_CYTHON is set.
    extensions = cythonize(extensions, force=True)

setup(name='scikit-bio',
      version=version,
      license='BSD-3-Clause',
      description=description,
      long_description=long_description,
      author="scikit-bio development team",
      author_email="gregcaporaso@gmail.com",
      maintainer="scikit-bio development team",
      maintainer_email="gregcaporaso@gmail.com",
      url='http://scikit-bio.org',
      packages=find_packages(),
      ext_modules=extensions,
      include_dirs=[np.get_include()],
      tests_require=['pytest', 'coverage'],
      install_requires=[
          'lockfile >= 0.10.2',  # req'd for our usage of CacheControl
          'CacheControl >= 0.11.5',
          'decorator >= 3.4.2',
          'IPython >= 3.2.0',
          'matplotlib >= 1.4.3',
          'natsort >= 4.0.3',
          'numpy >= 1.9.2',
          'pandas >= 1.5.0',
          'scipy >= 1.9.0',
          'h5py >= 3.6.0',
          'hdmedians >= 0.14.1',
      ],
      classifiers=classifiers,
      package_data={
          'skbio.diversity.alpha.tests': ['data/qiime-191-tt/*'],
          'skbio.diversity.beta.tests': ['data/qiime-191-tt/*'],
          'skbio.io.tests': ['data/*'],
          'skbio.io.format.tests': ['data/*'],
          'skbio.stats.tests': ['data/*'],
          'skbio.stats.distance.tests': ['data/*'],
          'skbio.stats.ordination.tests': ['data/*']
          }
      )
