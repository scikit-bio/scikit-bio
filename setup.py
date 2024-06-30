#!/usr/bin/env python

"""Setup script for scikit-bio installation.

----------------------------------------------------------------------------
Copyright (c) 2013--, scikit-bio development team.

Distributed under the terms of the Modified BSD License.

The full license is in the file LICENSE.txt, distributed with this software.
----------------------------------------------------------------------------
"""

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
from Cython.Build import cythonize


if sys.version_info.major != 3:
    sys.exit(
        "scikit-bio can only be used with Python 3. You are currently "
        "running Python %d." % sys.version_info.major
    )


def check_bin(ccbin, source, allow_dash):
    """Check if a given compiler matches the specified name."""
    # remove any parameters (e.g. gcc -I /a/b/c -> gcc)
    source0 = source.split()[0]
    # remove any path
    bsource = os.path.basename(source0)
    # now let's go search for ccbin
    if allow_dash:
        found = False
        # allow for the ccbin to be between - (e.g. gcc-1.2)
        for el in bsource.split("-"):
            if el == ccbin:
                found = True
                break
    else:
        found = bsource == ccbin
    return found


# Note: We are looking for Apple/MacOS clang, which does not support omp
#       Will treat "real clang" (e.g. llvm based) same as gcc
clang = False

# icc uses slightly different omp cmdline arguments
icc = False

# Are we using the default gcc as the compiler?
gcc = True
try:
    if os.environ["CC"] == "gcc":
        gcc = True
    elif os.environ["CC"] != "":
        gcc = False
except KeyError:
    pass


if not gcc:
    try:
        if check_bin("clang", os.environ["CC"], False):
            # note, the conda provideed clang is not detected here
            # and this is on purpose, as MacOS clang is very different
            # than conda-provised one (which is llvm based)
            # so do not look for substrings
            # (e.g. do not match x86_64-apple-darwin13.4.0-clang)
            clang = True
        elif check_bin("icc", os.environ["CC"], True):
            icc = True
    except KeyError:
        pass
else:
    try:
        if check_bin("clang", sysconfig.get_config_vars()["CC"], False):
            # as above
            clang = True
            gcc = False
        elif check_bin("icc", sysconfig.get_config_vars()["CC"], True):
            icc = True
            gcc = False
    except KeyError:
        pass

if gcc:
    # check if the default gcc is just a wrapper around clang
    try:
        if (
            subprocess.check_output(["gcc", "--version"], universal_newlines=True).find(
                "clang"
            )
            != -1
        ):
            clang = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r"__version__\s+=\s+(.*)")

with open("skbio/__init__.py", "rb") as f:
    hit = _version_re.search(f.read().decode("utf-8")).group(1)
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
    Programming Language :: Python :: 3.12
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
    Operating System :: Microsoft :: Windows
"""
classifiers = [s.strip() for s in classes.split("\n") if s]

description = (
    "Data structures, algorithms and educational " "resources for bioinformatics."
)

with open("README.rst") as f:
    long_description = f.read()


# Compile SSW module
ssw_extra_compile_args = ["-I."]

if platform.system() != "Windows":
    if icc:
        ssw_extra_compile_args.extend(["-qopenmp-simd", "-DSIMDE_ENABLE_OPENMP"])
    elif not clang:
        ssw_extra_compile_args.extend(["-fopenmp-simd", "-DSIMDE_ENABLE_OPENMP"])
elif platform.system() == "Windows":
    ssw_extra_compile_args.extend(["-openmp:experimental"])

stats_extra_compile_args = [] + ssw_extra_compile_args
stats_extra_link_args = []
if platform.system() != "Windows":
    if icc:
        stats_extra_compile_args.extend(["-qopenmp"])
        stats_extra_link_args.extend(["-qopenmp"])
    elif not clang:
        stats_extra_compile_args.extend(["-fopenmp"])
        stats_extra_link_args.extend(["-fopenmp"])

# Users with i686 architectures have reported that adding this flag allows
# SSW to be compiled. See https://github.com/scikit-bio/scikit-bio/issues/409
# and http://stackoverflow.com/q/26211814/3776794 for details.
if platform.machine() == "i686":
    ssw_extra_compile_args.append("-msse2")


# Cython modules (*.pyx). They will be compiled into C code (*.c) during build.
ext = ".pyx"
extensions = [
    Extension("skbio.metadata._intersection", ["skbio/metadata/_intersection" + ext]),
    Extension(
        "skbio.alignment._ssw_wrapper",
        ["skbio/alignment/_ssw_wrapper" + ext, "skbio/alignment/_lib/ssw.c"],
        extra_compile_args=ssw_extra_compile_args,
        include_dirs=[np.get_include()],
    ),
    Extension(
        "skbio.diversity._phylogenetic",
        ["skbio/diversity/_phylogenetic" + ext],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "skbio.stats.ordination._cutils",
        ["skbio/stats/ordination/_cutils" + ext],
        extra_compile_args=stats_extra_compile_args,
        extra_link_args=stats_extra_link_args,
    ),
    Extension(
        "skbio.stats.distance._cutils",
        ["skbio/stats/distance/_cutils" + ext],
        extra_compile_args=stats_extra_compile_args,
        extra_link_args=stats_extra_link_args,
    ),
]

extensions = cythonize(extensions, force=True)


setup(
    name="scikit-bio",
    version=version,
    license="BSD-3-Clause",
    description=description,
    long_description=long_description,
    author="scikit-bio development team",
    author_email="qiyunzhu@gmail.com",
    maintainer="scikit-bio development team",
    maintainer_email="qiyunzhu@gmail.com",
    url="https://scikit.bio",
    packages=find_packages(),
    ext_modules=extensions,
    include_dirs=[np.get_include()],
    tests_require=["pytest", "coverage"],
    install_requires=[
        "requests >= 2.20.0",
        "decorator >= 3.4.2",
        "natsort >= 4.0.3",
        "numpy >= 1.17.0",
        "pandas >= 1.5.0",
        "scipy >= 1.9.0",
        "h5py >= 3.6.0",
        "biom-format >= 2.1.16",
        "statsmodels >= 0.14.0",
    ],
    classifiers=classifiers,
    package_data={
        "skbio.diversity.alpha.tests": ["data/qiime-191-tt/*"],
        "skbio.diversity.beta.tests": ["data/qiime-191-tt/*"],
        "skbio.io.tests": ["data/*"],
        "skbio.io.format.tests": ["data/*"],
        "skbio.stats.tests": ["data/*"],
        "skbio.stats.distance.tests": ["data/*"],
        "skbio.stats.ordination.tests": ["data/*"],
        "skbio.metadata.tests": ["data/invalid/*", "data/valid/*"],
        "skbio.embedding.tests": ["data/*"],
    },
)
