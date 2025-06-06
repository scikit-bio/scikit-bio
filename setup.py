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

# Compile SSW module
extra_compile_args = ["-I."]

if platform.system() != "Windows":
    if icc:
        extra_compile_args.extend(["-qopenmp-simd", "-DSIMDE_ENABLE_OPENMP"])
    elif not clang:
        extra_compile_args.extend(["-fopenmp-simd", "-DSIMDE_ENABLE_OPENMP"])
elif platform.system() == "Windows":
    extra_compile_args.extend(["-openmp:experimental"])

stats_extra_compile_args = [] + extra_compile_args
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
    extra_compile_args.append("-msse2")


# Cython modules (*.pyx). They will be compiled into C code (*.c) during build.
ext = ".pyx"
extensions = [
    Extension("skbio.metadata._intersection", ["skbio/metadata/_intersection" + ext]),
    Extension(
        "skbio.tree._c_nj",
        ["skbio/tree/_c_nj" + ext],
        extra_compile_args=stats_extra_compile_args,
        extra_link_args=stats_extra_link_args,
    ),
    Extension(
        "skbio.tree._c_me",
        ["skbio/tree/_c_me" + ext],
        extra_compile_args=stats_extra_compile_args,
        extra_link_args=stats_extra_link_args,
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
    Extension(
        "skbio.alignment._cutils",
        ["skbio/alignment/_cutils" + ext],
        extra_compile_args=stats_extra_compile_args,
        extra_link_args=stats_extra_link_args,
        include_dirs=[np.get_include()],
    ),
]

extensions = cythonize(extensions, force=True)


setup(
    packages=find_packages(),
    ext_modules=extensions,
    include_dirs=[np.get_include()],
)
