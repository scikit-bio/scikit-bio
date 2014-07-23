#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import collections
import os
import os.path
import sys

import skbio as biolopy


def main():
    """Go on a power trip by nitpicking the scikit-bio repo.

    Attempts to find things that are wrong with the repo -- these are usually
    annoying details introduced by human error. The code goes out of its way
    to nitpick as much as possible in order to maximize the effectiveness of
    the power trip.

    Returns
    -------
    int
        Return code: 0 if there were no validation errors, 1 otherwise. Useful
        as an exit code (e.g. for use with ``sys.exit``).

    """
    root = 'skbio'
    validators = [InitValidator(), ExecPermissionValidator(),
                  GeneratedCythonValidator()]

    return_code = 0
    for validator in validators:
        success, msg = validator.validate(root)

        if not success:
            return_code = 1
            sys.stderr.write('\n'.join(msg))
            sys.stderr.write('\n\n')

    return return_code


class RepoValidator(object):
    """Abstract base class representing a repository validator.

    Subclasses must override and implement ``_validate`` (see its docstring for
    more details).

    Subclasses should also provide a ``reason``: this is a string describing
    the reason for a particular type of validation failure (see subclasses for
    examples). ``reason`` is included in the validation error message/report
    created by ``validate``.

    """
    reason = ''

    def validate(self, root):
        """Validate a directory tree recursively.

        Parameters
        ----------
        root : str
            Root directory to validate recursively.

        Returns
        -------
        tuple of (bool, list of str)
            First element is a ``bool`` indicating success status: ``True`` if
            `root` passed validation, ``False`` if there were any errors.
            Second element is a list of strings containing the validation error
            message.

        """
        invalids = []
        for root, dirs, files in os.walk(root):
            result = self._validate(root, dirs, files)
            invalids.extend(result)

        success = True
        msg = []
        if invalids:
            success = False
            msg.append(self.reason + ':')

            for invalid in invalids:
                msg.append("    %s" % invalid)

        return success, msg

    def _validate(self, root, dirs, files):
        """Validate a single directory.

        Subclasses must override and implement this method. The method is
        supplied with the three values yielded by ``os.walk``.

        Parameters
        ----------
        root : str
            Path to the current directory to be validated.
        dirs : list of str
            Directory names within `root`.
        files : list of str
            Filenames within `root`.

        Returns
        -------
        list of str
            List of filepaths or dirpaths to be considered invalid (i.e., that
            did not pass the validation checks).

        See Also
        --------
        os.walk

        """
        raise NotImplementedError("Subclasses must implement _validate.")


class InitValidator(RepoValidator):
    """Flag library code directories that are missing init files.

    This type of validation is important mainly because it is very easy to
    forget to add an __init__.py file to a new test directory. If this
    happens, nose will skip those tests unless it is run from the root of the
    source repository. Thus, the tests will be skipped if the package is
    pip-installed, e.g., as an end-user might install a release.

    Parameters
    ----------
    skip_dirs : iterable of str, optional
        Directory names to skip during validation. Defaults to skipping any
        directories named ``'data'`` or ``'__pycache__'`` (and anything
        contained within them).

    """
    reason = "Directories missing init files"

    def __init__(self, skip_dirs=None):
        if skip_dirs is None:
            skip_dirs = {'data', '__pycache__'}
        self.skip_dirs = set(skip_dirs)

    def _validate(self, root, dirs, files):
        # If any of the directories yet to be visited should be skipped, remove
        # them from ``dirs`` so that we don't visit them in a future iteration.
        # This guarantees that ``root`` is a valid directory that should not be
        # skipped (since we're doing a top-down walk).
        for skip_dir in self.skip_dirs:
            if skip_dir in dirs:
                dirs.remove(skip_dir)

        invalid_dirs = []
        if '__init__.py' not in files:
            invalid_dirs.append(root)
        return invalid_dirs


class ExecPermissionValidator(RepoValidator):
    """Flag code files that have execute permissions.

    Parameters
    ----------
    extensions : iterable of str, optional
        File extensions of files to validate. Defaults to Python, Cython, and
        C files (header and source files).

    """
    reason = "Library code with execute permissions"

    def __init__(self, extensions=None):
        if extensions is None:
            extensions = {'.py', '.pyx', '.h', '.c'}
        self.extensions = set(extensions)

    def _validate(self, root, dirs, files):
        invalid_fps = []
        for f in files:
            _, ext = os.path.splitext(f)

            if ext in self.extensions:
                fp = os.path.join(root, f)

                if os.access(fp, os.X_OK):
                    invalid_fps.append(fp)
        return invalid_fps


class GeneratedCythonValidator(RepoValidator):
    """Flag Cython files that are missing generated C files.

    Flags Cython files that aren't paired with a generated C file. The
    generated C file must be in the same directory as the Cython file, and its
    name (besides the file extension) must match. The validator also ensures
    that the generated C file is not empty.

    Parameters
    ----------
    cython_ext : str, optional
        File extension for Cython files.
    c_ext : str, optional
        File extension for generated C files.

    """
    reason = "Cython code missing generated C code"

    def __init__(self, cython_ext='.pyx', c_ext='.c'):
        self.cython_ext = cython_ext
        self.c_ext = c_ext

    def _validate(self, root, dirs, files):
        invalid_fps = []
        ext_to_base = collections.defaultdict(list)

        # Map from file extension to a list of basenames (without the
        # extension).
        for f in files:
            base, ext = os.path.splitext(f)
            ext_to_base[ext].append(base)

        # For each Cython file, try to find a matching C file. If we have a
        # match, make sure the C file isn't empty.
        for cython_base in ext_to_base[self.cython_ext]:
            cython_fp = os.path.join(root, cython_base + self.cython_ext)
            c_fp = os.path.join(root, cython_base + self.c_ext)

            if cython_base not in ext_to_base[self.c_ext]:
                invalid_fps.append(cython_fp)
            elif os.path.getsize(c_fp) <= 0:
                invalid_fps.append(cython_fp)

        return invalid_fps


if __name__ == '__main__':
    sys.exit(main())
