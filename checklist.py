#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import os
import os.path
import subprocess
import sys
import ast
import tokenize
import warnings

import dateutil.parser

if sys.version_info.major != 3:
    sys.exit("scikit-bio can only be used with Python 3. You are currently "
             "running Python %d." % sys.version_info.major)


class ChecklistWarning(Warning):
    """General warning class for warnings raised by checklist.py."""
    pass


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
    validators = [InitValidator(), CopyrightHeadersValidator(),
                  ExecPermissionValidator(), GeneratedCythonValidator(),
                  APIRegressionValidator()]

    return_code = 0
    for validator in validators:
        success, msg = validator.validate(root)

        if not success:
            return_code = 1
            sys.stderr.write('\n'.join(msg))
            sys.stderr.write('\n\n')

    return return_code


class RepoValidator:
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
            msg.append(self.reason)

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

    def _system_call(self, cmd):
        """Issue a system call, returning stdout, stderr, and return value.

        This code was taken from verman's
        ``verman.Version.verman_system_call``. See licenses/verman.txt and
        https://github.com/biocore/verman for more details.

        """
        proc = subprocess.Popen(cmd, shell=True, universal_newlines=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # communicate pulls all stdout/stderr from the PIPEs to
        # avoid blocking -- don't remove this line!
        stdout, stderr = proc.communicate()
        return_value = proc.returncode
        return stdout, stderr, return_value


class CopyrightHeadersValidator(RepoValidator):
    """Flag library files with non-standard copyright headers

    See the current standard for scikit-bio's copyright headers at
    ``http://scikit-bio.org/docs/latest/development/new_module.html``

    Individual files are ignored if the first line in the file is exactly:

    # checklist.py:CopyrightHeadersValidator IGNORE

    If a file is ignored, a ``ChecklistWarning`` is raised.

    Parameters
    ----------
    skip_dirs : iterable of str, optional
        Directory names to skip during validation. Defaults to skipping any
        directories named ``'data'`` or ``'__pycache__'`` (and anything
        contained within them).

    """

    reason = ("Files non-conforming to standard headers as described in\n"
              "http://scikit-bio.org/docs/latest/development/new_module.html:")

    COPYRIGHT_HEADER = """\
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
"""

    def __init__(self, skip_dirs=None):
        if skip_dirs is None:
            skip_dirs = {'data', '__pycache__'}
        self.skip_dirs = set(skip_dirs)

    def _validate(self, root, dirs, files):
        for skip_dir in self.skip_dirs:
            if skip_dir in dirs:
                dirs.remove(skip_dir)

        invalid_files = []
        for _file in files:
            if not _file.endswith('.py'):
                continue

            pos = 0
            filepath = os.path.join(root, _file)
            f = open(filepath)

            first_line = f.readline().rstrip('\n')
            if first_line == '# checklist.py:CopyrightHeadersValidator IGNORE':
                warnings.warn(
                    "File %s has IGNORE directive. Ignoring scikit-bio "
                    "copyright header validation." % filepath,
                    ChecklistWarning)
                continue

            f.seek(0)
            tokens = list(tokenize.generate_tokens(f.readline))

            # A module docstring is fully described using just two tokens: the
            # main string and a terminating newline. By convention, however, it
            # is always followed by a newline, and thus we advance by three
            # positions to get to the next logical line.
            if tokens[pos][0] == tokenize.STRING:
                pos += 3
            # copyright header consists of 7 lines, and by discussion in
            # preceding comment, spans through 14 tokens.
            cheader = ''.join(map(lambda x: x[1], tokens[pos:pos + 14]))
            # Ensure that there is no blank line at the end of the file
            if (cheader != self.COPYRIGHT_HEADER or
                    (tokens[pos + 14][0] != tokenize.NL and
                     tokens[pos + 14][0] != tokenize.ENDMARKER)):
                invalid_files.append(f.name)
            f.close()

        return invalid_files


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
    reason = "Directories missing init files:"

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
    reason = "Library code with execute permissions:"

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
    """Flag Cython files that have missing or outdated generated C files.

    Flags Cython files that aren't paired with an up-to-date generated C file.
    The generated C file must be in the same directory as the Cython file, and
    its name (besides the file extension) must match. The validator also
    ensures that the generated C file is not empty and that it was generated at
    the same time or later than the Cython file's timestamp.

    Parameters
    ----------
    cython_ext : str, optional
        File extension for Cython files.
    c_ext : str, optional
        File extension for generated C files.

    """
    reason = "Cython code with missing or outdated generated C code:"

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
        # match, make sure the C file isn't empty and that it was generated at
        # the same time or later than the Cython file.
        for cython_base in ext_to_base[self.cython_ext]:
            cython_fp = os.path.join(root, cython_base + self.cython_ext)
            c_fp = os.path.join(root, cython_base + self.c_ext)

            if cython_base not in ext_to_base[self.c_ext]:
                invalid_fps.append(cython_fp)
            elif os.path.getsize(c_fp) <= 0:
                invalid_fps.append(cython_fp)
            else:
                cython_ts = self._get_timestamp(cython_fp)
                c_ts = self._get_timestamp(c_fp)

                if c_ts < cython_ts:
                    invalid_fps.append(cython_fp)

        return invalid_fps

    def _get_timestamp(self, fp):
        cmd = 'git log -1 --format="%%ad" -- %s' % fp
        stdout, stderr, retval = self._system_call(cmd)

        if retval != 0:
            raise RuntimeError("Could not execute 'git log' command to "
                               "determine file timestamp.")
        return dateutil.parser.parse(stdout.strip())


class APIRegressionValidator(RepoValidator):
    """Flag tests that import from a non-minimized subpackage hierarchy.

    Flags tests that aren't imported from a minimally deep API target. (e.g.
    skbio.TabularMSA vs skbio.alignment.TabularMSA). This should prevent
    accidental regression in our API because tests will fail if any alias is
    removed, and this checklist will fail if any test doesn't import from the
    least deep API target.

    """
    reason = ("The following tests import `A` but should import `B`"
              " (file: A => B):")

    def __init__(self):
        self._imports = {}

    def _validate(self, root, dirs, files):
        errors = []
        test_imports = []
        for file in files:
            current_fp = os.path.join(root, file)
            package, ext = os.path.splitext(current_fp)
            if ext == ".py":
                imports = self._parse_file(current_fp, root)
                if os.path.split(root)[1] == "tests":
                    test_imports.append((current_fp, imports))

                temp = package.split(os.sep)
                # Remove the __init__ if it is a directory import
                if temp[-1] == "__init__":
                    temp = temp[:-1]
                    package = ".".join(temp)
                    self._add_imports(imports, package)
        for fp, imports in test_imports:
            for import_ in imports:
                substitute = self._minimal_import(import_)
                if substitute is not None:
                    errors.append("%s: %s => %s" %
                                  (fp, import_, substitute))

        return errors

    def _add_imports(self, imports, package):
        """Add the minimum depth import to our collection."""
        for import_ in imports:
            value = import_
            # The actual object imported will be the key.
            key = import_.split(".")[-1]
            # If package importing the behavior is shorter than its import:
            if len(package.split('.')) + 1 < len(import_.split('.')):
                value = ".".join([package, key])

            if key in self._imports:
                sub = self._imports[key]
                if len(sub.split('.')) > len(value.split('.')):
                    self._imports[key] = value
            else:
                self._imports[key] = value

    def _minimal_import(self, import_):
        """Given an normalized import, return a shorter substitute or None."""
        key = import_.split(".")[-1]
        if key not in self._imports:
            return None
        substitute = self._imports[key]
        if len(substitute.split('.')) == len(import_.split('.')):
            return None
        else:
            return substitute

    def _parse_file(self, fp, root):
        """Parse a file and return all normalized skbio imports."""
        imports = []
        with open(fp, 'U') as f:
            # Read the file and run it through AST
            source = ast.parse(f.read())
            # Get each top-level element, this is where API imports should be.
            for node in ast.iter_child_nodes(source):
                if isinstance(node, ast.Import):
                    # Standard imports are easy, just get the names from the
                    # ast.Alias list `node.names`
                    imports += [x.name for x in node.names]
                elif isinstance(node, ast.ImportFrom):
                    prefix = ""
                    # Relative import handling.
                    if node.level > 0:
                        prefix = root
                        extra = node.level - 1
                        while(extra > 0):
                            # Keep dropping...
                            prefix = os.path.split(prefix)[0]
                            extra -= 1
                        # We need this in '.' form not '/'
                        prefix = prefix.replace(os.sep, ".") + "."
                    # Prefix should be empty unless node.level > 0
                    if node.module is None:
                        node.module = ""
                    imports += [".".join([prefix + node.module, x.name])
                                for x in node.names]
        skbio_imports = []
        for import_ in imports:
            # Filter by skbio
            if import_.split(".")[0] == "skbio":
                skbio_imports.append(import_)
        return skbio_imports


if __name__ == '__main__':
    sys.exit(main())
