#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import os
import os.path
import sys


def main():
    root = 'skbio'
    validators = [TestInitValidator(), ExecPermissionValidator()]

    return_code = 0
    for validator in validators:
        success, msg = validator.validate(root)

        if not success:
            return_code = 1
            sys.stderr.write('\n'.join(msg))
            sys.stderr.write('\n\n')

    return return_code


class Validator(object):
    def validate(self, root):
        raise NotImplementedError("Validator subclasses must implement "
                                  "validate.")

class TestInitValidator(Validator):
    def __init__(self, test_dir_names=('test', 'tests'),
                 init_name='__init__.py'):
        self.test_dir_names = test_dir_names
        self.init_name = init_name

    def validate(self, root):
        missing_inits = []
        for root, dirs, files in os.walk(root):
            if (os.path.basename(root) in self.test_dir_names and
                self.init_name not in files):
                    missing_inits.append(root)

        success = True
        msg = []
        if missing_inits:
            success = False
            msg.append("Missing %s files inside test directories:" %
                       self.init_name)

            for missing_init in missing_inits:
                msg.append("    %s" % missing_init)

        return success, msg


class ExecPermissionValidator(Validator):
    def __init__(self, extensions=('.py', '.pyx', '.h', '.c')):
        self.extensions = extensions

    def validate(self, root):
        invalid_perms = []
        for root, dirs, files in os.walk(root):
            for file_ in files:
                if os.path.splitext(file_)[1] in self.extensions:
                    fp = os.path.join(root, file_)

                    if os.access(fp, os.X_OK):
                        invalid_perms.append(fp)

        success = True
        msg = []
        if invalid_perms:
            success = False
            msg.append("Library code with execute permissions:")

            for invalid_perm in invalid_perms:
                msg.append("    %s" % invalid_perm)

        return success, msg


if __name__ == '__main__':
    sys.exit(main())
