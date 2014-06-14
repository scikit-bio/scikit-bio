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
    validators = [TestInitValidator()]

    return_code = 0
    for validator in validators:
        success, msg = validator.validate(root)

        if not success:
            return_code = 1
            sys.stderr.write('\n'.join(msg))
            sys.stderr.write('\n')

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
            msg.append("Missing __init__.py files inside test directories:")

            for missing_init in missing_inits:
                msg.append("    %s" % missing_init)

        return success, msg


if __name__ == '__main__':
    sys.exit(main())
