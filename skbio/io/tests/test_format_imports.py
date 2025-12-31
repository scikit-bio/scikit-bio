# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import inspect
import pkgutil
from unittest import TestCase, main

import skbio.io


class ImportModuleTests(TestCase):
    def test_format_modules_imported_in_init(self):
        """Verify all format modules are imported in __init__.py"""
        init_file = os.path.join(os.path.dirname(skbio.io.__file__), '__init__.py')
        
        with open(init_file, 'r') as f:
            init_content = f.read()
        
        fmt_dir = os.path.dirname(skbio.io.format.__file__)
        
        # Skip these special cases
        skip_modules = {'emptyfile', '_test'}
        
        for _, name, ispkg in pkgutil.iter_modules([fmt_dir]):
            if not ispkg and not name.startswith("_") and name not in skip_modules:
                expected_import = f'import_module("skbio.io.format.{name}")'
                self.assertIn(
                    expected_import,
                    init_content,
                    f"\n\n\nModule '{name}' not imported in skbio.io.__init__.py"
                )

    def test_format_exceptions_importable(self):
        """Verify all format-specific exceptions are importable from skbio.io"""
        # Get all format exception classes defined in the io module
        io_exceptions = {
            name: obj
            for name, obj in vars(skbio.io).items()
            if inspect.isclass(obj)
            and issubclass(obj, Exception)
            and name.endswith("FormatError")
            and name != "FileFormatError"
        }
        
        failed_exceptions = []
        for name, exception_cls in io_exceptions.items():
            try:
                # Test instantiation
                ex = exception_cls()
                # Test raising
                raise ex
            except exception_cls:
                # This is expected - exception raised and caught correctly
                pass
            except Exception as e:
                # Any other exception means something is wrong
                failed_exceptions.append((name, str(e)))
        
        if failed_exceptions:
            self.fail(f"Failed to instantiate or raise exceptions: {failed_exceptions}")


if __name__ == "__main__":
    main()
