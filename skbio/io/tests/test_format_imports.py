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
from skbio.io import _exception


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

    def test_all_io_exceptions_importable(self):
        """Verify all I/O exceptions are importable and functional.

        This test covers ALL exception classes exported from skbio.io,
        not just those ending with 'FormatError'. This includes:
        - IOSourceError
        - FileFormatError (base class)
        - All format-specific exceptions (*FormatError)
        """
        # All exception classes exported from skbio.io
        io_exceptions = {
            name: obj
            for name, obj in vars(skbio.io).items()
            if inspect.isclass(obj) and issubclass(obj, Exception)
        }

        self.assertGreater(
            len(io_exceptions), 0,
            "No exception classes found in skbio.io"
        )

        failed_exceptions = []
        for name, exception_cls in io_exceptions.items():
            try:
                # Test instantiation with a message
                ex = exception_cls("test message")
                # Test raising and catching
                try:
                    raise ex
                except exception_cls as caught:
                    # Verify the message is preserved
                    self.assertIn("test message", str(caught))
            except AssertionError:
                raise
            except Exception as e:
                failed_exceptions.append((name, str(e)))

        if failed_exceptions:
            self.fail(
                f"Failed to instantiate or raise exceptions: {failed_exceptions}"
            )

    def test_exception_module_exports_match_definitions(self):
        """Verify public exceptions from _exception.py are exported in __all__.

        This ensures that exception classes intended for public use
        are properly exported from skbio.io.
        """
        # Get all public exception classes defined in _exception module
        defined_exceptions = {
            name
            for name, obj in vars(_exception).items()
            if inspect.isclass(obj)
            and issubclass(obj, Exception)
            and not name.startswith("_")
        }

        # Get exceptions listed in skbio.io.__all__
        exported_exceptions = {
            name for name in skbio.io.__all__
            if name.endswith("Error")
        }

        # Internal exceptions that are intentionally not exported
        internal_exceptions = {"InvalidRegistrationError", "DuplicateRegistrationError"}

        # Check that all non-internal exceptions are exported
        expected_exports = defined_exceptions - internal_exceptions
        missing_exports = expected_exports - exported_exceptions

        self.assertEqual(
            missing_exports, set(),
            f"Exceptions defined but not exported in __all__: {missing_exports}"
        )

    def test_exception_inheritance_hierarchy(self):
        """Verify format-specific exceptions inherit from FileFormatError.

        This ensures that users can catch all format-related errors
        by catching FileFormatError.
        """
        format_exceptions = {
            name: obj
            for name, obj in vars(skbio.io).items()
            if inspect.isclass(obj)
            and issubclass(obj, Exception)
            and name.endswith("FormatError")
            and name != "FileFormatError"
        }

        # These exceptions have different base classes by design
        exceptions_with_different_base = {"BIOMFormatError", "EmbedFormatError"}

        for name, exception_cls in format_exceptions.items():
            if name not in exceptions_with_different_base:
                self.assertTrue(
                    issubclass(exception_cls, skbio.io.FileFormatError),
                    f"{name} should inherit from FileFormatError"
                )


if __name__ == "__main__":
    main()
