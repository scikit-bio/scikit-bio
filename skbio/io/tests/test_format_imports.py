# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import importlib
import inspect
import pkgutil
from unittest import TestCase, main

import skbio.io


class ImportModuleTests(TestCase):
    def setUp(self):
        # Need this to handle cases where module name does not match format name.
        self.fmt_mod_map = {
            "blast6": "blast+6",
            "blast7": "blast+7",
            "emptyfile": "<emptyfile>",
        }

        self.registered_fmts = [
            fmt.name for fmt in skbio.io.io_registry._binary_formats.values()
        ] + [fmt.name for fmt in skbio.io.io_registry._text_formats.values()]

    def test_format_imports(self):
        fmt_dir = os.path.dirname(skbio.io.format.__file__)
        fmt_modules = []

        for _, name, ispkg in pkgutil.iter_modules([fmt_dir]):
            if not ispkg and not name.startswith("_"):
                fmt_modules.append(name)

        failed_imports = []
        for mod_name in fmt_modules:
            try:
                full_name = f"skbio.io.format.{mod_name}"
                module = importlib.import_module(full_name)
                format_name = self.fmt_mod_map.get(mod_name, mod_name)
                self.assertIn(
                    format_name,
                    self.registered_fmts,
                    f"Format '{format_name} (from module '{mod_name}') not registered in io_registry",
                )
            except Exception as e:
                failed_imports.append((mod_name, str(e)))

        if failed_imports:
            self.fail(f"Failed to import format modules: {failed_imports}")

    def test_format_exceptions(self):
        # Create dictionary of all existing format errors defined in the io module.
        # name is a string which contains the name of the error, whereas obj is
        # the actual object itself.
        io_exceptions = {
            name: obj
            for name, obj in vars(skbio.io).items()
            if inspect.isclass(obj)
            and issubclass(obj, Exception)
            and name.endswith("FormatError")
            and name != "FileFormatError"
        }

        for name, exception_cls in io_exceptions.items():
            # First, try to instantiate and raise the exception object itself.
            # If during these two steps (instantiation and raising), any other type
            # of exception is raised besides the one we want to be raised
            # (exception_cls), the test will fail.
            try:
                ex = exception_cls()
                raise ex
            except exception_cls:
                # Pass if the exception_cls is raised, this is exactly what we are
                # testing for.
                pass
            except Exception as e:
                # Fail if any other type of exception is raised.
                self.fail(f"Failed to instantiate or raise {name}: {e}")


if __name__ == "__main__":
    main()
