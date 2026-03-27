# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import types

import skbio.io


class Read:
    """A descriptor class to generate read methods for scikit-bio objects."""

    def __get__(self, instance, cls):
        if "_read_method" not in cls.__dict__:

            def _read_method(file, format=None, **kwargs):
                return skbio.io.read(file, into=cls, format=format, **kwargs)

            _read_method.__doc__ = self._make_docstring(cls)
            cls._read_method = _read_method
        return cls._read_method

    def _make_docstring(self, cls):
        name, supported_fmts, default, see = _docstring_vars(cls, "read")
        constructor_param = ""
        if name == "TabularMSA":
            constructor_param = """constructor : subclass of ``GrammaredSequence``
    Sequence constructor (e.g., ``DNA``, ``RNA``, ``Protein``). This parameter
    is required when reading ``TabularMSA`` from FASTA, FASTQ, Clustal,
    PHYLIP, and Stockholm formats.
"""
        return f"""Create a new ``{name}`` instance from a file.

This is a convenience method for :func:`skbio.io.registry.read`. For more information
about the I/O system in scikit-bio, please see :mod:`skbio.io`.

{supported_fmts}

Parameters
----------
file : openable (filepath, URL, filehandle, etc.)
    The location to read the given `format` into. Something that is understood by
    :func:`skbio.io.util.open`. Filehandles are not automatically closed, it is the
    responsibility of the caller.
format : str, optional
    The format of the file. The format must be a format name with a reader for
    ``{name}``. If None, the format will be inferred.
{constructor_param}\
kwargs : dict, optional
    Additional arguments passed to :func:`skbio.io.registry.read()` and the reader for
    ``{name}``.

Returns
-------
``{name}``
    A new instance.

See Also
--------
write
skbio.io.registry.read
skbio.io.util.open
{see}

"""


class Write:
    """A descriptor class to generate write methods for scikit-bio objects."""

    def __get__(self, instance, cls):
        """This gets called when any skbio object accesses its ``write`` attribute."""
        if instance is None:
            if "_write_method" not in cls.__dict__:
                cls._write_method = self._generate_write_method(cls)
            return cls._write_method
        if "_write_method" not in instance.__dict__:
            instance._write_method = types.MethodType(
                self._generate_write_method(cls), instance
            )
        return instance._write_method

    def _generate_write_method(self, cls):
        def _write_method(self, file, format=None, **kwargs):
            if format is None:
                if hasattr(cls, "default_write_format"):
                    format = cls.default_write_format
                else:
                    raise ValueError(f"{cls.__name__} has no default write format.")
            return skbio.io.write(self, into=file, format=format, **kwargs)

        _write_method.__doc__ = self._make_docstring(cls)

        return _write_method

    def _make_docstring(self, cls):
        name, supported_fmts, default, see = _docstring_vars(cls, "write")
        return f"""Write an instance of ``{name}`` to a file.

This is a convenience method for :func:`skbio.io.registry.write()`. For more
information about the I/O system in scikit-bio, please see :mod:`skbio.io`.

{supported_fmts}

Parameters
----------
file : openable (filepath, filehandle, etc.)
    The location to write the given `format` into. Something that is understood by
    :func:`skbio.io.util.open()`. Filehandles are not automatically closed, it is the
    responsibility of the caller.
format : str, optional
    The format to write the ``{name}`` object as. The format must be a registered
    format name with a writer for ``{name}``. Default is ``'{default}'``.
kwargs : dict, optional
    Additional arguments passed to the writer for ``{name}``.

See Also
--------
read
skbio.io.registry.write
skbio.io.util.open
{see}

"""


def _docstring_vars(cls, func):
    """Generate variables for dynamically generated docstrings."""
    from skbio.io.registry import io_registry

    if func == "write":
        formats = io_registry.list_write_formats(cls)
    elif func == "read":
        formats = io_registry.list_read_formats(cls)
    else:
        raise ValueError("'func' parameter must be 'read' or 'write'.")

    imports = io_registry._import_paths(formats)
    formats = io_registry._formats_for_docs(formats, imports)
    if formats:
        supported_fmts = f"Supported file formats include:\n\n{formats}"
    else:
        supported_fmts = ""
    name = cls.__name__
    default = getattr(cls, "default_write_format", "None")
    see = "\n".join(imports)

    return name, supported_fmts, default, see
