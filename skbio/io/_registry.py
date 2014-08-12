# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn

from .util import open_file, _is_string_or_bytes
from skbio.io import (FormatIdentificationError, FileFormatError,
                      DuplicateRegistrationError, UnprovenFormatWarning)

_formats = {}
_identifiers = {}


def register_identifier(fmt):
    """Return a decorator for an identifier function.

    A decorator factory for identifier functions.

    An identifier function should have at least the following signature:
    ``<format_name>_identifier(fp)``. `fp` is either a filepath or
    an open fileobject.

    **The identifier must not close an open fileobject**, cleanup must be
    handled external to the identifier and is not it's concern. If it is a
    filepath, then opening and closing the file will be the responsibility of
    the identifier.

    Any additional `*args` and `**kwargs` will be passed to the identifier and
    may be used if necessary.

    The identifier **must** return an True if it believes `fh` is a given
    `fmt`. Otherwise it should return False.

    The identifier may determine membership of a file in as many or as few
    lines of the file as it deems necessary.

    Parameters
    ----------
    fmt : str
        A format name which a decorated identifier will be bound to.

    Returns
    -------
    function
        A decorator to be used on a identifer. The decorator will raise a
        ``skbio.io.DuplicateRegistrationError`` if there already exists an
        *identifier* bound to the `fmt`.

    Note
    -----
        Failure to adhere to the above interface specified for an identifier
        will result in unintended side-effects.

        The returned decorator does not mutate the decorated function in any
        way, it only adds the function to a global registry for use with
        ``skbio.io.guess_format``

    See Also
    --------
    skbio.io.guess_format

    """
    def decorator(identifier):
        if fmt in _identifiers:
            raise DuplicateRegistrationError("'%s' already has an identifier."
                                             % fmt)

        def wrapped_identifier(fp, mode='U', **kwargs):
            with open_file(fp, mode) as fh:
                orig_pos = fh.tell()
                fh.seek(0)
                result = identifier(fh, **kwargs)
                fh.seek(orig_pos)
                return result

        wrapped_identifier.__doc__ = identifier.__doc__
        wrapped_identifier.__name__ = identifier.__name__

        _identifiers[fmt] = wrapped_identifier
        return wrapped_identifier
    return decorator


def register_reader(fmt, *cls):
    """Return a decorator for a reader function.

    A decorator factory for reader functions.

    A reader function should have at least the following signature:
    ``<format_name>_to_<class_name_or_generator>(fp)``. `fp` is either a
    filepath or an open fileobject.

    **The reader must not close an open fileobject**, cleanup must be
    handled external to the reader and is not it's concern. If it is a
    filepath, then opening and closing the file will be the responsibility of
    the reader.

    Any additional `*args` and `**kwargs` will be passed to the reader and may
    be used if necessary.

    The reader **must** return an instance of `cls` if `cls` is not None.
    Otherwise the reader must return a generator. The generator need not deal
    with closing the `fh` this is the responsibility of the caller and is
    handled for you in ``skbio.io.read``.


    Parameters
    ----------
    fmt : str
        A format name which a decorated reader will be bound to.
    cls : type, optional
        Positional argument.
        The class which a decorated reader will be bound to. If not provided
        or is None, the decorated reader will be bound as a generator.
        Default is None.

    Returns
    -------
    function
        A decorator to be used on a reader. The decorator will raise a
        ``skbio.io.DuplicateRegistrationError`` if there already exists a
        *reader* bound to the same permutation of `fmt` and `cls`.

    Raises
    ------
    TypeError

    Note
    -----
        Failure to adhere to the above interface specified for a reader will
        result in unintended side-effects.

        The returned decorator does not mutate the decorated function in any
        way, it only adds the function to a global registry for use with
        ``skbio.io.read``

    See Also
    --------
    skbio.io.read

    """
    arg_len = len(cls)
    if arg_len > 1:
        raise TypeError("register_reader takes 1 or 2 arguments (%d given)"
                        % (arg_len + 1))

    cls = None if arg_len == 0 else cls[0]

    def decorator(reader):
        format_class = _formats.setdefault(fmt, {}).setdefault(cls, {})

        if 'reader' in format_class:
            raise DuplicateRegistrationError("'%s' already has a %s for %s."
                                             % (fmt, 'reader', cls.__name__))

        if cls is None:
            def wrapped_reader(fp, mode='U', mutate_fh=False, **kwargs):
                with open_file(fp, mode) as fh:
                    generator = reader(fh, **kwargs)
                    if not mutate_fh and not _is_string_or_bytes(fp):
                        orig_pos = fh.tell()
                        read_pos = orig_pos
                        try:
                            while True:
                                orig_pos = fh.tell()

                                fh.seek(read_pos)
                                next_result = next(generator)
                                read_pos = fh.tell()

                                fh.seek(orig_pos)

                                yield next_result
                        finally:
                            fh.seek(orig_pos)
                    else:
                        while True:
                            yield next(generator)

        else:
            def wrapped_reader(fp, mode='U', mutate_fh=False, **kwargs):
                with open_file(fp, mode) as fh:
                    orig_pos = fh.tell()
                    result = reader(fh, **kwargs)
                    if not mutate_fh:
                        fh.seek(orig_pos)
                    return result

        wrapped_reader.__doc__ = reader.__doc__
        wrapped_reader.__name__ = reader.__name__

        format_class['reader'] = wrapped_reader
        return wrapped_reader
    return decorator


def register_writer(fmt, *cls):
    """Return a decorator for a writer function.

    A decorator factory for writer functions.

    A writer function should have at least the following signature:
    ``<class_name_or_generator>_to_<format_name>(obj, fp)`` where `obj` will
    be either an instance of <class_name> or a generator that is *identical* to
    the result of calling ``get_reader(<format>, None)``. `fp` is either a
    filepath or an open fileobject.

    **The writer must not close an open fileobject**, cleanup must be
    handled external to the writer and is not it's concern. If it is a
    filepath, then opening and closing the file will be the responsibility of
    the writer.

    Any additional `*args` and `**kwargs` will be passed to the writer and may
    be used if necessary.

    The writer must not return a value. Instead it should only mutate the `fh`
    in a way consistent with it's purpose.

    If the writer accepts a generator, it should exhaust the generator to
    ensure that the potentially open fileobject backing said generator is
    closed.

    Parameters
    ----------
    fmt : str
        A format name which a decorated writer will be bound to.
    cls : type, optional
        Positional argument.
        The class which a decorated writer will be bound to. If not provided
        or is None, the decorated writer will be bound as a generator.
        Default is None.

    Returns
    -------
    function
        A decorator to be used on a writer. The decorator will raise a
        ``skbio.io.DuplicateRegistrationError`` if there already exists a
        *writer* bound to the same permutation of `fmt` and `cls`.

    Raises
    ------
    TypeError

    Note
    -----
        Failure to adhere to the above interface specified for a writer will
        result in unintended side-effects.

        The returned decorator does not mutate the decorated function in any
        way, it only adds the function to a global registry for use with
        ``skbio.io.write``

    See Also
    --------
    skbio.io.write
    skbio.io.get_reader

    """
    arg_len = len(cls)
    if arg_len > 1:
        raise TypeError("register_writer takes 1 or 2 arguments (%d given)"
                        % (arg_len + 1))

    cls = None if arg_len == 0 else cls[0]

    def decorator(writer):
        format_class = _formats.setdefault(fmt, {}).setdefault(cls, {})

        if 'writer' in format_class:
            raise DuplicateRegistrationError("'%s' already has a %s for %s."
                                             % (fmt, 'writer', cls.__name__))

        def wrapped_writer(obj, fp, mode='w', **kwargs):
            with open_file(fp, mode) as fh:
                writer(obj, fh, **kwargs)

        wrapped_writer.__doc__ = writer.__doc__
        wrapped_writer.__name__ = writer.__name__

        format_class['writer'] = wrapped_writer
        return wrapped_writer
    return decorator


def list_read_formats(cls):
    """Return a list of available read formats for a given `cls` type.

    Parameters
    ----------
    cls : type
        The class which will be used to determine what read formats exist for
        an instance of `cls`.

    Returns
    -------
    list
        A list of available read formats for an instance of `cls`. List may be
        empty.

    See Also
    --------
    skbio.io.register_reader

    """
    return _rw_list_formats('reader', cls)


def list_write_formats(cls):
    """Return a list of available write formats for a given `cls` instance.

    Parameters
    ----------
    cls : type
        The class which will be used to determine what write formats exist for
        an instance of `cls`.

    Returns
    -------
    list
        A list of available write formats for an instance of `cls`. List may be
        empty.

    See Also
    --------
    skbio.io.register_writer

    """
    return _rw_list_formats('writer', cls)


def _rw_list_formats(name, cls):
    formats = []
    for fmt in _formats:
        if cls in _formats[fmt]:
            if name in _formats[fmt][cls]:
                formats.append(fmt)
    return formats


def get_identifier(fmt):
    """Return an identifier for a format.

    Parameters
    ----------
    fmt : str
        A format string which has a registered identifier.

    Returns
    -------
    function or None
        Returns an identifier function if one exists for the given `fmt`.
        Otherwise it will return None.

    See Also
    --------
    skbio.io.register_identifier

    """

    if fmt in _identifiers:
        return _identifiers[fmt]
    return None


def get_reader(fmt, *cls):
    """Return a reader for a format.

    Parameters
    ----------
    fmt : str
        A registered format string.
    cls : type, optional
        Positional argument.
        The class which the reader will return an instance of. If not provided
        or is None, the reader will return a generator.
        Default is None.

    Returns
    -------
    function or None
        Returns a reader function if one exists for a given `fmt` and `cls`.
        Otherwise it will return None.

    See Also
    --------
    skbio.io.register_reader

    """

    return _rw_getter('reader', fmt, *cls)


def get_writer(fmt, *cls):
    """Return a writer for a format.

    Parameters
    ----------
    fmt : str
        A registered format string.
    cls : type, optional
        Positional argument.
        The class which the writer will expect an instance of. If not provided
        or is None, the writer will expect a generator that identical to what
        is returned by ``get_reader(<some_format>, None)``.
        Default is None.

    Returns
    -------
    function or None
        Returns a writer function if one exists for a given `fmt` and `cls`.
        Otherwise it will return None.

    See Also
    --------
    skbio.io.register_writer

    """

    return _rw_getter('writer', fmt, *cls)


def _rw_getter(name, fmt, *args):
    cls = None
    arg_len = len(args)
    if arg_len > 1:
        raise TypeError("get_%s takes 1 or 2 arguments (%d given)"
                        % (name, arg_len + 1))
    if arg_len == 1:
        cls = args[0]

    if fmt in _formats:
        if cls in _formats[fmt]:
            if name in _formats[fmt][cls]:
                return _formats[fmt][cls][name]
    return None


def guess_format(fp, cls=None):
    """Attempt to guess the format of a file and return format str.

    Parameters
    ----------
    fp : filepath or fileobject
        The provided file to guess the format of. Filepaths are automatically
        closed; fileobjects are the responsibility of the caller.
    cls : type, optional
        A provided class that restricts the search for the format. Only formats
        which have a registered reader or writer for the given `cls` will be
        tested.
        Default is None.

    Returns
    -------
    str
        A registered format name.

    Raises
    ------
    FormatIdentificationError

    Note
    -----
        If a fileobject is provided, the current read offset will be reset.

        If the file is 'claimed' by multiple identifiers, or no identifier
        'claims' the file, an ``skbio.io.FormatIdentificationError`` will be
        raised.

    See Also
    --------
    skbio.io.register_identifier

    """
    possibles = []
    for fmt in _identifiers:
        if cls is not None and (fmt not in _formats or
                                cls not in _formats[fmt]):
            continue
        format_identifier = _identifiers[fmt]
        if format_identifier(fp, mode='U'):
            possibles.append(fmt)

    if not possibles:
        raise FormatIdentificationError("Cannot guess the format for %s."
                                        % str(fp))
    if len(possibles) > 1:
        raise FormatIdentificationError("File format is ambiguous, may be"
                                        " one of %s." % str(possibles))

    return possibles[0]


def read(fp, format=None, into=None, verify=True, mode='U', *args, **kwargs):
    """Generalized read function: multiplex read functionality in skbio.

    This function is able to reference and execute all *registered* read
    operations in skbio.

    Parameters
    ----------
    fp : filepath or fileobject
        The location to read the given `format` `into`. Filepaths are
        automatically closed when read; fileobjects are the responsibility
        of the caller. In the case of a generator, a filepath will be closed
        when ``StopIteration`` is raised; fileobjects are still the
        responsibility of the caller.
    format : str, optional
        The format must be a reigstered format name with a reader for the given
        `into` class. If a `format` is not provided or is None, all registered
        identifiers for the provied `into` class will be evaluated to attempt
        to guess the format. Will raise an
        ``skbio.io.FormatIdentificationError`` if it is unable to guess.
        Default is None.
    into : type, optional
        A class which has a registered reader for a given `format`. If `into`
        is not provided or is None, read will return a generator.
        Default is None.
    verify : bool, optional
        Whether or not to confirm the format of a file if `format` is provided.
        Will raise a ``skbio.io.UnprovenFormatWarning`` if the identifier of
        `format` returns False.
        Default is True.
    mode : str, optional
        The read mode. This is passed to `open(fp, mode)` internally.
        Default is 'U'.
    args : tuple, optional
        Will be passed directly to the appropriate reader.
    kwargs : dict, optional
        Will be passed directly to the appropriate reader.

    Returns
    -------
    object or generator
        If `into` is not None, an instance of the `into` class will be
        provided with internal state consistent with the provided file.
        If `into` is None, a generator will be returned.

    Raises
    ------
    ValueError
    skbio.io.FileFormatError
    skbio.io.FormatIdentificationError
    skbio.io.UnprovenFormatWarning

    Note
    -----
        Will raise a ``ValueError`` if `format` and `into` are both be None.

    See Also
    --------
    skbio.io.register_reader
    skbio.io.register_identifier

    """
    if format is None and into is None:
        raise ValueError("`format` and `into` cannot both be None.")

    if format is None:
        format = guess_format(fp, cls=into)
    elif verify:
        identifier = get_identifier(format)
        if identifier is not None:
            if not identifier(fp):
                warn("%s could not be positively identified as a %s file." %
                     (str(fp), format),
                     UnprovenFormatWarning)

    reader = get_reader(format, into)
    if reader is None:
        raise FileFormatError("Cannot read %s into %s, no reader found."
                              % (format, into.__name__
                                 if into is not None
                                 else 'generator'))

    return reader(fp, mode=mode, **kwargs)


def write(obj, format=None, into=None, mode='w', *args, **kwargs):
    """Generalized write function: multiplex write functionality in skbio.

    This function is able to reference and execute all *registered* write
    operations in skbio.

    Parameters
    ----------
    obj : object
        The object must have a registered writer for a provided `format`.
    format : str
        The format must be a reigstered format name with a writer for the given
        `obj`
    into : filepath or fileobject
        The location to write the given `format` from `obj` into. Filepaths are
        automatically closed when written; fileobjects are the responsibility
        of the caller.
    mode : str, optional
        The write mode. This is passed to `open(fp, mode)` internally.
        Default is 'w'.
    args : tuple, optional
        Will be passed directly to the appropriate writer.
    kwargs : dict, optional
        Will be passed directly to the appropriate writer.

    Raises
    ------
    ValueError
    skbio.io.FileFormatError

    See Also
    --------
    skbio.io.register_writer

    """
    if format is None:
        raise ValueError("Must specify a `format` to write out as.")
    if into is None:
        raise ValueError("Must provide a filepath or filehandle for `into`")

    writer = get_writer(format, obj.__class__)
    if writer is None:
        raise FileFormatError("Cannot write %s into %s, no writer found."
                              % (format, str(into)))

    writer(obj, into, mode=mode, **kwargs)


@register_identifier('<empty file>')
def empty_file_identifier(fh):
    for line in fh:
        if line.strip():
            return False
    return True
