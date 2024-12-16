# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import wraps
from collections import namedtuple

from ._exception import OverrideError
from ._docstring import (
    _note_into_doc,
    _note_into_doc_param,
    _deprecation_note,
    _renaming_note,
)
from ._warning import _warn_deprecated, _warn_renamed, _warn_param_renamed


def overrides(interface_class):
    """Indicate that a member is being overridden from a specific parent class.

    Decorator for class-level members. Used to indicate that a member is being
    overridden from a specific parent class. If the member does not have a docstring,
    it will pull one from the parent class. When chaining decorators, this should be
    first as it is relatively nondestructive.

    Parameters
    ----------
    interface_class : class
        The class which has a member overridden by the decorated member.

    Returns
    -------
    function
        The function is not changed or replaced.

    Raises
    ------
    OverrideError
        If the `interface_class` does not possess a member of the same name
        as the decorated member.

    """
    # Adapted from http://stackoverflow.com/a/8313042/579416

    def overrider(method):
        if method.__name__ not in dir(interface_class):
            raise OverrideError(
                f"{method.__name__} is not present in parent "
                f"class: {interface_class.__name__}."
            )
        backup = classproperty.__get__
        classproperty.__get__ = lambda x, y, z: x
        if method.__doc__ is None:
            method.__doc__ = getattr(interface_class, method.__name__).__doc__
        classproperty.__get__ = backup
        return method

    return overrider


class classproperty(property):
    """Decorator for class-level properties.

    Supports read access only. The property will be read-only within an
    instance. However, the property can always be redefined on the class, since
    Python classes are mutable.

    Parameters
    ----------
    func : function
        Method to make a class property.

    Returns
    -------
    property
        Decorated method.

    Raises
    ------
    AttributeError
        If the property is set on an instance.

    """

    def __init__(self, func):
        name = func.__name__
        doc = func.__doc__
        super(classproperty, self).__init__(classmethod(func))
        self.__name__ = name
        self.__doc__ = doc

    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()

    def __set__(self, obj, value):
        raise AttributeError("can't set attribute")


class classonlymethod(classmethod):
    """Just like `classmethod`, but it can't be called on an instance."""

    def __get__(self, obj, cls=None):
        if obj is not None:
            raise TypeError(
                f"Class-only method called on an instance. Use "
                f"{cls.__name__}.{self.__func__.__name__} "
                "instead."
            )
        return super().__get__(obj, cls)


def deprecated(ver, msg=None, append=True):
    """Mark a function or a method as deprecated.

    Parameters
    ----------
    ver : str
        Version when deprecation became effective.
    msg : str, optional
        A custom warning message.
    append : bool, optional
        Append the custom message to the end of the default message (True, default),
        or replace the entire default message with the custom message (False).

    Notes
    -----
    This decorator has two effects: 1) A note will be added to the docstring of the
    function below the title line to indicate the deprecated status. 2) When the
    function is called at the first time, a deprecation warning will be displayed.

    """

    def decorator(func):
        note = _deprecation_note(ver, msg)
        func.__doc__ = _note_into_doc(note, func.__doc__)

        @wraps(func)
        def wrapper(*args, **kwargs):
            _warn_deprecated(func, ver, msg, append)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def aliased(name, ver=None, warn=False):
    """Create an alias for a function or method.

    Parameters
    ----------
    name : str
        Alias name of the function or method.
    ver : str, optional
        Version when alias was created.
    warn : bool, optional
        If True, raise a deprecation warning when alias is called.

    See Also
    --------
    deprecated
    register_aliases
    params_aliased

    Notes
    -----
    This is a decorator that can be applied to a function or method to indicate its
    alias status. To make it effective, :func:`register_aliases` must be applied to the
    parent module or class.

    """

    def decorator(func):
        msg = _renaming_note(name, ver, warn)
        func.__doc__ = _note_into_doc(msg, func.__doc__)

        @wraps(func)
        def wrapper(*args, **kwargs):
            if warn:
                _warn_renamed(func, name, ver)
            return func(*args, **kwargs)

        wrapper._skipdoc = True
        func._alias = (name, wrapper)
        return func

    return decorator


def register_aliases(cls):
    """Register aliases of members of a module or class.

    See Also
    --------
    aliased

    Notes
    -----
    This function may serve as a decorator applied to a class to register aliases for
    its member methods:

    .. code-block:: python

       @register_aliases
       class MyClass():

    Alternatively, this function may be called at the end of a module to register
    aliases for its member functions.

    .. code-block:: python

       register_aliases(sys.modules[__name__])

    """
    aliases = []
    for func in cls.__dict__.values():
        if hasattr(func, "_alias"):
            aliases.append(func._alias)
    for name, wrapper in aliases:
        setattr(cls, name, wrapper)
    return cls


def params_aliased(params=[]):
    r"""Create aliases for parameters of a function or method.

    Parameters
    ----------
    params : list of tuple
        Aliases of parameters. Each tuple has four elements: param, alias, ver, and
        warn. Refer to :func:`aliased` for definitions.

    See Also
    --------
    aliased

    Notes
    -----
    This is a decorator that can be applied to a function or method to create aliases
    for its parameters. Unlike :func:`aliased`, this decorator does not require
    :func:`register_aliases` added to the module or class.

    """

    def decorator(func):
        for param, alias, ver, warn in params:
            msg = _renaming_note(alias, ver, warn)
            if doc := func.__doc__:
                func.__doc__ = _note_into_doc_param(msg, doc, param)

        @wraps(func)
        def wrapper(*args, **kwargs):
            for param, alias, ver, warn in params:
                if alias in kwargs:
                    kwargs[param] = kwargs.pop(alias)
                    if warn:
                        _warn_param_renamed(func, param, alias, ver)
            return func(*args, **kwargs)

        return wrapper

    return decorator
