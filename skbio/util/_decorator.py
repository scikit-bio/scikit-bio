# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
from functools import wraps
from collections import namedtuple

from ._exception import OverrideError
from ._warning import _warn_renamed


def _alias_msg(name, since=None, until=None, warn=False):
    """Create a message indicating the alias status of a function or a parameter."""
    if not since:
        msg = f"Alias: ``{name}``"
    else:
        msg = (
            f".. versionchanged:: {since} "
            f"Renamed from ``{name}``. The old name is kept as an alias"
        )
        if warn:
            msg += " but is deprecated"
            if until:
                msg += f" and will be removed in {until}"
        msg += "."
    return msg


def _msg_into_doc(msg, doc):
    """Insert a message into a docstring under the description.

    Parameters
    ----------
    msg : str
        Message to insert.
    doc : str
        Docstring to modify.

    Returns
    -------
    str
        Modified docstring.

    """
    # no docstring: message becomes the entire docstring.
    if doc is None:
        return msg + "\n"

    lines = doc.splitlines()
    n = len(lines)

    # find indentation size, which is important for the docstring to be rendered
    indent = 0
    for line in lines[1:]:
        if line:
            indent = len(line) - len(line.lstrip())
            break

    indent = " " * indent

    # docstring has at least two paragraphs: insert message after the first paragraph
    for i in range(1, n):
        if not lines[i]:
            return "\n".join(lines[:i] + ["", indent + msg] + lines[i:])

    # docstring has only one paragraph: append message to the end of docstring
    return doc + "\n" + indent + msg + "\n"


def aliased(name, since=None, until=None, warn=False):
    """Create an alias for a function or method.

    Parameters
    ----------
    name : str
        Alias name of the function or method.
    since : str, optional
        Version when alias was created.
    until : str, optional
        Version when alias will be removed.
    warn : bool, optional
        Raise a deprecation warning when alias is called.

    See Also
    --------
    add_aliases

    Notes
    -----
    This is a decorator that can be applied to a function or method to indicate its
    alias status.

    """

    def decorator(func):
        func._alias = (name, since, until, warn)
        return func

    return decorator


def add_aliases(cls):
    """Register aliases of members of a module or class.

    For each member that has an alias, create a wrapper named as the alias, and add a
    "skip" flag which instructs Sphinx autodoc to skip it when rendering the
    documentation.

    See Also
    --------
    aliased

    """
    toadd = []
    for func in cls.__dict__.values():
        if hasattr(func, "_alias"):
            toadd.append(func)

    for func in toadd:
        name, since, until, warn = func._alias

        @wraps(func)
        def wrapper(*args, **kwargs):
            if warn:
                _warn_renamed(func, name, since, until)
            return func(*args, **kwargs)

        wrapper._skip = True
        delattr(wrapper, "_alias")
        setattr(cls, name, wrapper)

        msg = _alias_msg(name, since, until, warn)
        func.__doc__ = _msg_into_doc(msg, func.__doc__)

    return cls


# parameter alias
ParamAlias = namedtuple(
    "ParamAlias",
    ["param", "alias", "since", "until", "warn"],
    defaults=[None, None, False],
)


def _param_msg_into_doc(param, msg, doc):
    """Insert a message into a docstring under the description of a parameter.

    Parameters
    ----------
    param : str
        Target parameter.
    msg : str
        Message to insert.
    doc : str
        Docstring to modify.

    Returns
    -------
    str
        Modified docstring.

    Raises
    ------
    ValueError
        Parameter is missing or its format is invalid.

    """
    # Find the header line of the parameter.
    match = re.search(rf"(^\s*{param}\s*:.*?\n)", doc, re.MULTILINE)
    if not match:
        raise ValueError(
            f'Parameter "{param}" is missing from the docstring or its format is '
            "invalid."
        )
    header = match.group(1)

    # Determine the indentation of the header line.
    indent = len(header) - len(header.lstrip())

    # Find the next line with the same or less indentation, which indicates the next
    # parameter or the end of the "Parameters" section.
    end = match.end()
    after = doc[end:]
    match = re.search(rf"(^\s{{0,{indent}}}[^\s].*?$)", after, re.MULTILINE)

    # Determine the insertion point.
    insert = end + (match.start() if match else len(after))

    # Insert the message. It should have 1+ indentation level (i.e., 4 spaces) than
    # the header line.
    return doc[:insert] + "\n" + " " * (indent + 4) + msg + "\n\n" + doc[insert:]


def params_aliased(params=[]):
    r"""Create aliases for parameters of a function or method.

    Parameters
    ----------
    params : list of ParamAlias
        Aliases of parameters.

    See Also
    --------
    aliased

    Notes
    -----
    This is a decorator that can be applied to a function or method to create an alias
    for one of its parameters. It can be applied multiple times to create aliases for
    multiple parameters. Unlike :func:`aliased`, this decorator does not require
    :func:`add_aliases` added to the module or class.

    """

    def decorator(func):
        for param, alias, since, until, warn in params:
            msg = _alias_msg(alias, since, until, warn)
            func.__doc__ = _param_msg_into_doc(param, msg, func.__doc__)

        @wraps(func)
        def wrapper(*args, **kwargs):
            for param, alias, since, until, warn in params:
                if alias in kwargs:
                    kwargs[param] = kwargs.pop(alias)
                    if warn:
                        _warn_renamed(func, alias, since, until, param=param)
            return func(*args, **kwargs)

        return wrapper

    return decorator


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
