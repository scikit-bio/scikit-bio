# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from ._exception import OverrideError


# Adapted from http://stackoverflow.com/a/8313042/579416
def overrides(interface_class):
    """Decorator for class-level members.

    Used to indicate that a member is being overridden from a specific parent
    class. If the member does not have a docstring, it will pull one from the
    parent class. When chaining decorators, this should be first as it is
    relatively nondestructive.

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
    def overrider(method):
        if method.__name__ not in dir(interface_class):
            raise OverrideError("%r is not present in parent class: %r." %
                                (method.__name__, interface_class.__name__))
        if method.__doc__ is None:
            method.__doc__ = getattr(interface_class, method.__name__).__doc__
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
