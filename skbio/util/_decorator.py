# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from warnings import warn
from textwrap import wrap

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


class _state_decorator(object):
    """ Base class for decorators of all public functionality.

    **This documentation will be ported to an easily accessible location before
    this code is merged. This location may be linked from the docstrings.**

    All public functionality in scikit-bio has a defined stability state.
    These states inform users and developers to what extent they can rely on
    different APIs in the package. Definitions of the stability states and the
    information associated with each follow.

    stable
    ------
    Functionality defined as stable is part of scikit-bio's backward-
    compatible API. Users can be confident that the API will not change without
    first passing through the deprecated state, typically for at least two
    release cycles. We make every effort to maintain the API of this code.

    The docstrings of stable functionality will indicate the first scikit-bio
    version where the functionality was considered stable.

    experimental
    ------------
    Functionality defined as experimental is being considered for addition to
    scikit-bio's stable API. Users are encouraged to use this code, but to be
    aware that its API may change or be removed. Experimental functionality
    will typically pass through the deprecated state before it is removed, but
    in rare cases it may be removed directly (for example, if a serious
    methodological flaw is discovered that makes the functionality
    scientifically invalid).

    The docstrings of experimental functionality will indicate the first
    scikit-bio version where the functionality was considered experimental.

    We aim to move functionality through the experimental phase quickly (for
    example, two releases before moving to stable), but we don't make specific
    promises about when experimental functionality will become stable. This
    aligns with our philosophy that we don't make promises about experimental
    APIs, only about stable APIs.

    deprecated
    ----------
    Functionality defined as deprecated is targeted for removal from
    scikit-bio. Users should transition away from using it.

    The docstrings of deprecated functionality will indicate the first version
    of scikit-bio where the functionality was deprecated, the version of
    scikit-bio when the functionality will be removed, and the reason for
    deprecation of the code (for example, because a function was determined to
    be scientifically invalid, or because the API was adapted, and users should
    be using a different version of the function).

    Using deprecated functionality will raise a DeprecationWarning.

    """

    _line_prefix = '\n        '

    def _update_doc_string(self, func, state_desc):
        doc_lines = func.__doc__.split('\n')
        # wrap lines at 79 characters, accounting for the length of
        # self._line_prefix
        state_desc_lines = wrap('State: ' + state_desc,
                                79 - len(self._line_prefix))
        doc_lines.insert(
            1, self._line_prefix + self._line_prefix.join(state_desc_lines))
        return '\n'.join(doc_lines)


class stable(_state_decorator):
    """ State decorator indicating stable functionality.

    Used to indicate that public functionality is considered ``stable``,
    meaning that its API will be backward compatible unless it is deprecated.
    Decorating functionality as stable will update its doc string to indicate
    the first version of scikit-bio when the functionality was considered
    stable.

    Parameters
    ----------
    as_of : str
        First release version where functionality is considered to be stable.

    See Also
    --------
    experimental
    deprecated

    Examples
    --------
    >>> @stable(as_of='0.3.0')
    ... def f_stable():
    ...     \"\"\" An example stable function.
    ...     \"\"\"
    ...     pass
    >>> help(f_stable)
    Help on function f_stable in module skbio.util._decorator:
    <BLANKLINE>
    f_stable()
        An example stable function.
    <BLANKLINE>
        State: Stable as of 0.3.0.
    <BLANKLINE>
    """

    def __init__(self, **kwargs):
        self.as_of = kwargs['as_of']

    def __call__(self, func):
        state_desc = 'Stable as of %s.' % self.as_of
        func.__doc__ = self._update_doc_string(func, state_desc)
        return func


class experimental(_state_decorator):
    """ State decorator indicating experimental functionality.

    Used to indicate that public functionality is considered experimental,
    meaning that its API is subject to change or removal with little or
    (rarely) no warning. Decorating functionality as experimental will update
    its doc string to indicate the first version of scikit-bio when the
    functionality was considered experimental.

    Parameters
    ----------
    as_of : str
        First release version where feature is considered to be experimental.

    See Also
    --------
    stable
    deprecated

    Examples
    --------
    >>> @experimental(as_of='0.3.0')
    ... def f_experimental():
    ...     \"\"\" An example experimental function.
    ...     \"\"\"
    ...     pass
    >>> help(f_experimental)
    Help on function f_experimental in module skbio.util._decorator:
    <BLANKLINE>
    f_experimental()
        An example experimental function.
    <BLANKLINE>
        State: Experimental as of 0.3.0.
    <BLANKLINE>

    """

    def __init__(self, **kwargs):
        self.as_of = kwargs['as_of']

    def __call__(self, func):
        state_desc = 'Experimental as of %s.' % self.as_of
        func.__doc__ = self._update_doc_string(func, state_desc)
        return func


class deprecated(_state_decorator):
    """ State decorator indicating deprecated functionality.

    Used to indicate that a public class or function is deprecated, meaning
    that its API will be removed in a future version of scikit-bio. Decorating
    functionality as experimental will update its doc string to indicate the
    first version of scikit-bio when the functionality was deprecated, the
    first version of scikit-bio when the functionality will no longer exist,
    and the reason for deprecation of the API. It will also cause calls to the
    API to raise a ``DeprecationWarning``.

    Parameters
    ----------
    as_of : str
        First development version where feature is considered to be deprecated.
    until : str
        First release version where feature will no longer exist.
    reason : str
        Brief description of why the API is deprecated.

    See Also
    --------
    stable
    experimental

    Examples
    --------
    >>> @deprecated(as_of='0.3.0', until='0.3.3',
    ...             reason='Users should now use skbio.g().')
    ... def f_deprecated(x, verbose=False):
    ...     \"\"\" An example deprecated function.
    ...     \"\"\"
    ...     pass
    >>> help(f_deprecated)
    Help on function f_deprecated in module skbio.util._decorator:
    <BLANKLINE>
    f_deprecated(x, verbose=False)
        An example deprecated function.
    <BLANKLINE>
        State: Deprecated as of 0.3.0 for removal in 0.3.3. Users should now
        use skbio.g().
    <BLANKLINE>

    """

    def __init__(self, **kwargs):
        self.as_of = kwargs['as_of']
        self.until = kwargs['until']
        self.reason = kwargs['reason']

    def __call__(self, func, *args, **kwargs):
        state_desc = 'Deprecated as of %s for removal in %s. %s' %\
         (self.as_of, self.until, self.reason)
        func.__doc__ = self._update_doc_string(func, state_desc)

        def wrapped_f(*args, **kwargs):
            warn('%s is deprecated as of scikit-bio version %s, and will be'
                 ' removed in version %s. %s' %
                 (func.__name__, self.as_of, self.until, self.reason),
                 DeprecationWarning)
            return func(*args, **kwargs)
        # Currently func's signature is lost (and this is caught by the
        # doctests). It looks like the solution is probably to use the
        # [decorator module](https://pypi.python.org/pypi/decorator) which
        # @ebolyen pointed me to. I'm out of time for experimentation now
        # though. This needs to be resolved before this code is merged.
        # For testing purposes, we can test the signature of a deprecated
        # function ``f`` with:
        # import inspect
        # inspect.getargspec(f)
        wrapped_f.__name__ = func.__name__
        wrapped_f.__doc__ = func.__doc__

        return wrapped_f
