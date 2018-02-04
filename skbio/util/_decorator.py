# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import warnings
import textwrap

import decorator

from ._exception import OverrideError
from ._warning import DeprecationWarning as SkbioDeprecationWarning


class _state_decorator:
    """ Base class for decorators of all public functionality.
    """

    _required_kwargs = ()

    def _get_indentation_level(self, docstring_lines,
                               default_existing_docstring=4,
                               default_no_existing_docstring=0):
        """ Determine the level of indentation of the docstring to match it.

            The indented content after the first line of a docstring can
            differ based on the nesting of the functionality being documented.
            For example, a top-level function may have its "Parameters" section
            indented four-spaces, but a method nested under a class may have
            its "Parameters" section indented eight spaces. This function
            determines the indentation level of the first non-whitespace line
            following the initial summary line.
        """
        # if there is no existing docstring, return the corresponding default
        if len(docstring_lines) == 0:
            return default_no_existing_docstring

        # if there is an existing docstring with only a single line, return
        # the corresponding default
        if len(docstring_lines) == 1:
            return default_existing_docstring

        # find the first non-blank line (after the initial summary line) and
        # return the number of leading spaces on that line
        for line in docstring_lines[1:]:
            if len(line.strip()) == 0:
                # ignore blank lines
                continue
            else:
                return len(line) - len(line.lstrip())

        # if there is an existing docstring with only a single non-whitespace
        # line, return the corresponding default
        return default_existing_docstring

    def _update_docstring(self, docstring, state_desc,
                          state_desc_prefix='State: '):
        # Hande the case of no initial docstring
        if docstring is None:
            return "%s%s" % (state_desc_prefix, state_desc)

        docstring_lines = docstring.split('\n')
        docstring_content_indentation = \
            self._get_indentation_level(docstring_lines)

        # wrap lines at 79 characters, accounting for the length of
        # docstring_content_indentation and start_desc_prefix
        len_state_desc_prefix = len(state_desc_prefix)
        wrap_at = 79 - (docstring_content_indentation + len_state_desc_prefix)
        state_desc_lines = textwrap.wrap(state_desc, wrap_at)
        # The first line of the state description should start with
        # state_desc_prefix, while the others should start with spaces to align
        # the text in this section. This is for consistency with numpydoc
        # formatting of deprecation notices, which are done using the note
        # Sphinx directive.
        state_desc_lines[0] = '%s%s%s' % (' ' * docstring_content_indentation,
                                          state_desc_prefix,
                                          state_desc_lines[0])
        header_spaces = ' ' * (docstring_content_indentation +
                               len_state_desc_prefix)
        for i, line in enumerate(state_desc_lines[1:], 1):
            state_desc_lines[i] = '%s%s' % (header_spaces, line)

        new_doc_lines = '\n'.join(state_desc_lines)
        docstring_lines[0] = '%s\n\n%s' % (docstring_lines[0], new_doc_lines)
        return '\n'.join(docstring_lines)

    def _validate_kwargs(self, **kwargs):
        for required_kwarg in self._required_kwargs:
            if required_kwarg not in kwargs:
                raise ValueError('%s decorator requires parameter: %s' %
                                 (self.__class__, required_kwarg))


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

    _required_kwargs = ('as_of', )

    def __init__(self, *args, **kwargs):
        self._validate_kwargs(**kwargs)
        self.as_of = kwargs['as_of']

    def __call__(self, func):
        state_desc = 'Stable as of %s.' % self.as_of
        func.__doc__ = self._update_docstring(func.__doc__, state_desc)
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

    _required_kwargs = ('as_of', )

    def __init__(self, *args, **kwargs):
        self._validate_kwargs(**kwargs)
        self.as_of = kwargs['as_of']

    def __call__(self, func):
        state_desc = 'Experimental as of %s.' % self.as_of
        func.__doc__ = self._update_docstring(func.__doc__, state_desc)
        return func


class deprecated(_state_decorator):
    """ State decorator indicating deprecated functionality.

    Used to indicate that a public class or function is deprecated, meaning
    that its API will be removed in a future version of scikit-bio. Decorating
    functionality as deprecated will update its doc string to indicate the
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
    ...             reason='Use skbio.g().')
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
        .. note:: Deprecated as of 0.3.0 for removal in 0.3.3. Use skbio.g().
    <BLANKLINE>

    """

    _required_kwargs = ('as_of', 'until', 'reason')

    def __init__(self, *args, **kwargs):
        self._validate_kwargs(**kwargs)
        self.as_of = kwargs['as_of']
        self.until = kwargs['until']
        self.reason = kwargs['reason']

    def __call__(self, func, *args, **kwargs):
        state_desc = 'Deprecated as of %s for removal in %s. %s' %\
            (self.as_of, self.until, self.reason)
        func.__doc__ = self._update_docstring(func.__doc__, state_desc,
                                              state_desc_prefix='.. note:: ')

        def wrapped_f(*args, **kwargs):
            warnings.warn('%s is deprecated as of scikit-bio version %s, and '
                          'will be removed in version %s. %s' %
                          (func.__name__, self.as_of, self.until, self.reason),
                          SkbioDeprecationWarning)
            # args[0] is the function being wrapped when this is called
            # after wrapping with decorator.decorator, but why???
            return func(*args[1:], **kwargs)

        return decorator.decorator(wrapped_f, func)


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
            raise TypeError("Class-only method called on an instance. Use"
                            " '%s.%s' instead."
                            % (cls.__name__, self.__func__.__name__))
        return super().__get__(obj, cls)
