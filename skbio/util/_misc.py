# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import hashlib
import inspect
from types import FunctionType
from numbers import Integral

import numpy as np


def resolve_key(obj, key):
    """Resolve key given an object and key."""
    if callable(key):
        return key(obj)
    elif hasattr(obj, "metadata"):
        return obj.metadata[key]
    raise TypeError(
        "Could not resolve key %r. Key must be callable or %s must"
        " have `metadata` attribute." % (key, obj.__class__.__name__)
    )


def make_sentinel(name):
    return type(
        name,
        (),
        {"__repr__": lambda s: name, "__str__": lambda s: name, "__class__": None},
    )()


def find_sentinels(function, sentinel):
    params = inspect.signature(function).parameters
    return [name for name, param in params.items() if param.default is sentinel]


class MiniRegistry(dict):
    def __call__(self, name):
        """Act as a decorator to register functions with self."""

        def decorator(func):
            self[name] = func
            return func

        return decorator

    def copy(self):
        """Use for inheritance."""
        return self.__class__(super(MiniRegistry, self).copy())

    def formatted_listing(self):
        """Produce an RST list with descriptions."""
        if len(self) == 0:
            return "\tNone"
        else:
            return "\n".join(
                [
                    "\t%r\n\t  %s" % (name, self[name].__doc__.split("\n")[0])
                    for name in sorted(self)
                ]
            )

    def interpolate(self, obj, name):
        """Inject the formatted listing in the second blank line of `name`."""
        f = getattr(obj, name)
        f2 = FunctionType(
            f.__code__,
            f.__globals__,
            name=f.__name__,
            argdefs=f.__defaults__,
            closure=f.__closure__,
        )

        # Conveniently the original docstring is on f2, not the new ones if
        # inheritance is happening. I have no idea why.
        t = f2.__doc__.split("\n\n")
        t.insert(2, self.formatted_listing())
        f2.__doc__ = "\n\n".join(t)

        setattr(obj, name, f2)


def chunk_str(s, n, char):
    """Insert `char` character every `n` characters in string `s`.

    Canonically pronounced "chunkster".

    """
    # Modified from http://stackoverflow.com/a/312464/3776794
    if n < 1:
        raise ValueError(
            "Cannot split string into chunks with n=%d. n must be >= 1." % n
        )
    return char.join((s[i : i + n] for i in range(0, len(s), n)))


def cardinal_to_ordinal(n):
    """Return ordinal string version of cardinal int `n`.

    Parameters
    ----------
    n : int
        Cardinal to convert to ordinal. Must be >= 0.

    Returns
    -------
    str
        Ordinal version of cardinal `n`.

    Raises
    ------
    ValueError
        If `n` is less than 0.

    Notes
    -----
    This function can be useful when writing human-readable error messages.

    Examples
    --------
    >>> from skbio.util import cardinal_to_ordinal
    >>> cardinal_to_ordinal(0)
    '0th'
    >>> cardinal_to_ordinal(1)
    '1st'
    >>> cardinal_to_ordinal(2)
    '2nd'
    >>> cardinal_to_ordinal(3)
    '3rd'

    """
    # Taken and modified from http://stackoverflow.com/a/20007730/3776794
    # Originally from http://codegolf.stackexchange.com/a/4712 by Gareth
    if n < 0:
        raise ValueError("Cannot convert negative integer %d to ordinal " "string." % n)
    return "%d%s" % (n, "tsnrhtdd"[(n // 10 % 10 != 1) * (n % 10 < 4) * n % 10 :: 4])


def safe_md5(open_file, block_size=2**20):
    """Compute an md5 sum without loading the file into memory.

    Parameters
    ----------
    open_file : file object
        open file handle to the archive to compute the checksum. It
        must be open as a binary file
    block_size : int, optional
        size of the block taken per iteration

    Returns
    -------
    md5 : md5 object from the hashlib module
        object with the loaded file

    Notes
    -----
    This method is based on the answers given in:
    http://stackoverflow.com/a/1131255/379593

    Examples
    --------
    >>> from io import BytesIO
    >>> from skbio.util import safe_md5
    >>> fd = BytesIO(b"foo bar baz") # open file like object
    >>> x = safe_md5(fd)
    >>> x.hexdigest()
    'ab07acbb1e496801937adfa772424bf7'
    >>> fd.close()

    """
    md5 = hashlib.md5()
    data = True
    while data:
        data = open_file.read(block_size)
        if data:
            md5.update(data)
    return md5


def find_duplicates(iterable):
    """Find duplicate elements in an iterable.

    Parameters
    ----------
    iterable : iterable
        Iterable to be searched for duplicates (i.e., elements that are
        repeated).

    Returns
    -------
    set
        Repeated elements in `iterable`.

    """
    # modified from qiita.qiita_db.util.find_repeated
    # https://github.com/biocore/qiita
    # see licenses/qiita.txt
    seen, repeated = set(), set()
    for e in iterable:
        if e in seen:
            repeated.add(e)
        else:
            seen.add(e)
    return repeated


def get_rng(seed=None):
    r"""Get a random generator.

    .. versionchanged:: 0.6.3
        Added legacy random generator support.

    Parameters
    ----------
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance.

    Returns
    -------
    np.random.Generator
        Random generator instance.

    See Also
    --------
    numpy.random.Generator
    numpy.random.RandomState

    Notes
    -----
    A random generator ensures reproducibility of outputs. scikit-bio utilizes NumPy's
    new random generator (``Generator``) [1]_, which was introduced in version 1.17.
    See NEP 19 [3]_ for an introduction to this change.

    The following code demonstrates the recommended usage of the random generator with
    various scikit-bio functions that are stochastic. With a random generator created
    in advance, you can plug it into multiple function calls. The results of the entire
    code will be reproducible.

    (42 is an arbitrarily chosen random seed. It can be any non-negative integer.)

    .. code-block:: python

       rng = np.random.default_rng(42)
       skbio_func1(..., seed=rng)
       skbio_func2(..., seed=rng)
       ...

    Alternatively, you may specify an integer seed to make the result of each function
    call reproducible, such as:

    .. code-block:: python

       skbio_func1(..., seed=42)
       skbio_func2(..., seed=42)
       ...

    Meanwhile, scikit-bio respects the legacy random generator (``RandomState``) [4]_.
    If ``np.random.seed`` has been called, or a ``RandomState`` instance is provided.
    This ensures reproducibility of legacy code and compatibility with packages that
    use the legacy mechanism. For example:

    .. code-block:: python

       np.random.seed(42)
       skbio_func1(...)
       skbio_func2(...)
       ...

    References
    ----------
    .. [1] https://numpy.org/devdocs/reference/random/generator.html

    .. [2] https://numpy.org/doc/stable/reference/random/legacy.html

    .. [3] https://numpy.org/neps/nep-0019-rng-policy.html

    .. [4] https://numpy.org/doc/stable/reference/random/legacy.html

    """
    # Seed is a new Generator: directly return it.
    # This code will work for legacy NumPy (<1.17) without Generator.
    has_rng = hasattr(np.random, "Generator")
    if has_rng and isinstance(seed, np.random.Generator):
        return seed

    # Seed is a legacy RandomState: directly return it.
    # This code may work for future NumPy without RandomState.
    has_rs = hasattr(np.random, "RandomState")
    if has_rs and isinstance(seed, np.random.RandomState):
        return seed

    # Seed is an integer: create a generator from it.
    # Integral is compatible with NumPy integers.
    if isinstance(seed, Integral):
        # Return a new generator if available, otherwise a legacy one.
        if has_rng:
            return np.random.default_rng(seed)
        if has_rs:
            return np.random.RandomState(seed)
        raise ValueError("Supported random generators are not available.")

    # Seed is not provided: create a generator.
    if seed is None:
        # Recover the global random generator, which should be a legacy RandomState
        # instance and could be set by np.random.seed.
        # This behavior is consistent with SciPy and scikit-learn's check_random_state.
        if has_rs:
            return np.random.mtrand._rand

        # However, a potential method to adopt the new generator stochasticity is:
        #   seed = np.random.randint(np.iinfo(np.int32).max + 1)
        #   return np.random.default_rng(seed)
        # Here, the result of np.iinfo(np.int32).max + 1 is usually 2,147,483,648, or
        # 1 << 31, which is also the result of np.random.get_state()[1][0] if np.random
        # hasn't been called.

        # Otherwise (for future NumPy), return a new generator without seed, which is
        # stochastic.
        if has_rng:
            return np.random.default_rng()
        raise ValueError("Supported random generators are not available.")

    # Seed is invalid.
    raise ValueError(
        "Invalid seed. It must be an integer or a random generator instance."
    )
