# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Randomization."""

import numpy as np


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
    See NEP 19 [2]_ for an introduction to this change.

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

    Meanwhile, scikit-bio respects the legacy random generator (``RandomState``) [3]_.
    If ``np.random.seed`` has been called, or a ``RandomState`` instance is provided,
    scikit-bio will create a new random generator from a seed selected by the legacy
    random generator. This ensures reproducibility of legacy code and compatibility
    with packages that use the legacy mechanism. For example:

    .. code-block:: python

       np.random.seed(42)
       skbio_func1(...)
       skbio_func2(...)
       ...

    It should be noted that the legacy random generator will not be directly used by
    scikit-bio functions for random number generation. Only the new random generator
    will be exposed to the functions.

    References
    ----------
    .. [1] https://numpy.org/devdocs/reference/random/generator.html

    .. [2] https://numpy.org/neps/nep-0019-rng-policy.html

    .. [3] https://numpy.org/doc/stable/reference/random/legacy.html

    """
    # seed is new Generator: directly return it
    try:
        if isinstance(seed, np.random.Generator):
            return seed
    except AttributeError:
        raise ValueError(
            "The installed NumPy version does not support random.Generator. "
            "Please use NumPy >= 1.17."
        )

    # no seed: create a generator from a seed selected by the legacy random state
    # this ensures reproducibility of legacy code starting with `np.random.seed(x)`
    if seed is None:
        # the result of `np.iinfo(np.int32).max + 1` is usually 2,147,483,648, or
        # 1 << 31, # which is also the result of `np.random.get_state()[1][0]` if
        # `np.random` hasn't been called
        seed = np.random.randint(np.iinfo(np.int32).max + 1)

    # seed is legacy RandomState: create a generator from a seed selected by it
    elif isinstance(seed, np.random.RandomState):
        seed = seed.randint(np.iinfo(np.int32).max + 1)

    # seed is integer: create a generator from this random state
    elif not np.issubdtype(type(seed), np.integer):
        raise ValueError(
            "Invalid seed. It must be an integer or a random generator instance."
        )

    return np.random.default_rng(seed)
