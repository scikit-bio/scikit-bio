Adding a new module to skbio
############################

Each module needs an `__init__.py` file and a `tests` folder that also
contains an `__init__.py` file. For a module, a simple one may look
like this::

  r"""
  A module (:mod:`skbio.module`)
  ==============================

  .. currentmodule:: skbio.module

  Documentation for this module.
  """

  # ----------------------------------------------------------------------------
  # Copyright (c) 2013--, scikit-bio development team.
  #
  # Distributed under the terms of the Modified BSD License.
  #
  # The full license is in the file COPYING.txt, distributed with this software.
  # ----------------------------------------------------------------------------

  from skbio.util import TestRunner
  test = TestRunner(__file__).test

Usually, some functionality from the module will be made accessible by
importing it in `__init__.py`. It's convenient to use explicit
relative imports (`from .implementation import compute`), so that
functionality can be neatly separated in different files but the user
doesn't face a deeply nested package: `from skbio.module import
compute` instead of `from skbio.module.implementation import compute`.

Inside the tests folder, a simpler `__init__.py` works fine (it is
necessary so that all tests can be run after installation)::

  # ----------------------------------------------------------------------------
  # Copyright (c) 2013--, scikit-bio development team.
  #
  # Distributed under the terms of the Modified BSD License.
  #
  # The full license is in the file COPYING.txt, distributed with this software.
  # ----------------------------------------------------------------------------

Finally, remember to also follow the `documentation guidelines
<https://github.com/biocore/scikit-bio/blob/master/doc/README.md#documenting-a-module-in-scikit-bio>`_.
