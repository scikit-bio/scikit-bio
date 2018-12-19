Coding guidelines
=================

As project size increases, consistency of the code base and documentation becomes more important. We therefore provide guidelines for code and documentation that is contributed to scikit-bio. Our goal is to create a consistent code base where:

* it is easy to find relevant functionality (and to determine when functionality that you're looking for doesn't exist),
* you can trust that the code that you're working with is sufficiently tested, and
* names and interfaces are intuitive.

**As scikit-bio is in beta, our coding guidelines are presented here as a working draft. These guidelines are requirements for all code submitted to scikit-bio, but at this stage the guidelines themselves are malleable. If you disagree with something, or have a suggestion for something new to include, you should** `create an issue`_ **to initiate a discussion.**

.. _`create an issue`: https://github.com/biocore/scikit-bio/issues

What are the naming conventions? and How should I format my code?
-----------------------------------------------------------------

We adhere to the `PEP 8`_ python coding guidelines for code and documentation standards. Before submitting any code to scikit-bio, you should read these carefully and apply the guidelines in your code.

.. _`PEP 8`: http://legacy.python.org/dev/peps/pep-0008/


What should I call my variables?
--------------------------------

- *Choose the name that people will most likely guess.* Make it descriptive, but not too long: ``curr_record`` is better than ``c``, or ``curr``, or ``current_genbank_record_from_database``.

- *Good names are hard to find.* Don't be afraid to change names except when they are part of interfaces that other people are also using. It may take some time working with the code to come up with reasonable names for everything: if you have unit tests, it's easy to change them, especially with global search and replace.

- *Use singular names for individual things, plural names for collections.* For example, you'd expect ``self.name`` to hold something like a single string, but ``self.names`` to hold something that you could loop through like a list or dictionary. Sometimes the decision can be tricky: is ``self.index`` an integer holding a positon, or a dictionary holding records keyed by name for easy lookup? If you find yourself wondering these things, the name should probably be changed to avoid the problem: try ``self.position`` or ``self.look_up``.

- *Don't make the type part of the name.* You might want to change the implementation later. Use ``Records`` rather than ``RecordDict`` or ``RecordList``, etc. Don't use Hungarian Notation either (i.e. where you prefix the name with the type).

- *Make the name as precise as possible.* If the variable is the path of the input file, call it ``input_fp``, not ``input`` or ``file`` (which you shouldn't use anyway, since they're keywords), and not ``infile`` (because that looks like it should be a file object, not just its name).

- *Use* ``result`` *to store the value that will be returned from a method or function.* Use ``data`` for input in cases where the function or method acts on arbitrary data (e.g. sequence data, or a list of numbers, etc.) unless a more descriptive name is appropriate.

- *One-letter variable names should only occur in math functions or as loop iterators with limited scope.* Limited scope covers things like ``for k in keys: print k``, where ``k`` survives only a line or two. Loop iterators should refer to the variable that they're looping through: ``for k in keys, i in items``, or ``for key in keys, item in items``. If the loop is long or there are several 1-letter variables active in the same scope, rename them.

- *Limit your use of abbreviations.* A few well-known abbreviations are OK, but you don't want to come back to your code in 6 months and have to figure out what ``sptxck2`` is. It's worth it to spend the extra time typing ``species_taxon_check_2``, but that's still a horrible name: what's check number 1? Far better to go with something like ``taxon_is_species_rank`` that needs no explanation, especially if the variable is only used once or twice.

Acceptable abbreviations
^^^^^^^^^^^^^^^^^^^^^^^^

The following list of abbreviations can be considered well-known and used with impunity within mixed name variables, but some should not be used by themselves as they would conflict with common functions, python built-in's, or raise an exception. Do not use the following by themselves as variable names: ``dir``,  ``exp`` (a common ``math`` module function), ``in``, ``max``, and ``min``. They can, however, be used as part of a name, eg ``matrix_exp``.

+--------------------+--------------+
|        Full        |  Abbreviated |
+====================+==============+
|          alignment |          aln |
+--------------------+--------------+
|           archaeal |         arch |
+--------------------+--------------+
|          auxiliary |          aux |
+--------------------+--------------+
|          bacterial |         bact |
+--------------------+--------------+
|           citation |         cite |
+--------------------+--------------+
|            current |         curr |
+--------------------+--------------+
|           database |           db |
+--------------------+--------------+
|         dictionary |         dict |
+--------------------+--------------+
|          directory |          dir |
+--------------------+--------------+
|    distance matrix |           dm |
+--------------------+--------------+
|        end of file |          eof |
+--------------------+--------------+
|         eukaryotic |          euk |
+--------------------+--------------+
|          filepath  |           fp |
+--------------------+--------------+
|          frequency |         freq |
+--------------------+--------------+
|           expected |          exp |
+--------------------+--------------+
|              index |          idx |
+--------------------+--------------+
|              input |           in |
+--------------------+--------------+
|            maximum |          max |
+--------------------+--------------+
|            minimum |          min |
+--------------------+--------------+
|      mitochondrial |           mt |
+--------------------+--------------+
|             number |          num |
+--------------------+--------------+
|        observation |          obs |
+--------------------+--------------+
|           observed |          obs |
+--------------------+--------------+
|           original |         orig |
+--------------------+--------------+
|             output |          out |
+--------------------+--------------+
|          parameter |        param |
+--------------------+--------------+
|          phylogeny |        phylo |
+--------------------+--------------+
|           previous |         prev |
+--------------------+--------------+
|        probability |         prob |
+--------------------+--------------+
|            protein |         prot |
+--------------------+--------------+
|             record |          rec |
+--------------------+--------------+
|          reference |          ref |
+--------------------+--------------+
|           sequence |          seq |
+--------------------+--------------+
| standard deviation |        stdev |
+--------------------+--------------+
|         statistics |        stats |
+--------------------+--------------+
|             string |          str |
+--------------------+--------------+
|          structure |       struct |
+--------------------+--------------+
|          temporary |         temp |
+--------------------+--------------+
|               taxa |          tax |
+--------------------+--------------+
|              taxon |          tax |
+--------------------+--------------+
|          taxonomic |          tax |
+--------------------+--------------+
|           taxonomy |          tax |
+--------------------+--------------+
|           variance |          var |
+--------------------+--------------+

How do I organize my modules (source files)?
--------------------------------------------

- *Have a docstring with a description of the module's functions*. If the description is long, the first line should be a short summary that makes sense on its own, separated from the rest by a newline.

- *All code, including import statements, should follow the docstring.* Otherwise, the docstring will not be recognized by the interpreter, and you will not have access to it in interactive sessions (i.e. through ``obj.__doc__``) or when generating documentation with automated tools.

- *Import built-in modules first, followed by third-party modules, followed by any changes to the path and your own modules.* Especially, additions to the path and names of your modules are likely to change rapidly: keeping them in one place makes them easier to find.

- *Don't use* ``from module import *``, *instead use* ``from module import Name, Name2, Name3...`` *or possibly* ``import module``. This makes it *much* easier to see name collisions and to replace implementations.

- If you are importing `NumPy`_, `Matplotlib`_, or another package that encourages a standard style for their import statements use them as needed for example:

::

    import numpy as np
    import numpy.testing as npt
    import pandas as pd

    from matplotlib import pyplot as plt

.. _`NumPy`: http://www.numpy.org/
.. _`Matplotlib`: http://matplotlib.org/

Example of module structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The structure of your module should be similar to the example below. scikit-bio uses the `NumPy doc`_ standard for documentation. Our `doc/README.md`_ explains how to write your docstrings using the `NumPy doc`_ standards for scikit-bio:

.. _`NumPy doc`: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
.. _`doc/README.md`: https://github.com/biocore/scikit-bio/blob/master/doc/README.md

.. code-block:: python

    r"""
    Numbers (:mod:`skbio.numbers`)
    ==============================

    .. currentmodule:: skbio.numbers

    Numbers holds a sequence of numbers, and defines several statistical
    operations (mean, stdev, etc.) FrequencyDistribution holds a mapping from
    items (not necessarily numbers) to counts, and defines operations such as
    Shannon entropy and frequency normalization.


    Classes
    -------

    .. autosummary::
       :toctree: generated/

       Numbers

    """

    # ----------------------------------------------------------------------------
    # Copyright (c) 2013--, scikit-bio development team.
    #
    # Distributed under the terms of the Modified BSD License.
    #
    # The full license is in the file COPYING.txt, distributed with this software.
    # ----------------------------------------------------------------------------

    from random import choice, random

    import numpy as np
    from utils import indices


    class Numbers(list):
        pass


    class FrequencyDistribution(dict):
        pass


How should I write comments?
----------------------------

- *Always update the comments when the code changes.* Incorrect comments are far worse than no comments, since they are actively misleading.

- *Comments should say more than the code itself.* Examine your comments carefully: they may indicate that you'd be better off rewriting your code (especially if *renaming your variables* would allow you to get rid of the comment.) In particular, don't scatter magic numbers and other constants that have to be explained through your code. It's far better to use variables whose names are self-documenting, especially if you use the same constant more than once. Also, think about making constants into class or instance data, since it's all too common for 'constants' to need to change or to be needed in several methods.

    +-------+------------------------------------------------------------+
    | Wrong |       ``win_size -= 20        # decrement win_size by 20`` |
    +-------+------------------------------------------------------------+
    |    OK | ``win_size -= 20        # leave space for the scroll bar`` |
    +-------+------------------------------------------------------------+
    | Right |                             ``self._scroll_bar_size = 20`` |
    +-------+------------------------------------------------------------+
    |       |                      ``win_size -= self._scroll_bar_size`` |
    +-------+------------------------------------------------------------+


- *Use comments starting with #, not strings, inside blocks of code.*
- *Start each method, class and function with a docstring using triple double quotes (""").* Make sure the docstring follows the `NumPy doc`_ standard.

- *Always update the docstring when the code changes.* Like outdated comments, outdated docstrings can waste a lot of time. "Correct examples are priceless, but incorrect examples are worse than worthless." `Jim Fulton`_.

.. _`Jim Fulton`: http://www.python.org/pycon/dc2004/papers/4/PyCon2004DocTestUnit.pdf

How should I test my code ?
---------------------------

There are several different approaches for testing code in python: ``nose``, ``unittest`` and ``numpy.testing``. Their purpose is the same, to check that execution of code given some input produces a specified output. The cases to which the approaches lend themselves are different.

Whatever approach is employed, the general principle is every line of code should be tested. It is critical that your code be fully tested before you draw conclusions from results it produces. For scientific work, bugs don't just mean unhappy users who you'll never actually meet: **they may mean retracted publications**.

Tests are an opportunity to invent the interface(s) you want. Write the test for a method before you write the method: often, this helps you figure out what you would want to call it and what parameters it should take. It's OK to write the tests a few methods at a time, and to change them as your ideas about the interface change. However, you shouldn't change them once you've told other people what the interface is. In the spirit of this, your tests should also import the functionality that they test from the shortest alias possible. This way any change to the API will cause your tests to break, and rightly so!

Never treat prototypes as production code. It's fine to write prototype code without tests to try things out, but when you've figured out the algorithm and interfaces you must rewrite it *with tests* to consider it finished. Often, this helps you decide what interfaces and functionality you actually need and what you can get rid of.

"Code a little test a little". For production code, write a couple of tests, then a couple of methods, then a couple more tests, then a couple more methods, then maybe change some of the names or generalize some of the functionality. If you have a huge amount of code where all you have to do is write the tests', you're probably closer to 30% done than 90%. Testing vastly reduces the time spent debugging, since whatever went wrong has to be in the code you wrote since the last test suite. And remember to use python's interactive interpreter for quick checks of syntax and ideas.

Run the test suite when you change `anything`. Even if a change seems trivial, it will only take a couple of seconds to run the tests and then you'll be sure. This can eliminate long and frustrating debugging sessions where the change turned out to have been made long ago, but didn't seem significant at the time. **Note that tests are executed using Travis CI**, see `this document's section`_ for further discussion.

.. _`this document's section`: https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md#testing-guidelines

Some pointers
^^^^^^^^^^^^^

- *Use the* ``unittest`` *or the* ``nose`` *framework with tests in a separate file for each module.* Name the test file ``test_module_name.py`` and include it inside the tests folder of the module. Keeping the tests separate from the code reduces the temptation to change the tests when the code doesn't work, and makes it easy to verify that a completely new implementation presents the same interface (behaves the same) as the old.

- *Always include an* ``__init__.py`` *file in your tests directory*. This is required for the module to be included when the package is built and installed via ``setup.py``.

- *Always import from a minimally deep API target*. That means you would use ``from skbio import DistanceMatrix`` instead of ``from skbio.stats.distance import DistanceMatrix``. This allows us prevent most cases of accidental regression in our API.

- *Use* ``numpy.testing`` *if you are doing anything with floating point numbers, arrays or permutations* (use ``numpy.testing.assert_almost_equal``). Do *not* try to compare floating point numbers using ``assertEqual`` if you value your sanity.

- *Test the interface of each class in your code by defining at least one* ``TestCase`` *with the name* ``ClassNameTests``. This should contain tests for everything in the public interface.

- *If the class is complicated, you may want to define additional tests with names* ``ClassNameTests_test_type``. These might subclass ``ClassNameTests`` in order to share ``setUp`` methods, etc.

- *Tests of private methods should be in a separate* ``TestCase`` *called* ``ClassNameTests_private``. Private methods may change if you change the implementation. It is not required that test cases for private methods pass when you change things (that's why they're private, after all), though it is often useful to have these tests for debugging.

- *Test `all` the methods in your class.* You should assume that any method you haven't tested has bugs. The convention for naming tests is ``test_method_name``. Any leading and trailing underscores on the method name can be ignored for the purposes of the test; however, *all tests must start with the literal substring* ``test`` *for* ``unittest`` and ``nose`` *to find them.* If the method is particularly complex, or has several discretely different cases you need to check, use ``test_method_name_suffix``, e.g. ``test_init_empty``, ``test_init_single``, ``test_init_wrong_type``, etc. for testing ``__init__``.

- *Docstrings for testing methods should be considered optional*, instead the description of what the method does should be included in the name itself, therefore the name should be descriptive enough such that when running the tests in verbose mode you can immediately see the file and test method that's failing.

.. code-block:: none

    $ python -c "import skbio; skbio.test(verbose=True)"
    skbio.maths.diversity.alpha.tests.test_ace.test_ace ... ok
    test_berger_parker_d (skbio.maths.diversity.alpha.tests.test_base.BaseTests) ... ok

    ----------------------------------------------------------------------
    Ran 2 tests in 0.1234s

    OK

- *Module-level functions should be tested in their own* ``TestCase``\ *, called* ``modulenameTests``. Even if these functions are simple, it's important to check that they work as advertised.

- *It is much more important to test several small cases that you can check by hand than a single large case that requires a calculator.* Don't trust spreadsheets for numerical calculations -- use R instead!

- *Make sure you test all the edge cases: what happens when the input is None, or '', or 0, or negative?* What happens at values that cause a conditional to go one way or the other? Does incorrect input raise the right exceptions? Can your code accept subclasses or superclasses of the types it expects? What happens with very large input?

- *To test permutations, check that the original and shuffled version are different, but that the sorted original and sorted shuffled version are the same.* Make sure that you get *different* permutations on repeated runs and when starting from different points.

- *To test random choices, figure out how many of each choice you expect in a large sample (say, 1000 or a million) using the binomial distribution or its normal approximation.* Run the test several times and check that you're within, say, 3 standard deviations of the mean.

- All tests that depend on a random value should be seeded, for example if using NumPy, `numpy.random.seed(0)` should be used, in any other case the appropriate API should be used to create consistent outputs between runs. It is preferable that you do this for each test case instead of doing it in the `setUp` function/method (if any exists).

- Stochastic failures should occur less than 1/10,1000 times, otherwise you risk adding a significant amount of time to the total running time of the test suite.

Example of a ``nose`` test module structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # ----------------------------------------------------------------------------
    # Copyright (c) 2013--, scikit-bio development team.
    #
    # Distributed under the terms of the Modified BSD License.
    #
    # The full license is in the file COPYING.txt, distributed with this software.
    # ----------------------------------------------------------------------------

    import numpy as np
    from nose.tools import assert_almost_equal, assert_raises

    from skbio.math.diversity.alpha.ace import ace


    def test_ace():
        assert_almost_equal(ace(np.array([2, 0])), 1.0)
        assert_almost_equal(ace(np.array([12, 0, 9])), 2.0)
        assert_almost_equal(ace(np.array([12, 2, 8])), 3.0)
        assert_almost_equal(ace(np.array([12, 2, 1])), 4.0)
        assert_almost_equal(ace(np.array([12, 1, 2, 1])), 7.0)
        assert_almost_equal(ace(np.array([12, 3, 2, 1])), 4.6)
        assert_almost_equal(ace(np.array([12, 3, 6, 1, 10])), 5.62749672)

        # Just returns the number of OTUs when all are abundant.
        assert_almost_equal(ace(np.array([12, 12, 13, 14])), 4.0)

        # Border case: only singletons and 10-tons, no abundant OTUs.
        assert_almost_equal(ace([0, 1, 1, 0, 0, 10, 10, 1, 0, 0]), 9.35681818182)


    def test_ace_only_rare_singletons():
        with assert_raises(ValueError):
            ace([0, 0, 43, 0, 1, 0, 1, 42, 1, 43])


    if __name__ == '__main__':
        import nose
        nose.runmodule()

Git pointers
------------

Commit messages are a useful way to document the changes being made to a project, it additionally documents who is making these changes and when are these changes being made, all of which are relevant when tracing back problems.

Authoring a commit message
^^^^^^^^^^^^^^^^^^^^^^^^^^

The most important metadata in a commit message is (arguably) the author's name and the author's e-mail. GitHub uses this information to attribute your contributions to a project, see for example the `scikit-bio list of contributors`_.

.. _`scikit-bio list of contributors`: https://github.com/biocore/scikit-bio/graphs/contributors

Follow `this guide`_ to set up your system and **make sure the e-mail you use in this step is the same e-mail associated to your GitHub account**.

.. _`this guide`: http://git-scm.com/book/en/Getting-Started-First-Time-Git-Setup

After doing this you should see your name and e-mail when you run the following commands:

.. code-block:: none

    $ git config --global user.name
    Yoshiki Vázquez Baeza
    $ git config --global user.email
    yoshiki89@gmail.com

Writing a commit message
^^^^^^^^^^^^^^^^^^^^^^^^

In general the writing of a commit message should adhere to `NumPy's guidelines`_ which if followed correctly will help you structure your changes better i. e. bug fixes will be in a commit followed by a commit updating the test suite and with one last commit that update the documentation as needed.

GitHub provides a set of handy features that will link together a commit message to a ticket in the issue tracker, this is specially helpful because you can `close an issue automatically`_ when the change is merged into the main repository, this reduces the amount of work that has to be done making sure outdated issues are not open.

.. _`NumPy's guidelines`: http://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html#writing-the-commit-message
.. _`close an issue automatically`: https://help.github.com/articles/closing-issues-via-commit-messages
