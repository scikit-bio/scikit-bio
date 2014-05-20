Coding guidelines
=================

As project size increases, consistency increases in importance. Unit testing and a consistent style are critical to having trusted code to integrate. Also, guesses about names and interfaces will be correct more often.

What should I call my variables?
--------------------------------

- *Choose the name that people will most likely guess.* Make it descriptive, but not too long: ``curr_record`` is better than ``c``, or ``curr``, or ``current_genbank_record_from_database``.

- *Good names are hard to find.* Don't be afraid to change names except when they are part of interfaces that other people are also using. It may take some time working with the code to come up with reasonable names for everything: if you have unit tests, it's easy to change them, especially with global search and replace.

- *Use singular names for individual things, plural names for collections.* For example, you'd expect ``self.name`` to hold something like a single string, but ``self.names`` to hold something that you could loop through like a list or dict. Sometimes the decision can be tricky: is ``self.index`` an int holding a positon, or a dict holding records keyed by name for easy lookup? If you find yourself wondering these things, the name should probably be changed to avoid the problem: try ``self.position`` or ``self.look_up``.

- *Don't make the type part of the name.* You might want to change the implementation later. Use ``Records`` rather than ``RecordDict`` or ``RecordList``, etc. Don't use Hungarian Notation either (i.e. where you prefix the name with the type).

- *Make the name as precise as possible.* If the variable is the name of the input file, call it ``infile_name``, not ``input`` or ``file`` (which you shouldn't use anyway, since they're keywords), and not ``infile`` (because that looks like it should be a file object, not just its name).

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
|          auxillary |          aux |
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
|          taxonomic |          tax |
+--------------------+--------------+
|           variance |          var |
+--------------------+--------------+

What are the naming conventions?
--------------------------------

We adhere to the `PEP 8`_ python coding guidelines for code and documentation standards. Before submitting any code to scikit-bio, you should read these carefully and apply the guidelines in your code.

.. _`PEP 8`: http://legacy.python.org/dev/peps/pep-0008/

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

    from matplotlib import pyplot as plt

.. _`NumPy`: http://www.numpy.org/
.. _`Matplotlib`: http://matplotlib.org/

Example of module structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The structure of your module should be similar to the example bellow, note that scikit-bio uses the `NumPy doc`_ standard for documentation, `this document explains`_ how to do this:

.. _`NumPy doc`: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt
.. _`this document explains`: https://github.com/biocore/scikit-bio/blob/master/doc/README.md

.. code-block:: python

    r"""
    Number List (:mod:`skbio.core.numbers`)
    =======================================

    .. currentmodule:: skbio.core.numbers

    NumberList holds a sequence of numbers, and defines several statistical
    operations (mean, stdev, etc.) FrequencyDistribution holds a mapping from
    items (not necessarily numbers) to counts, and defines operations such as
    Shannon entropy and frequency normalization.


    Classes
    -------

    .. autosummary::
       :toctree: generated/

       NumberList

    """
    # ----------------------------------------------------------------------------
    # Copyright (c) 2013--, scikit-bio development team.
    #
    # Distributed under the terms of the Modified BSD License.
    #
    # The full license is in the file COPYING.txt, distributed with this software.
    # ----------------------------------------------------------------------------

    from __future__ import absolute_import, division, print_function

    import numpy as np
    from random import choice, random
    from utils import indices

    class NumberList(list):
        pass    # much code deleted
    class FrequencyDistribution(dict):
        pass    # much code deleted


How should I write comments?
----------------------------

- *Always update the comments when the code changes.* Incorrect comments are far worse than no comments, since they are actively misleading.

- *Comments should say more than the code itself.* Examine your comments carefully: they may indicate that you'd be better off rewriting your code (especially, *renaming your variables* and getting rid of the comment.) In particular, don't scatter magic numbers and other constants that have to be explained through your code. It's far better to use variables whose names are self-documenting, especially if you use the same constant more than once. Also, think about making constants into class or instance data, since it's all too common for 'constants' to need to change or to be needed in several methods.

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
- *Start each method, class and function with a docstring using triple double quotes (""").* The docstring should start with a 1-line description that makes sense by itself (many automated formatting tools, and the IDE, use this). This should be followed by a blank line, followed by descriptions of the parameters (if any). Finally, add any more detailed information, such as a longer description, notes about the algorithm, detailed notes about the parameters, etc. If there is a usage example, it should appear at the end. Make sure any descriptions of parameters have the correct spelling, case, etc. For example: ::

.. code-block:: python

    class BiologicalSequence(Sequence):
        """Base class for biological sequences.

        Parameters
        ----------
        sequence : python Sequence (e.g., str, list or tuple)
            The biological sequence.
        id : str, optional
            The sequence id (e.g., an accession number).
        description : str, optional
            A description or comment about the sequence (e.g., "green
            fluorescent protein").
        validate : bool, optional
            If True, runs the `is_valid` method after construction and raises
            BiologicalSequenceError if ``is_valid == False``.

        Attributes
        ----------
        description
        id

        Raises
        ------
        skbio.core.exception.BiologicalSequenceError
          If ``validate == True`` and ``is_valid == False``.

        See Also
        --------
        NucleotideSequence
        DNASequence
        RNASequence

        Notes
        -----
        `BiologicalSequence` objects are immutable. Where applicable, methods
        return a new object of the same class.
        Subclasses are typically defined by methods relevant to only a specific
        type of biological sequence, and by containing characters only contained in
        the IUPAC standard character set [1]_ for that molecule type.

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> t = BiologicalSequence('GGUCCUGAAGGU')

        References
        ----------
        .. [1] Nomenclature for incompletely specified bases in nucleic acid
           sequences: recommendations 1984.
           Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
           A Cornish-Bowden

        """

For more information refer to the `NumPy doc`_ standard.

- *Always update the docstring when the code changes.* Like outdated comments, outdated docstrings can waste a lot of time. "Correct examples are priceless, but incorrect examples are worse than worthless." `Jim Fulton`_.

How should I format my code?
----------------------------

- *Use 4 spaces for indentation.* Do not use tabs (set your editor to convert tabs to spaces). The behaviour of tabs is not predictable across platforms, and will cause syntax errors. If we all use the same indentation, collaboration is much easier.

- *Lines should not be longer than 79 characters.* Long lines are inconvenient in some editors. Use \\ for line continuation. Note that there cannot be whitespace after the \\.

- *Blank lines should be used to highlight class and method definitions.* Separate class definitions by two blank lines. Separate methods by one blank line.

- *Be consistent with the use of whitespace around operators.* Inconsistent whitespace makes it harder to see at a glance what is grouped together.

    +------+--------------------------+
    | Good |        ``((a+b)*(c+d))`` |
    +------+--------------------------+
    |   OK |  ``((a + b) * (c + d))`` |
    +------+--------------------------+
    |  Bad | ``( (a+ b)  *(c +d  ))`` |
    +------+--------------------------+

- *Don't put whitespace after delimiters or inside slicing delimiters.* Whitespace here makes it harder to see what's associated.

    +------+-------------+------------------+
    | Good |   ``(a+b)`` |         ``d[k]`` |
    +------+-------------+------------------+
    |  Bad | ``( a+b )`` | ``d [k], d[ k]`` |
    +------+-------------+------------------+

How should I test my code ?
---------------------------

There are two basic approaches for testing code in python: unit testing and doc testing. Their purpose is the same, to check that execution of code given some input produces a specified output. The cases to which the two approaches lend themselves are different.

An excellent discourse on testing code and the pros and cons of these alternatives is provided in a presentation by `Jim Fulton`_, which is recommended reading. A significant change since that presentation is that ``doctest`` can now read content that is not contained within docstrings. A another comparison of these two approaches, along with a third (``py.test``) is also available_. To see examples of both styles of testing look in ``PyCogent/tests``: files ending in .rst are using ``doctest``, those ending in .py are using ``unittest``.

.. _`Jim Fulton`: http://www.python.org/pycon/dc2004/papers/4/PyCon2004DocTestUnit.pdf
.. _available: http://agiletesting.blogspot.com/2005/11/articles-and-tutorials-page-updated.html

In general, it's easier to start writing ``doctest``'s, as you don't need to learn the ``unittest`` API but the latter give's much greater control.

Whatever approach is employed, the general principle is every line of code should be tested. It is critical that your code be fully tested before you draw conclusions from results it produces. For scientific work, bugs don't just mean unhappy users who you'll never actually meet: they may mean retracted publications.

Tests are an opportunity to invent the interface(s) you want. Write the test for a method before you write the method: often, this helps you figure out what you would want to call it and what parameters it should take. It's OK to write the tests a few methods at a time, and to change them as your ideas about the interface change. However, you shouldn't change them once you've told other people what the interface is.

Never treat prototypes as production code. It's fine to write prototype code without tests to try things out, but when you've figured out the algorithm and interfaces you must rewrite it *with tests* to consider it finished. Often, this helps you decide what interfaces and functionality you actually need and what you can get rid of.

"Code a little test a little". For production code, write a couple of tests, then a couple of methods, then a couple more tests, then a couple more methods, then maybe change some of the names or generalize some of the functionality. If you have a huge amount of code where 'all you have to do is write the tests', you're probably closer to 30% done than 90%. Testing vastly reduces the time spent debugging, since whatever went wrong has to be in the code you wrote since the last test suite. And remember to use python's interactive interpreter for quick checks of syntax and ideas.

Run the test suite when you change `anything`. Even if a change seems trivial, it will only take a couple of seconds to run the tests and then you'll be sure. This can eliminate long and frustrating debugging sessions where the change turned out to have been made long ago, but didn't seem significant at the time.

Some ``unittest`` pointers
^^^^^^^^^^^^^^^^^^^^^^^^^^

- *Use the* ``unittest`` *framework with tests in a separate file for each module.* Name the test file ``test_module_name.py``. Keeping the tests separate from the code reduces the temptation to change the tests when the code doesn't work, and makes it easy to verify that a completely new implementation presents the same interface (behaves the same) as the old.

- *Use* ``evo.unit_test`` *if you are doing anything with floating point numbers or permutations* (use ``assertFloatEqual``). Do *not* try to compare floating point numbers using ``assertEqual`` if you value your sanity. ``assertFloatEqualAbs`` and ``assertFloatEqualRel`` can specifically test for absolute and relative differences if the default behavior is not giving you what you want. Similarly, ``assertEqualItems``, ``assertSameItems``, etc. can be useful when testing permutations.

- *Test the interface of each class in your code by defining at least one* ``TestCase`` *with the name* ``ClassNameTests``. This should contain tests for everything in the public interface.

- *If the class is complicated, you may want to define additional tests with names* ``ClassNameTests_test_type``. These might subclass ``ClassNameTests`` in order to share ``setUp`` methods, etc.

- *Tests of private methods should be in a separate* ``TestCase`` *called* ``ClassNameTests_private``. Private methods may change if you change the implementation. It is not required that test cases for private methods pass when you change things (that's why they're private, after all), though it is often useful to have these tests for debugging.

- *Test `all` the methods in your class.* You should assume that any method you haven't tested has bugs. The convention for naming tests is ``test_method_name``. Any leading and trailing underscores on the method name can be ignored for the purposes of the test; however, *all tests must start with the literal substring* ``test`` *for* ``unittest`` *to find them.* If the method is particularly complex, or has several discretely different cases you need to check, use ``test_method_name_suffix``, e.g. ``test_init_empty``, ``test_init_single``, ``test_init_wrong_type``, etc. for testing ``__init__``.

- *Write good docstrings for all your test methods.* When you run the test with the ``-v`` command-line switch for verbose output, the docstring for each test will be printed along with ``...OK`` or ``...FAILED`` on a single line. It is thus important that your docstring is short and descriptive, and makes sense in this context.

    **Good docstrings:** ::

        NumberList.var should raise ValueError on empty or 1-item list
        NumberList.var should match values from R if list has >2 items
        NumberList.__init__ should raise error on values that fail float()
        FrequencyDistribution.var should match corresponding NumberList var

    **Bad docstrings:** ::

        var should calculate variance           # lacks class name, not descriptive
        Check initialization of a NumberList    # doesn't say what's expected
        Tests of the NumberList initialization. # ditto

- *Module-level functions should be tested in their own* ``TestCase``\ *, called* ``modulenameTests``. Even if these functions are simple, it's important to check that they work as advertised.

- *It is much more important to test several small cases that you can check by hand than a single large case that requires a calculator.* Don't trust spreadsheets for numerical calculations -- use R instead!

- *Make sure you test all the edge cases: what happens when the input is None, or '', or 0, or negative?* What happens at values that cause a conditional to go one way or the other? Does incorrect input raise the right exceptions? Can your code accept subclasses or superclasses of the types it expects? What happens with very large input?

- *To test permutations, check that the original and shuffled version are different, but that the sorted original and sorted shuffled version are the same.* Make sure that you get *different* permutations on repeated runs and when starting from different points.

- *To test random choices, figure out how many of each choice you expect in a large sample (say, 1000 or a million) using the binomial distribution or its normal approximation.* Run the test several times and check that you're within, say, 3 standard deviations of the mean.

Example of a ``unittest`` test module structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    #!/usr/bin/env python

    """Tests NumberList and FrequencyDistribution, classes for statistics."""

    from cogent.util.unit_test import TestCase, main # for floating point test use unittestfp
    from statistics import NumberList, FrequencyDistribution

    class NumberListTests(TestCase): # remember to subclass TestCase
        """Tests of the NumberList class."""
        def setUp(self):
            """Define a few standard NumberLists."""
            self.Null = NumberList()            # test empty init
            self.Empty = NumberList([])         # test init with empty sequence
            self.Single = NumberList([5])       # single item
            self.Zero = NumberList([0])         # single, False item
            self.Three = NumberList([1,2,3])    # multiple items
            self.ZeroMean = NumberList([1,-1])  # items nonzero, mean zero
            self.ZeroVar = NumberList([1,1,1])  # items nonzero, mean nonzero, variance zero
            # etc. These objects shared by all tests, and created new each time a method
            # starting with the string 'test' is called (i.e. the same object does not
            # persist between tests: rather, you get separate copies).

            def test_mean_empty(self):
                """NumberList.mean() should raise ValueError on empty object"""
                for empty in (self.Null, self.Empty):
                    self.assertRaises(ValueError, empty.mean)
            def test_mean_single(self):
                """NumberList.mean() should return item if only 1 item in list"""
                for single in (self.Single, self.Zero):
                    self.assertEqual(single.mean(), single[0])
            # other tests of mean
            def test_var_failures(self):
                """NumberList.var() should raise ZeroDivisionError if <2 items"""
                for small in (self.Null, self.Empty, self.Single, self.Zero):
                    self.assertRaises(ZeroDivisionError, small.var)
            # other tests of var
            # tests of other methods

    class FrequencyDistributionTests(TestCase):
        pass    # much code deleted
    # tests of other classes

    if __name__ == '__main__':    # run tests if called from command-line
        main()
