::

               _ _    _ _          _     _
              (_) |  (_) |        | |   (_)
      ___  ___ _| | ___| |_ ______| |__  _  ___
     / __|/ __| | |/ / | __|______| '_ \| |/ _ \
     \__ \ (__| |   <| | |_       | |_) | | (_) |
     |___/\___|_|_|\_\_|\__|      |_.__/|_|\___/


           Opisthokonta
                   \  Amoebozoa
                    \ /
                     *    Euryarchaeota
                      \     |_ Crenarchaeota
                       \   *
                        \ /
                         *
                        /
                       /
                      /
                     *
                    / \
                   /   \
        Proteobacteria  \
                       Cyanobacteria

|Build Status| |Coverage Status|

scikit-bio is an open-source, BSD-licensed python package providing data structures, algorithms and educational resources for bioinformatics.

To view scikit-bio's documentation, visit `scikit-bio.org
<http://scikit-bio.org>`__.

scikit-bio is currently in alpha. We are very actively developing it, and **backwards-incompatible interface changes can and will arise**. Once the API has started to solidify, we will strive to maintain backwards compatibility. We will provide deprecation warnings wherever possible in the scikit-bio code, documentation, and CHANGELOG.md.

**Note:** Deprecation warnings will be issued using Python's ``DeprecationWarning`` class. Since Python 2.7, these types of warnings are **silenced by default**. When developing a tool that uses scikit-bio, we recommend enabling the display of deprecation warnings to be informed of upcoming API changes. For details on how to display deprecation warnings, see `Python's deprecation warning docs <https://docs.python.org/3/whatsnew/2.7.html#changes-to-the-handling-of-deprecation-warnings>`_.

Installation of release version (recommended for most users)
------------------------------------------------------------

To install the latest release version of scikit-bio you should run::

    pip install scikit-bio

Equivalently, you can use the ``conda`` package manager available in `Anaconda <http://continuum.io/downloads>`_ or `miniconda <http://conda.pydata.org/miniconda.html>`_ to install scikit-bio and all its dependencies, without having to compile them::

     conda install scikit-bio

Finally, most scikit-bio's dependencies (in particular, the ones that are trickier to build) are also available, albeit only for Python 2, in `Canopy Express <https://www.enthought.com/canopy-express/>`_.

You can verify your installation by running the scikit-bio unit tests as follows::

    nosetests --with-doctest skbio

Installation of development version
-----------------------------------

If you're interested in working with the latest development release of scikit-bio (recommended for developers only, as the development code can be unstable and less documented than the release code), you can clone the repository and install as follows. This will require that you have ``git`` installed.
::

    git clone git@github.com:biocore/scikit-bio.git
    cd scikit-bio
    pip install .

After this completes, you can run the scikit-bio unit tests from with a Python or IPython session as follows::

    import skbio
    skbio.test()

To only run tests for a specific scikit-bio module, for example skbio.sequence, you can do::

    skbio.sequence.test()

For developers of scikit-bio, if you don't want to be forced to re-install after every change, you can modify the above ``pip install`` command to::

    pip install -e .

This will build scikit-bio's Cython extensions, and will create a link in the ``site-packages`` directory to the scikit-bio source directory. When you then make changes to code in the source directory, those will be used (e.g., by the unit tests) without re-installing.

Finally, if you don't want to use ``pip`` to install scikit-bio, and prefer to just put ``scikit-bio`` in your ``$PYTHONPATH``, at the minimum you should run::

    python setup.py build_ext --inplace

This will build scikit-bio's Cython extensions, but not create a link to the scikit-bio source directory in ``site-packages``. If this isn't done, using certain components of scikit-bio will be inefficient and will produce an ``EfficiencyWarning``.

Getting help
------------

To get help with scikit-bio, you should use the `skbio <http://stackoverflow.com/questions/tagged/skbio>`_ tag on StackOverflow (SO). Before posting a question, check out SO's guide on how to `ask a question <http://stackoverflow.com/questions/how-to-ask>`_. The scikit-bio developers regularly monitor the skbio SO tag.

Licensing
---------

scikit-bio is available under the new BSD license. See
`COPYING.txt <https://github.com/biocore/scikit-bio/blob/master/COPYING.txt>`__ for scikit-bio's license, and the
`licenses directory <https://github.com/biocore/scikit-bio/tree/master/licenses>`_ for the licenses of third-party software that is
(either partially or entirely) distributed with scikit-bio.

Projects using scikit-bio
-------------------------

Some of the projects that we know of that are using scikit-bio are:

-  `QIIME <http://qiime.org/>`__
-  `Emperor <http://biocore.github.io/emperor/>`__
-  `An Introduction to Applied
   Bioinformatics <http://caporasolab.us/An-Introduction-To-Applied-Bioinformatics/>`__
-  `tax2tree <https://github.com/biocore/tax2tree>`__

If you're using scikit-bio in your own projects, you can issue a
pull request to add them to this list.

scikit-bio development
----------------------

If you're interested in getting involved in or learning about
scikit-bio development, see `CONTRIBUTING.md <https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md>`__.

See the `list of all of scikit-bio's contributors
<https://github.com/biocore/scikit-bio/graphs/contributors>`__.

Summaries of our weekly developer meetings are posted on
HackPad. Click `here
<https://hackpad.com/2014-scikit-bio-developer-meeting-notes-1S2RbMqy0iM>`__
to view the meeting notes for 2014.

The pre-history of scikit-bio
-----------------------------

scikit-bio began from code derived from `PyCogent
<http://www.pycogent.org>`__ and `QIIME <http://www.qiime.org>`__, and
the contributors and/or copyright holders have agreed to make the code
they wrote for PyCogent and/or QIIME available under the BSD
license. The contributors to PyCogent and/or QIIME modules that have
been ported to scikit-bio are: Rob Knight (`@rob-knight
<https://github.com/rob-knight>`__), Gavin Huttley (`@gavin-huttley
<https://github.com/gavin-huttley>`__), Daniel McDonald (`@wasade
<https://github.com/wasade>`__), Micah Hamady, Antonio Gonzalez
(`@antgonza <https://github.com/antgonza>`__), Sandra Smit, Greg
Caporaso (`@gregcaporaso <https://github.com/gregcaporaso>`__), Jai
Ram Rideout (`@ElBrogrammer <https://github.com/ElBrogrammer>`__),
Cathy Lozupone (`@clozupone <https://github.com/clozupone>`__), Mike Robeson
(`@mikerobeson <https://github.com/mikerobeson>`__), Marcin Cieslik,
Peter Maxwell, Jeremy Widmann, Zongzhi Liu, Michael Dwan, Logan Knecht
(`@loganknecht <https://github.com/loganknecht>`__), Andrew Cochran,
Jose Carlos Clemente (`@cleme <https://github.com/cleme>`__), Damien
Coy, Levi McCracken, Andrew Butterfield, Will Van Treuren (`@wdwvt1
<https://github.com/wdwvt1>`__), Justin Kuczynski (`@justin212k
<https://github.com/justin212k>`__), Jose Antonio Navas Molina
(`@josenavas <https://github.com/josenavas>`__), Matthew Wakefield
(`@genomematt <https://github.com/genomematt>`__) and Jens Reeder
(`@jensreeder <https://github.com/jensreeder>`__).

Logo
----

scikit-bio's logo was created by `Alina Prassas <http://cargocollective.com/alinaprassas>`_.
scikit-bio's ASCII art tree was created by `@gregcaporaso
<https://github.com/gregcaporaso>`_. Our text logo was created at `patorjk.com
<http://patorjk.com/software/taag/>`__.

.. |Build Status| image:: https://travis-ci.org/biocore/scikit-bio.svg?branch=master
   :target: https://travis-ci.org/biocore/scikit-bio
.. |Coverage Status| image:: https://coveralls.io/repos/biocore/scikit-bio/badge.png
   :target: https://coveralls.io/r/biocore/scikit-bio
