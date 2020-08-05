
.. image:: http://scikit-bio.org/assets/logo.svg
   :target: http://scikit-bio.org
   :alt: scikit-bio logo

|Build Status| |Coverage Status| |ASV Benchmarks| |Gitter Badge| |Depsy Badge| |Anaconda Build Platforms| |Anaconda Build Version| |License| |Downloads| |Install|

scikit-bio is an open-source, BSD-licensed Python 3 package providing data structures, algorithms and educational resources for bioinformatics.

To view scikit-bio's documentation, visit `scikit-bio.org
<http://scikit-bio.org>`__.

**Note:** scikit-bio is no longer compatible with Python 2. scikit-bio is compatible with Python 3.6 and later.

scikit-bio is currently in beta. We are very actively developing it, and **backward-incompatible interface changes can and will arise**. To avoid these types of changes being a surprise to our users, our public APIs are decorated to make it clear to users when an API can be relied upon (stable) and when it may be subject to change (experimental). See the `API stability docs <https://github.com/biocore/scikit-bio/blob/master/doc/source/user/api_stability.rst>`_ for more details, including what we mean by *stable* and *experimental* in this context.

Installing
----------

The recommended way to install scikit-bio is via the ``conda`` package manager available in `Anaconda <http://continuum.io/downloads>`_ or `miniconda <http://conda.pydata.org/miniconda.html>`_.

To install the latest release of scikit-bio::

    conda install -c conda-forge scikit-bio

Alternatively, you can install scikit-bio using ``pip``::

    pip install scikit-bio

You can verify your installation by running the scikit-bio unit tests::

    python -m skbio.test

For users of Debian, ``skbio`` is in the Debian software distribution and may
be installed using::

    sudo apt-get install python3-skbio python-skbio-doc


Getting help
------------

To get help with scikit-bio, you should use the `skbio <http://stackoverflow.com/questions/tagged/skbio>`_ tag on StackOverflow (SO). Before posting a question, check out SO's guide on how to `ask a question <http://stackoverflow.com/questions/how-to-ask>`_. The scikit-bio developers regularly monitor the ``skbio`` SO tag.

Projects using scikit-bio
-------------------------

Some of the projects that we know of that are using scikit-bio are:

- `QIIME <http://qiime.org/>`__
- `Emperor <http://biocore.github.io/emperor/>`__
- `An Introduction to Applied
  Bioinformatics <http://readIAB.org>`__
- `tax2tree <https://github.com/biocore/tax2tree>`__
- `Qiita <http://qiita.microbio.me>`__
- `ghost-tree <https://github.com/JTFouquier/ghost-tree>`__
- `Platypus-Conquistador <https://github.com/biocore/Platypus-Conquistador>`__

If you're using scikit-bio in your own projects, feel free to issue a pull request to add them to this list.

scikit-bio development
----------------------

If you're interested in getting involved in scikit-bio development, see `CONTRIBUTING.md <https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md>`__.

See the list of `scikit-bio's contributors
<https://github.com/biocore/scikit-bio/graphs/contributors>`__.

Licensing
---------

scikit-bio is available under the new BSD license. See
`COPYING.txt <https://github.com/biocore/scikit-bio/blob/master/COPYING.txt>`__ for scikit-bio's license, and the
`licenses directory <https://github.com/biocore/scikit-bio/tree/master/licenses>`_ for the licenses of third-party software that is
(either partially or entirely) distributed with scikit-bio.

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
Ram Rideout (`@jairideout <https://github.com/jairideout>`__),
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

.. |Build Status| image:: https://travis-ci.org/biocore/scikit-bio.svg?branch=master
   :target: https://travis-ci.org/biocore/scikit-bio
.. |Coverage Status| image:: https://coveralls.io/repos/biocore/scikit-bio/badge.png
   :target: https://coveralls.io/r/biocore/scikit-bio
.. |ASV Benchmarks| image:: http://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat
   :target: https://s3-us-west-2.amazonaws.com/scikit-bio.org/benchmarks/master/index.html
.. |Gitter Badge| image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/biocore/scikit-bio
   :target: https://gitter.im/biocore/scikit-bio?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge
.. |Depsy Badge| image:: http://depsy.org/api/package/pypi/scikit-bio/badge.svg
   :target: http://depsy.org/package/python/scikit-bio
.. |Anaconda Build Platforms| image:: https://anaconda.org/conda-forge/scikit-bio/badges/platforms.svg
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |Anaconda Build Version| image:: https://anaconda.org/conda-forge/scikit-bio/badges/version.svg
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |License| image:: https://anaconda.org/conda-forge/scikit-bio/badges/license.svg
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |Downloads| image:: https://anaconda.org/conda-forge/scikit-bio/badges/downloads.svg
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |Install| image:: https://anaconda.org/conda-forge/scikit-bio/badges/installer/conda.svg
   :target: https://conda.anaconda.org/conda-forge
