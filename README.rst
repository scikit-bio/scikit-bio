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

Core objects, functions and statistics for working with biological
data in Python.

To view scikit-bio's documentation, visit `scikit-bio.org
<http://scikit-bio.org>`__.

**Note:** scikit-bio is currently under active development and its
API is not stable. Major compatibility-breaking API changes will
likely happen as development progresses. Once the API has started to
solidify and we have an official release, we will strive to maintain
backwards compatibility, provide deprecation warnings, etc. wherever
possible. In the mean time, feel free to try out scikit-bio and let us
know what you think!


Licensing
---------

scikit-bio is available under the new BSD license. See
`COPYING.txt <https://github.com/biocore/scikit-bio/blob/master/COPYING.txt>`__ for scikit-bio's license, and the
``licenses`` directory for the licenses of other software that is
(either partially or entirely) distributed with scikit-bio.

Installation
------------

In order to install scikit-bio::

    pip install numpy
    pip install scikit-bio

If you'd like to install the dependencies manually (or some other way
than using pip), you can find those here:

-  `Python <http://www.python.org/>`__ 2.7
-  `numpy <http://www.numpy.org/>`__ >= 1.7
-  `scipy <http://www.scipy.org/>`__ >= 0.13.0
-  `matplotlib <http://www.matplotlib.org/>`__ >= 1.1.0
-  `pandas <http://pandas.pydata.org/>`__
-  `future <https://pypi.python.org/pypi/future>`__

Projects using scikit-bio
-------------------------

-  `QIIME <http://qiime.org/>`__
-  `EMPeror <http://biocore.github.io/emperor/>`__
-  `An Introduction to Applied
   Bioinformatics <http://caporasolab.us/An-Introduction-To-Applied-Bioinformatics/>`__

If you're using scikit-bio in your own projects, you should issue a
pull request adding them to this list.

scikit-bio development
----------------------

If you're interested in getting involved in or learning about
scikit-bio development, see `CONTRIBUTING.md <https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md>`__, in
this directory.

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
<https://github.com/justin212k>`__), and Jose Antonio Navas Molina
(`@josenavas <https://github.com/josenavas>`__).

Logo
----

scikit-bio ASCII art created by `@gregcaporaso
<https://github.com/gregcaporaso>`_. Text logo created at `patorjk.com
<http://patorjk.com/software/taag/>`__.

.. |Build Status| image:: https://travis-ci.org/biocore/scikit-bio.svg?branch=master
   :target: https://travis-ci.org/biocore/scikit-bio
.. |Coverage Status| image:: https://coveralls.io/repos/biocore/scikit-bio/badge.png
   :target: https://coveralls.io/r/biocore/scikit-bio
