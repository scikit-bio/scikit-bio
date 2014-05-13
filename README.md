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
    
    
---
[![Build Status](https://travis-ci.org/biocore/scikit-bio.png?branch=master)](https://travis-ci.org/biocore/scikit-bio)

Core objects, functions and statistics for working with biological data in Python.

To view scikit-bio's documentation, visit [scikit-bio.org](http://scikit-bio.org).

**Note:** scikit-bio is currently under active development and its API is not
stable. Major compatibility-breaking API changes will likely happen as
development progresses. Once the API has started to solidify and we have an
official release, we will strive to maintain backwards compatibility, provide
deprecation warnings, etc. wherever possible. In the meantime, feel free to try
out scikit-bio and let us know what you think!

Licensing
---------

scikit-bio is available under the new BSD license. See
[COPYING.txt](COPYING.txt) for scikit-bio's license, and the ```licenses```
directory for the licenses of other software that is (either partially or
entirely) distributed with scikit-bio.

See the [list of all of scikit-bio's contributors](https://github.com/biocore/scikit-bio/graphs/contributors).

The pre-history of scikit-bio
-----------------------------

scikit-bio began from code derived from [PyCogent](http://www.pycogent.org) and [QIIME](http://www.qiime.org), and the contributors and/or copyright holders have agreed to make the code they wrote for PyCogent and/or QIIME available under the BSD license. The contributors to PyCogent and/or QIIME modules that have been ported to scikit-bio are: Rob Knight (@rob-knight), Gavin Huttley (@gavin-huttley), Daniel McDonald (@wasade), Micah Hamady, Antonio Gonzalez (@antgonza), Sandra Smit, Greg Caporaso (@gregcaporaso), Jai Ram Rideout (@ElBrogrammer), Cathy Lozupone (@clozupone), Mike Robeson (@mikerobeson), Marcin Cieslik, Peter Maxwell, Jeremy Widmann, Zongzhi Liu, Michael Dwan, Logan Knecht (@loganknecht), Andrew Cochran, Jose Carlos Clemente (@cleme), Damien Coy, Levi McCracken, and Andrew Butterfield.

Installation
------------

In order to install sci-kit bio first install numpy

    pip install numpy

Then install scikit-bio

    pip install scikit-bio

If you're using a newer version of pip and scikit-bio fails to install, you may need to modify the pip command to be:

    pip install --allow-all-external --allow-unverified scikit-bio --process-dependency-links .

If you'd like to install the book's dependencies manually (or some other way than using pip), here's what you'll need:

- [Python](http://www.python.org/) 2.7
- [numpy](http://www.numpy.org/) >= 1.7

scikit-bio development
----------------------

If you're interested in getting involved in or learning about scikit-bio development, see [CONTRIBUTING.md](CONTRIBUTING.md), in this directory.


Logo
----
scikit-bio ASCII art created by [@gregcaporaso](https://github.com/gregcaporaso). Text logo created at [patorjk.com](http://patorjk.com/software/taag/).
