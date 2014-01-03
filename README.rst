ordination
==========

A Python package to perform ordination analysis. The techniques that
are currently implemented are

* Correspondence Analysis
* Redundancy Analysis
* Canonical Correspondence Analysis

Dependencies
============

It depends on `numpy` and `matplotlib`, either on Python 2.7 (adding
compatibility with 2.6 would be easy) or 3.

Installation
============

As usual, running::

  $ python setup.py install

ought to work. You probably want to install it in some prefix, for
example by running `python setup.py install --user`.

Testing
=======

After installing the package and `nose`, run::

  $ nosetests --exe ordination

Usage
=====

The easiest way to get pretty pictures is to run it as a module::

  $ python -m ordination.base
