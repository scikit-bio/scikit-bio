==========
ordination
==========

.. image:: https://travis-ci.org/Jorge-C/ordination.png?branch=master
    :target: https://travis-ci.org/Jorge-C/ordination

A Python package to perform ordination analysis. The techniques that
are currently implemented are

* Correspondence Analysis
* Redundancy Analysis
* Canonical Correspondence Analysis

Dependencies
============

It depends on ``numpy`` and ``matplotlib``, either on Python 2.7 (adding
compatibility with 2.6 would be easy) or 3.

Installation
============

As usual, running::

  $ python setup.py install

ought to work. You probably want to install it in some prefix, for
example by running ``python setup.py install --user``.

Testing
=======

After installing the package, run::

  $ nosetests --exe ordination

if you have ``nose`` installed or::

  $ py.test --pyargs ordination

if ``py.test`` is installed.

Usage
=====

The easiest way to get some pretty pictures is to run one of the
scripts in the examples folder::

  $ python simple.py
