bipy documentation
==================

This guide contains instructions for building the bipy documentation, as well
as guidelines for contributing to the documentation. Finally, lower-level
details are included for maintainers of the documentation system.

**Note:** If you're only interested in viewing the bipy documentation, visit
[bipy.org][http://bipy.org].

Building the documentation
--------------------------

To build the documentation, you'll need the following Python packages
installed:

- [Sphinx](http://sphinx-doc.org/) >= 1.1
- [numpydoc](https://github.com/numpy/numpydoc) (latest master branch)
- [sphinx-bootstrap-theme](https://pypi.python.org/pypi/sphinx-bootstrap-theme/)

**Note:** At the time of this writing, the numpydoc version (0.4) on PyPI is a
bit outdated and some things don't link up correctly in the generated
documentation. You'll need to install the latest development version of
numpydoc that's available on GitHub.

An easy way to install the dependencies is via pip:

    pip install Sphinx sphinx-bootstrap-theme
    pip install git+git://github.com/numpy/numpydoc.git

Finally, you will need to install bipy. **IMPORTANT:** The documentation will
be built for whatever version of bipy is *currently installed* on your system
(i.e., the module imported by ```import bipy```). This may not match the code
located in this repository. You will need to either install this version of
bipy somewhere (e.g., in a virtualenv) or point your ```PYTHONPATH```
environment variable to this code, *before* building the documentation.

To build the documentation, assuming you are at the root of the bipy
repository:

    cd doc
    make clean
    make html

The built HTML documentation will be at ```build/html/index.html```.

Contributing to the documentation
---------------------------------

If you would like to contribute to the documentation, whether to add something
entirely new or to modify existing documentation, please first review our [bipy
contribution guide](../CONTRIBUTING.md).
