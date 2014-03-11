bipy documentation
==================

This guide contains instructions for building the bipy documentation, as well
as guidelines for contributing to the documentation. Finally, lower-level
details are included for maintainers of the documentation system.

**Note:** If you're only interested in viewing the bipy documentation, visit
[bipy.org](http://bipy.org).

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

### Documentation guidelines

Most of bipy's API documentation is automatically generated from docstrings.
The advantage to this approach is that users can access the documentation in an
interactive Python session or from our website as HTML. Other output forms are
also possible, such as PDF.

bipy docstrings follow the [numpydoc conventions](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt). This ensures that the docstrings are easily readable both from the interpreter and HTML, PDF, etc.. Please read the numpydoc guidelines before continuing.

### Documenting a module in bipy

In addition to following the numpydoc conventions for docstrings, we have a few
more conventions that will ensure your documentation is correctly built and linked
within our website, and that it maintains consistency with the rest of the bipy
docs.

We'll take a top-down approach by discussing how to document a new
module that you'd like to add to bipy (let's call it ```bipy/core/example.py```).

The easiest way to get started with documenting your code is to look at the
docstrings in existing bipy modules. A couple of modules to start with are
```bipy.core.sequence``` and ```bipy.core.distance```. We recommend looking
through those now.

We highly recommend adding a module-level docstring.
