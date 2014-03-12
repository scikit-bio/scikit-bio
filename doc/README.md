bipy documentation
==================

This guide contains instructions for building the bipy documentation, as well
as guidelines for contributing to the documentation.

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

Finally, you will need to install bipy.

**Important:** The documentation will be built for whatever version of bipy is
*currently installed* on your system (i.e., the version imported by
```import bipy```). This may not match the code located in this repository. You
will need to either install this version of bipy somewhere (e.g., in a
virtualenv) or point your ```PYTHONPATH``` environment variable to this code,
*before* building the documentation.

To build the documentation, assuming you are at the root of the bipy
repository:

    cd doc
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

bipy docstrings follow the [numpydoc conventions](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt).
This ensures that the docstrings are easily readable both from the interpreter
and HTML, PDF, etc.. Please read the numpydoc guidelines before continuing.

### Documenting a module in bipy

In addition to following the numpydoc conventions for docstrings, we have a few
more conventions that will ensure your documentation is correctly built and
linked within our website, and that it maintains consistency with the rest of
the bipy docs.

The easiest way to get started with documenting your code is to look at the
docstrings in existing bipy modules. A couple of modules to start with are
```bipy.core.sequence``` and ```bipy.core.distance```. Go ahead and look
through those now. We've structured our docs in a similar way to
[SciPy's documentation](http://docs.scipy.org/doc/scipy/reference/), so that
may be another good place to look for examples.

We'll take a top-down approach by discussing how to document a new module that
you'd like to add to bipy (let's call it ```bipy/core/example.py```).

#### Module docstring

The first thing you'll need to add is a docstring for the module. The docstring
should be the first thing in the file following the ```#!``` line. It should
start with a title for the module:

    #!/usr/bin/env python
    """
    Documentation examples (:mod:`bipy.core.example`)
    =================================================

It is important to include the ```:mod:``` Sphinx directive in the title, as
this title will be included in the table of contents. Also make sure that the
title underline is the same length as the title.

We also need to include another Sphinx directive below this:

    .. currentmodule:: bipy.core.example

This directive tells Sphinx that other classes, functions, etc. that we will
reference are located in the ```bipy.core.example``` module.

Next, include a more detailed description of the module. For example:

    This module consists of several example classes and functions to illustrate
    the bipy documentation system.

Following that, list any classes, functions, and exceptions that you'd like
documentation generated for. Note that you do *not* need to include every
single class, function, or exception that is defined in the module. Also, you
do not need to list class methods, as those will be automatically included in
the generated class documentation. Only include objects that should be exposed
as part of the public API.

For example:

    Classes
    -------

    .. autosummary::
       :toctree: generated/

       ExampleClass1
       ExampleClass2

    Functions
    ---------

    .. autosummary::
       :toctree: generated/

       example_function1
       example_function2

    Exceptions
    ----------

    .. autosummary::
       :toctree: generated/

       ExampleError

The ```autosummary``` directives are important as they generate RST files in
the ```generated/``` directory for each object. A single-line summary and link
to each object is inserted into the page for you.

After listing public module members, we encourage a usage example section
showing how to use some of the module's functionality. Examples should be
written in [doctest](http://docs.python.org/2/library/doctest.html) format so
that they can be automatically tested (e.g., using ```nosetests
--with-doctest``` or ```make doctest```).

    Examples
    --------

    Run the ``example_function1`` function:

    >>> from bipy.core.example import example_function1
    >>> example_function1("hello", "world")
    hello world!

You're now ready to document the members of your module.

#### Documenting module members

When documenting the members of a module (e.g., classes, methods, attributes,
functions, and exceptions), follow the numpydoc conventions. In addition to
these conventions, there are a few things to keep in mind:

- When documenting a class, only public methods and attributes are included in
  the built documentation by default. If a method or attribute starts with an
  underscore, it is assumed to be private. If you want a private method to be
  included in the built documentation, add the following line to the method's
  docstring:

    ```
    .. shownumpydoc
    ```

  For example, you might want to document "special" methods such as
  ```__getitem__```, ```__str__```, etc., which would be ignored by default. We
  recommend placing this at the end of the docstring for consistency. Note that
  this will only work for methods; private attributes will *always* be ignored.

- When documenting a class, include the ```Parameters``` and ```Attributes```
  sections in the class docstring, instead of in the ```__init__``` docstring.
  While numpydoc technically supports either form,
  ```__init__``` is not included in the list of methods and thus should have
  its documentation included in the class docstring.

#### Including the module in the docs

Until now, we've only been editing docstrings, which are attached to Python
code. The final step is to hook up this new module's docstrings to the
documentation build system:

1. Make sure you're within the ```bipy/doc``` directory.
2. Create a new file with the same name as your module under the ```source```
   directory. Do not include ```bipy``` as part of the name, and use ```.rst```
   as the suffix. For example, ```source/core.example.rst```.
3. Add the following line to ```source/core.example.rst``` to have your
   module's docstring pulled into the document:

    ```
    .. automodule:: bipy.core.example
    ```

4. Add the following line to ```source/index.rst``` to add the new page to the
   top-level table of contents:

    ```
    bipy.core.example
    ```

That's it! You can now try building the documentation, which should include the
documentation for your new module!

### Documenting a subpackage in bipy

The process of documenting a subpackage is very similar to documenting a module
in bipy. The only difference is that the module docstring goes in the
subpackage's ```__init__.py```.

### Troubleshooting

If things aren't working correctly, try running ```make clean``` and then
rebuild the docs. If things still aren't working, try building the docs
*without* your changes, and see if there are any Sphinx errors or warnings.
Make note of these, and then see what new errors or warnings are generated when
you add your changes again.
