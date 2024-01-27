# Configuration file for the Sphinx documentation builder.

import os
import sys

import skbio

# -- Project information -----------------------------------------------------

project = 'scikit-bio'
author = f'{project} development team'
copyright = f'2014--, {author}'
version = skbio.__version__
release = skbio.__version__


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.linkcode',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'matplotlib.sphinxext.plot_directive'
]

root_doc = 'index'
templates_path = ['_templates']
exclude_patterns = []

# -- Options for manual page output ------------------------------------------

# startdocname, name, description, authors, section
man_pages = [
    (root_doc, project, f'{project} Documentation', author, 1),
]

# -- Options for Texinfo output ----------------------------------------------

# startdocname, targetname, title, author, dir_entry, description, category,
# toctree_only
texinfo_documents = [
    (root_doc, project, f'{project} Documentation', author, project,
     'Data structures, algorithms, and educational resources for working with '
     'biological data in Python.', 'Science'),
]


# -- Options for HTML output -------------------------------------------------

# html_theme = 'alabaster'
# html_theme = 'bootstrap'
html_title = f'{project} {version} documentation'
html_short_title = project
html_baseurl = 'scikit.bio'
html_logo = 'assets/logo.png'
html_favicon = 'assets/favicon.ico'
html_static_path = ['_static']
htmlhelp_basename = 'skbio-doc'


# -- PyData Theme configuration ----------------------------------------------

html_theme = 'pydata_sphinx_theme'
html_theme_options = {
   'logo': {
      'image_light': 'assets/logo.png',
      'image_dark': 'assets/logo-inverted.png',
   }
}


# -- Intersphinx configuration -----------------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'pandas': ('https://pandas.pydata.org/pandas-docs/stable/', None),
    'biom-format': ('https://biom-format.org/', None)
}


# -- Source code links --------------------------------------------------------

import inspect
from os.path import relpath, dirname


def linkcode_resolve(domain, info):
    """Determine the URL corresponding to Python object.
    """

    if domain != 'py':
        return None

    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        return None

    obj = submod
    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except:
            return None

    try:
        fn = inspect.getsourcefile(obj)
    except:
        fn = None
    if not fn:
        try:
            fn = inspect.getsourcefile(sys.modules[obj.__module__])
        except:
            fn = None
    if not fn:
        return None

    try:
        _, lineno = inspect.findsource(obj)
    except:
        lineno = None

    branch = 'master' if 'dev' in version else version

    fn = relpath(fn, start=dirname(skbio.__file__))

    linespec = f'#L{lineno + 1}' if lineno else ''

    return ('https://github.com/scikit-bio/scikit-bio/blob/'
            f'{branch}/skbio/{fn}{linespec}')


# You might see the following exception when building the documentation:
# TypeError: 'abstractproperty' object is not iterable

from skbio.util._decorator import classproperty


def _closure():
    def __get__(self, cls, owner):
        return self

    classproperty.__get__ = __get__
_closure()
