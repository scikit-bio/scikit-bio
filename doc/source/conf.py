# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Configuration file for the Sphinx documentation builder.

# Instructions:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime

import skbio


# -- Project information -----------------------------------------------------

project = 'scikit-bio'
author = f'{project} development team'
copyright = f'2014-{datetime.now().year}, {author}'
version = skbio.__version__
release = skbio.__version__


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.linkcode',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
    'numpydoc',
    'sphinx_design',
    'sphinx_copybutton',
    'matplotlib.sphinxext.plot_directive'
]

root_doc = 'index'
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for manual page output ------------------------------------------

# Elements: startdocname, name, description, authors, section
man_pages = [
    (root_doc, project, f'{project} Documentation', author, 1),
]


# -- Options for Texinfo output ----------------------------------------------

# Elements: startdocname, targetname, title, author, dir_entry, description,
# category, [toctree_only]
texinfo_documents = [
    (root_doc, project, f'{project} Documentation', author, project,
     'Data structures, algorithms, and educational resources for working with '
     'biological data in Python.', 'Science'),
]


# -- Options for HTML output -------------------------------------------------

html_title = f'{project} {version} documentation'
html_short_title = project
html_baseurl = 'https://scikit.bio'
html_logo = '_static/logo.svg'
html_favicon = '_static/favicon.ico'
htmlhelp_basename = 'skbio-doc'

# static files
html_static_path = ['_static']
html_css_files = ['css/style.css']

# do not show source links
html_show_sourcelink = False


# -- External links ----------------------------------------------------------

github_url = f'https://github.com/{project}/{project}'
twitter_url = 'https://twitter.com/scikitbio'

extlinks = {
    'home': (f'{html_baseurl}/%s', None),
    'repo': (f'{github_url}/%s', None),
}


# -- numpydoc configuration --------------------------------------------------

# References:
# https://numpydoc.readthedocs.io/en/latest/install.html#configuration

numpydoc_class_members_toctree = False
numpydoc_show_class_members = False
numpydoc_show_inherited_class_members = False


# -- PyData Theme configuration ----------------------------------------------

# References:
# https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html#
# references

html_theme = 'pydata_sphinx_theme'

html_theme_options = {

    # logo image for light/dark modes
    # image files must be placed under _static/
   'logo': {
      'link': html_baseurl,
      'alt_text': html_title,
      'image_light': '_static/logo.svg',
      'image_dark': '_static/logo_inv.svg',
    },

    # announcement banner on top of the screen
    # 'announcement': (
    #     f"{project} is back in active development! Check out our <a href='"
    #     f"{github_url}/discussions/1935'>announcement of revitalization</a>."
    # ),

    # social media links displayed as icons
    'github_url': github_url,
    'twitter_url': twitter_url,

    # show warning if not latest stable version
    # 'show_version_warning_banner': True,

    # version switcher
    'switcher': {
        'json_url': f'{html_baseurl}/versions.json',
        'version_match': 'dev' if version.endswith('-dev') else version,
    },

    # simplify section navigation
    'navigation_depth': 2,
    'collapse_navigation': True,

    # display all header links
    'header_links_before_dropdown': 7,

    # header layout
    'navbar_start': ['navbar-logo', 'version-switcher'],

    # footer layout
    'footer_start': ['copyright'],
    'footer_center': ['sphinx-version'],
    'footer_end': ['theme-version'],

    # google analytics
    'analytics': {
        'google_analytics_id': 'UA-6636235-9',
    }

}


# -- Intersphinx configuration -----------------------------------------------

# References:
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'pandas': ('https://pandas.pydata.org/pandas-docs/stable/', None),
    'biom-format': ('https://biom-format.org/', None)
}


# -- matplotlib.sphinxext.plot_directive -------------------------------------

# References:
# https://matplotlib.org/stable/api/sphinxext_plot_directive_api.html

plot_include_source = True
plot_html_show_source_link = False
plot_formats = ['png']
plot_html_show_formats = False


# -- Source code links --------------------------------------------------------

import inspect
from os.path import relpath, dirname


def linkcode_resolve(domain, info):
    """Determine the URL corresponding to a Python object.
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

    branch = 'main' if 'dev' in version else version

    fn = relpath(fn, start=dirname(skbio.__file__))

    linespec = f'#L{lineno + 1}' if lineno else ''

    return f'{github_url}/blob/{branch}/skbio/{fn}{linespec}'


# You might see the following exception when building the documentation:
# TypeError: 'abstractproperty' object is not iterable

from skbio.util._decorator import classproperty


def _closure():
    def __get__(self, cls, owner):
        return self

    classproperty.__get__ = __get__
_closure()
