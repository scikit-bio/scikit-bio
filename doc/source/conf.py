# NOTE: parts of this file were taken from scipy's doc/source/conf.py. See
# scikit-bio/licenses/scipy.txt for scipy's license.

import glob
import sys
import os
import types
import re

if sys.version_info.major != 3:
    raise RuntimeError("scikit-bio can only be used with Python 3. You are "
                       "currently running Python %d." % sys.version_info.major)

# Force matplotlib to not use any Xwindows backend.
import matplotlib
matplotlib.use('Agg')

import sphinx
import sphinx.ext.autosummary as autosummary

class NewAuto(autosummary.Autosummary):
    def get_items(self, names):
        # Camel to snake case from http://stackoverflow.com/a/1176023/579416
        first_cap_re = re.compile('(.)([A-Z][a-z]+)')
        all_cap_re = re.compile('([a-z0-9])([A-Z])')
        def fix_item(display_name, sig, summary, real_name):
            class_names = {
                'TreeNode': 'tree',
                'TabularMSA': 'msa'
            }

            class_name = real_name.split('.')[-2]
            if class_name in class_names:
                nice_name = class_names[class_name]
            else:
                s1 = first_cap_re.sub(r'\1_\2', class_name)
                nice_name = all_cap_re.sub(r'\1_\2', s1).lower()
                if len(nice_name) > 10:
                    nice_name = ''.join([e[0] for e in nice_name.split('_')])

            def fmt(string):
                count = string.count('%s')
                return string % tuple([nice_name] * count)

            specials = {
                '__eq__': fmt('%s1 == %s2'),
                '__ne__': fmt('%s1 != %s2'),
                '__gt__': fmt('%s1 > %s2'),
                '__lt__': fmt('%s1 < %s2'),
                '__ge__': fmt('%s1 >= %s2'),
                '__le__': fmt('%s1 <= %s2'),
                '__getitem__': fmt('%s[x]'),
                '__iter__': fmt('iter(%s)'),
                '__contains__': fmt('x in %s'),
                '__bool__': fmt('bool(%s)'),
                '__str__': fmt('str(%s)'),
                '__reversed__': fmt('reversed(%s)'),
                '__len__': fmt('len(%s)'),
                '__copy__': fmt('copy.copy(%s)'),
                '__deepcopy__': fmt('copy.deepcopy(%s)'),
            }
            if display_name in specials:
                prefixes = autosummary.get_import_prefixes_from_env(self.env)
                obj = autosummary.import_by_name(display_name,
                                                 prefixes=prefixes)
                # Filter out any slot_wrappers that work their way in (more below)
                if type(obj[1]).__name__ == 'wrapper_descriptor':
                    return None
                return specials[display_name], '', summary, real_name
            return display_name, sig, summary, real_name

        skip = ['__init_subclass__']

        items = []
        for item in super(NewAuto, self).get_items(names):
            if item[0] not in skip:
                temp_item = fix_item(*item)
                # Drop slot_wrappers (see above)
                if temp_item is not None:
                    items.append(temp_item)

        return items

autosummary.Autosummary = NewAuto

import sphinx_bootstrap_theme

# The extra_public_methods depends on what class we are looking at.

import skbio
from skbio.util._decorator import classproperty

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here:
#
#    sys.path.insert(0, os.path.abspath('../sphinxext/foo'))

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
# Using `sphinx_version` doesn't work, likely because Sphinx is expecting a
# version string of the form X.Y, not X.Y.Z.
needs_sphinx = '1.6'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.linkcode',
    'sphinx.ext.coverage',
    'sphinx.ext.doctest',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx'
]

# Determine if the matplotlib has a recent enough version of the
# plot_directive.

try:
    from matplotlib.sphinxext import plot_directive
except ImportError:
    use_matplotlib_plot_directive = False
else:
    try:
        use_matplotlib_plot_directive = (plot_directive.__version__ >= 2)
    except AttributeError:
        use_matplotlib_plot_directive = False

if use_matplotlib_plot_directive:
    extensions.append('matplotlib.sphinxext.plot_directive')
else:
    raise RuntimeError("You need a recent enough version of matplotlib")

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'scikit-bio'
copyright = u'2014--, scikit-bio development team'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = skbio.__version__
# The full version, including alpha/beta/rc tags.
release = skbio.__version__

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.

# Exclude this file since it is only used by autosummary to generate other RST
# files during the build process, and it will generate sphinx errors and
# warnings otherwise.
exclude_patterns = ['_templates/autosummary/*.rst']

# The reST default role (used for this markup: `text`) to use for all
# documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
#keep_warnings = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'bootstrap'
html_theme_path = sphinx_bootstrap_theme.get_html_theme_path()

# Theme options are theme-specific and customize the look and feel of a theme
# further. For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    # Navigation bar title. (Default: ``project`` value)
    'navbar_title': 'scikit-bio docs',

    # Render the next and previous page links in navbar. (Default: true)
    'navbar_sidebarrel': False,

    # Bootswatch (http://bootswatch.com/) theme.
    #
    # Options are nothing with "" (default) or the name of a valid theme
    # such as "amelia" or "cosmo".
    'bootswatch_theme': 'united',

    # Location of link to source.
    # Options are "nav" (default), "footer" or anything else to exclude.
    'source_link_position': False
}

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static/']

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
#html_extra_path = []

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
#html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'scikit-biodoc'


# -- Options for LaTeX output ---------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
  ('index', 'scikit-bio.tex', u'scikit-bio Documentation',
   u'scikit-bio development team', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'scikit-bio', u'scikit-bio Documentation',
     [u'scikit-bio development team'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'scikit-bio', u'scikit-bio Documentation',
   u'scikit-bio development team', 'scikit-bio',
   'Data structures, algorithms, and educational resources for working with '
   'biological data in Python.', 'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'

# If true, do not generate a @detailmenu in the "Top" node's menu.
#texinfo_no_detailmenu = False

# -- Options for autosummary ----------------------------------------------
autosummary_generate = glob.glob('*.rst')

#------------------------------------------------------------------------------
# Plot
#------------------------------------------------------------------------------
plot_pre_code = """
import numpy as np
import scipy as sp
np.random.seed(123)
"""
plot_include_source = True
plot_formats = [('png', 96), ]
#plot_html_show_formats = False

font_size = 13*72/96.0  # 13 px

plot_rcparams = {
    'font.size': font_size,
    'axes.titlesize': font_size,
    'axes.labelsize': font_size,
    'xtick.labelsize': font_size,
    'ytick.labelsize': font_size,
    'legend.fontsize': font_size,
    'figure.subplot.bottom': 0.2,
    'figure.subplot.left': 0.2,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.9,
    'figure.subplot.wspace': 0.2,
    'text.usetex': False,

    # Some of our figures have legends outside the axes area. When they're
    # rendered in an interactive context, nothing gets cut off, but when
    # rendered in a static context (e.g., with savefig, which the plot
    # directive uses), the legend can get cut off. Specifying 'tight' instead
    # of 'standard' fixes the issue. See http://stackoverflow.com/a/10154763
    'savefig.bbox': 'tight'
}

matplotlib.rcParams.update(plot_rcparams)

# -----------------------------------------------------------------------------
# Intersphinx configuration
# -----------------------------------------------------------------------------
intersphinx_mapping = {
        'http://docs.python.org/dev': None,
        'http://docs.scipy.org/doc/numpy': None,
        'http://docs.scipy.org/doc/scipy/reference': None,
        'http://matplotlib.org': None,
        'http://pandas.pydata.org/pandas-docs/stable': None,
        'http://www.biom-format.org':None
}

# -----------------------------------------------------------------------------
# Source code links
# -----------------------------------------------------------------------------

import inspect
from os.path import relpath, dirname

def linkcode_resolve(domain, info):
    """
    Determine the URL corresponding to Python object
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
        source, lineno = inspect.findsource(obj)
    except:
        lineno = None

    if lineno:
        linespec = "#L%d" % (lineno + 1)
    else:
        linespec = ""

    fn = relpath(fn, start=dirname(skbio.__file__))

    if 'dev' in skbio.__version__:
        return "http://github.com/biocore/scikit-bio/blob/master/skbio/%s%s" % (
           fn, linespec)
    else:
        return "http://github.com/biocore/scikit-bio/blob/%s/skbio/%s%s" % (
           skbio.__version__, fn, linespec)

#------------------------------------------------------------------------------
# linkcheck
#------------------------------------------------------------------------------

# Link-checking on Travis sometimes times out.
linkcheck_timeout = 30

# You might see the following exception when building the documentation:
# TypeError: 'abstractproperty' object is not iterable
def _closure():
    def __get__(self, cls, owner):
        return self

    classproperty.__get__ = __get__
_closure()

# Add the 'copybutton' javascript, to hide/show the prompt in code
# examples, originally taken from scikit-learn's doc/conf.py
def setup(app):
    app.add_javascript('copybutton.js')
    app.add_stylesheet('style.css')
