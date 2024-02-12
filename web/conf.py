# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Configuration file for the Sphinx documentation builder.

from datetime import datetime


# -- Project information -----------------------------------------------------

project = 'scikit-bio'
author = f'{project} development team'
copyright = f'2014-{datetime.now().year}, {author}'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.extlinks',
    'sphinx_design',
    'sphinx_copybutton',
    'sphinxcontrib.youtube',
]

root_doc = 'index'
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

html_title = project
html_short_title = project
html_baseurl = 'https://scikit.bio'
html_logo = '_static/img/logo.svg'
html_favicon = '_static/favicon.ico'

# static files
html_static_path = ['_static']
html_css_files = ['css/style.css']

# do not show side bar with section navigation
html_sidebars = {"**": []}

# do not show source links
html_show_sourcelink = False


# -- External links ----------------------------------------------------------

github_url = f'https://github.com/{project}/{project}'
twitter_url = 'https://twitter.com/scikitbio'

extlinks = {
    'home': (f'{html_baseurl}/%s', None),
    'repo': (f'{github_url}/%s', None),
    'docs': (f'{html_baseurl}/docs/dev/generated/skbio.%s.html', None),
}


# -- PyData Theme configuration ----------------------------------------------

# References:
# https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/layout.html#
# references

html_theme = 'pydata_sphinx_theme'

html_theme_options = {

    # logo image for light/dark modes
    # image files must be placed under _static/
   'logo': {
      'alt_text': html_title,
      'image_light': '_static/img/logo.svg',
      'image_dark': '_static/img/logo_inv.svg',
    },

    # announcement banner on top of the screen
    'announcement': (
        f"{project} is back in active development! Check out our <a href='"
        f"{github_url}/discussions/1935'>announcement of revitalization</a>."
    ),

    # social media links displayed as icons
    'github_url': github_url,
    'twitter_url': twitter_url,

    # disable prev & next buttons
    'show_prev_next': False,

    # disable search button
    'navbar_persistent': [],

    # display all header links
    'header_links_before_dropdown': 7,

    # footer layout
    'footer_start': ['copyright'],
    'footer_center': ['sphinx-version'],
    'footer_end': ['theme-version'],

    # google analytics
    'analytics': {
        'google_analytics_id': 'UA-6636235-9',
    }
}
