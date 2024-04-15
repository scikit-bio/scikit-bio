:html_theme.sidebar_secondary.remove:

.. page style and classes

.. raw:: html

   <style>

      /* hide prev/next button */

     .prev-next-footer {
       display: none;
     }

     /* arange level 2 entries in the same line */

     .toctree-l1 {
       margin-top: 20px;
       line-height: 2em;
     }

     .toctree-l1 > a:first-of-type {
       font-size: 1.5em;
       font-weight: bold;
     }

     .toctree-l2 {
       display: inline;
     }

     .toctree-l2:not(:last-of-type)::after {
       content: " - ";
     }

   </style>


.. hidden page title

.. title:: Home

.. toctree hidden from document but provides header links

.. toctree::
   :hidden:
   :maxdepth: 1

   Install <https://scikit.bio/install.html>
   Learn <https://scikit.bio/learn.html>
   Documentation <self>
   Contribute <https://scikit.bio/contribute.html>
   Community <https://github.com/scikit-bio/scikit-bio/discussions>
   Releases <https://github.com/scikit-bio/scikit-bio/blob/main/CHANGELOG.md>
   About <https://scikit.bio/about.html>


scikit-bio |version| documentation
==================================

scikit-bio (canonically pronounced *sigh-kit-buy-oh*) is a library for working with biological data in Python 3. scikit-bio is open source, BSD-licensed software that is currently under active development.

.. toctree::
   :maxdepth: 2

   io
   sequence
   alignment
   tree
   diversity
   stats
   table
   metadata
   workflow
   util
