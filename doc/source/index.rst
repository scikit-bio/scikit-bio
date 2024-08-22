:html_theme.sidebar_secondary.remove:

.. page style and classes

.. raw:: html

   <style>

      /* hide prev/next button */

     .prev-next-footer {
       display: none;
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
   :maxdepth: 1

   io
   sequence
   alignment
   tree
   diversity
   stats
   embedding
   table
   metadata
   workflow
   util
