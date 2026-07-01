{# Template for skbio.tree.BPTree, a succinct balanced-parentheses tree backend. #}

{% extends "autosummary/class.rst" %}

   {% block methods %}

   .. rubric:: Tree IO

   .. autosummary::
      :toctree:

      ~{{ name }}.read
      ~{{ name }}.write

   .. rubric:: Tree navigation

   .. autosummary::
      :toctree:

      ~{{ name }}.root
      ~{{ name }}.parent
      ~{{ name }}.fchild
      ~{{ name }}.lchild
      ~{{ name }}.nsibling
      ~{{ name }}.psibling
      ~{{ name }}.lca
      ~{{ name }}.isancestor
      ~{{ name }}.deepestnode

   .. rubric:: Tree traversal

   .. autosummary::
      :toctree:

      ~{{ name }}.preorder
      ~{{ name }}.postorder
      ~{{ name }}.preorderselect
      ~{{ name }}.postorderselect
      ~{{ name }}.levelancestor
      ~{{ name }}.levelnext

   .. rubric:: Tree analysis

   .. autosummary::
      :toctree:

      ~{{ name }}.isleaf
      ~{{ name }}.depth
      ~{{ name }}.height
      ~{{ name }}.subtree
      ~{{ name }}.ntips

   .. rubric:: Tree manipulation

   .. autosummary::

      ~{{ name }}.shear
      ~{{ name }}.collapse

   {% endblock %}
