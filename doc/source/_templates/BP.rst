{# Template for skbio.tree.BPTree, a succinct balanced-parentheses tree backend. #}

{% extends "autosummary/class.rst" %}

   {% block methods %}

   .. rubric:: Tree IO

   .. autosummary::
      :toctree:

      ~{{ name }}.read
      ~{{ name }}.write
      ~{{ name }}.to_npz
      ~{{ name }}.from_npz

   .. rubric:: Tree navigation

   .. autosummary::
      :toctree:

      ~{{ name }}.root
      ~{{ name }}.parent
      ~{{ name }}.first_child
      ~{{ name }}.last_child
      ~{{ name }}.next_sibling
      ~{{ name }}.previous_sibling
      ~{{ name }}.lca
      ~{{ name }}.is_ancestor
      ~{{ name }}.deepest_node

   .. rubric:: Tree traversal

   .. autosummary::
      :toctree:

      ~{{ name }}.preorder_rank
      ~{{ name }}.postorder_rank
      ~{{ name }}.preorder_select
      ~{{ name }}.postorder_select
      ~{{ name }}.level_ancestor
      ~{{ name }}.level_next

   .. rubric:: Tree analysis

   .. autosummary::
      :toctree:

      ~{{ name }}.is_tip
      ~{{ name }}.depth
      ~{{ name }}.height
      ~{{ name }}.count

   .. rubric:: Tree manipulation

   .. autosummary::
      :toctree:

      ~{{ name }}.shear
      ~{{ name }}.collapse

   .. rubric:: Tree conversion

   .. autosummary::
      :toctree:

      ~{{ name }}.from_treenode
      ~{{ name }}.to_array

   {% endblock %}
