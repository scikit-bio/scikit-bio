{# This template is specifically for skbio.tree.TreeNode. #}

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

      ~{{ name }}.ancestors
      ~{{ name }}.has_children
      ~{{ name }}.is_root
      ~{{ name }}.is_tip
      ~{{ name }}.lca
      ~{{ name }}.neighbors
      ~{{ name }}.path
      ~{{ name }}.root
      ~{{ name }}.siblings

   .. rubric:: Tree traversal

   .. autosummary::
      :toctree:

      ~{{ name }}.levelorder
      ~{{ name }}.non_tips
      ~{{ name }}.postorder
      ~{{ name }}.pre_and_postorder
      ~{{ name }}.preorder
      ~{{ name }}.tips
      ~{{ name }}.traverse

   .. rubric:: Tree copying

   .. autosummary::
      :toctree:

      ~{{ name }}.copy
      ~{{ name }}.deepcopy
      ~{{ name }}.subtree

   .. rubric:: Tree manipulation

   .. autosummary::
      :toctree:

      ~{{ name }}.append
      ~{{ name }}.bifurcate
      ~{{ name }}.extend
      ~{{ name }}.insert
      ~{{ name }}.pop
      ~{{ name }}.prune
      ~{{ name }}.remove
      ~{{ name }}.remove_by_func
      ~{{ name }}.remove_deleted
      ~{{ name }}.shear
      ~{{ name }}.shuffle
      ~{{ name }}.unpack
      ~{{ name }}.unpack_by_func

   .. rubric:: Tree rerooting

   .. autosummary::
      :toctree:

      ~{{ name }}.root_at
      ~{{ name }}.root_at_midpoint
      ~{{ name }}.root_by_outgroup
      ~{{ name }}.unroot
      ~{{ name }}.unrooted_copy
      ~{{ name }}.unrooted_deepcopy
      ~{{ name }}.unrooted_move

   .. rubric:: Tree searching

   .. autosummary::
      :toctree:

      ~{{ name }}.assign_ids
      ~{{ name }}.cache_attr
      ~{{ name }}.clear_caches
      ~{{ name }}.create_caches
      ~{{ name }}.has_caches
      ~{{ name }}.find
      ~{{ name }}.find_all
      ~{{ name }}.find_by_func
      ~{{ name }}.find_by_id
      ~{{ name }}.index_tree

   .. rubric:: Tree analysis

   .. autosummary::
      :toctree:

      ~{{ name }}.bipart
      ~{{ name }}.biparts
      ~{{ name }}.cophenet
      ~{{ name }}.count
      ~{{ name }}.depth
      ~{{ name }}.distance
      ~{{ name }}.height
      ~{{ name }}.is_bifurcating
      ~{{ name }}.maxdist
      ~{{ name }}.observed_node_counts
      ~{{ name }}.subset
      ~{{ name }}.subsets
      ~{{ name }}.total_length

   .. rubric:: Tree comparison

   .. autosummary::
      :toctree:

      ~{{ name }}.compare_biparts
      ~{{ name }}.compare_cophenet
      ~{{ name }}.compare_rfd
      ~{{ name }}.compare_subsets
      ~{{ name }}.compare_wrfd

   .. rubric:: Tree visualization

   .. autosummary::
      :toctree:

      ~{{ name }}.ascii_art

   .. rubric:: Format conversion

   .. autosummary::
      :toctree:

      ~{{ name }}.assign_supports
      ~{{ name }}.from_linkage_matrix
      ~{{ name }}.from_taxdump
      ~{{ name }}.from_taxonomy
      ~{{ name }}.to_array
      ~{{ name }}.to_taxonomy

   {% endblock %}
