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

      ~{{ name }}.is_tip
      ~{{ name }}.is_root
      ~{{ name }}.has_children
      ~{{ name }}.root
      ~{{ name }}.ancestors
      ~{{ name }}.siblings
      ~{{ name }}.neighbors
      ~{{ name }}.lowest_common_ancestor
      ~{{ name }}.lca
      ~{{ name }}.path

   .. rubric:: Tree traversal

   .. autosummary::
      :toctree:

      ~{{ name }}.traverse
      ~{{ name }}.preorder
      ~{{ name }}.postorder
      ~{{ name }}.pre_and_postorder
      ~{{ name }}.levelorder
      ~{{ name }}.tips
      ~{{ name }}.non_tips

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
      ~{{ name }}.extend
      ~{{ name }}.insert
      ~{{ name }}.pop
      ~{{ name }}.remove
      ~{{ name }}.remove_by_func
      ~{{ name }}.remove_deleted
      ~{{ name }}.prune
      ~{{ name }}.shear
      ~{{ name }}.unpack
      ~{{ name }}.unpack_by_func
      ~{{ name }}.bifurcate
      ~{{ name }}.shuffle

   .. rubric:: Tree rerooting

   .. autosummary::
      :toctree:

      ~{{ name }}.unroot
      ~{{ name }}.unrooted_copy
      ~{{ name }}.unrooted_deepcopy
      ~{{ name }}.unrooted_move
      ~{{ name }}.root_at
      ~{{ name }}.root_at_midpoint
      ~{{ name }}.root_by_outgroup

   .. rubric:: Tree searching

   .. autosummary::
      :toctree:

      ~{{ name }}.has_caches
      ~{{ name }}.clear_caches
      ~{{ name }}.invalidate_caches
      ~{{ name }}.cache_attr
      ~{{ name }}.assign_ids
      ~{{ name }}.index_tree
      ~{{ name }}.create_caches
      ~{{ name }}.find
      ~{{ name }}.find_all
      ~{{ name }}.find_by_id
      ~{{ name }}.find_by_func

   .. rubric:: Tree analysis

   .. autosummary::
      :toctree:

      ~{{ name }}.count
      ~{{ name }}.subset
      ~{{ name }}.subsets
      ~{{ name }}.assign_supports
      ~{{ name }}.is_bifurcating
      ~{{ name }}.observed_node_counts
      ~{{ name }}.accumulate_to_ancestor
      ~{{ name }}.descending_branch_length
      ~{{ name }}.distance
      ~{{ name }}.get_max_distance
      ~{{ name }}.tip_tip_distances
      ~{{ name }}.compare_rfd
      ~{{ name }}.compare_subsets
      ~{{ name }}.compare_tip_distances

   .. rubric:: Tree visualization

   .. autosummary::
      :toctree:

      ~{{ name }}.ascii_art

   .. rubric:: Format conversion

   .. autosummary::
      :toctree:

      ~{{ name }}.from_linkage_matrix
      ~{{ name }}.from_taxonomy
      ~{{ name }}.to_taxonomy
      ~{{ name }}.from_taxdump
      ~{{ name }}.to_array

   {% endblock %}
