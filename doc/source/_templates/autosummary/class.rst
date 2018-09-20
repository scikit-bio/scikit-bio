{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   .. automethod:: __init__

   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}

   {% endif %}


   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}

   .. rubric:: Built-ins

   .. autosummary::
      :toctree:
   {% for item in all_methods %}
      {# We want to build dunder methods if they exist, but not every kind of dunder. These are the dunders provided by default on `object` #}
      {%- if (item not in ['__class__',
                           '__delattr__',
                           '__getattribute__',
                           '__init__',
                           '__dir__',
                           '__new__',
                           '__reduce__',
                           '__reduce_ex__',
                           '__repr__',
                           '__setattr__',
                           '__sizeof__',
                           '__subclasshook__'] and item.startswith('__')) %}
      ~{{ name }}.{{ item }}
      {%- endif -%}
   {%- endfor %}

   {% endif %}
