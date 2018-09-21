{# This template was modified from autosummaries default format #}
{{ fullname | escape | underline}}

{# We need a list of the built-ins that we implemented, not the default ones #}
{% set built_in_methods = []  %}
{% for item in all_methods %}
   {% if (item not in ['__class__',
                        '__delattr__',
                        '__getattribute__',
                        '__init__',
                        '__dir__',
                        '__format__',
                        '__new__',
                        '__reduce__',
                        '__reduce_ex__',
                        '__repr__',
                        '__setattr__',
                        '__sizeof__',
                        '__subclasshook__'] and item.startswith('__')) %}
      {{ built_in_methods.append(item) or '' }}
   {% endif %}
{% endfor %}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}

   {% if built_in_methods %}
   .. rubric:: Built-ins

   .. autosummary::
      :toctree:
   {% for item in built_in_methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}

   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree:
   {% for item in methods %}
      {% if item != '__init__' %}
      ~{{ name }}.{{ item }}
      {% endif %}
   {%- endfor %}
   {% endif %}
