{# This template was modified from autosummary's default template. #}
{{ fullname | escape | underline}}


{# Identify special methods (i.e., dunder or magic methods) that are relevant. #}

{% set methods_to_ignore = [
   '__class__',
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
   '__subclasshook__',
   '__init_subclass__',
   '__class_getitem__',
] %}
{% set special_methods = [] %}
{% for item in all_methods %}
   {% if item.startswith('__') and item not in methods_to_ignore %}
      {{ special_methods.append(item) or '' }}
   {% endif %}
{% endfor %}


{# Distinguish inherited and own members. #}

{% set inherited_attributes = [] %}
{% set own_attributes = [] %}
{% for item in attributes %}
   {% if item in inherited_members %}
      {{ inherited_attributes.append(item) or '' }}
   {% else %}
      {{ own_attributes.append(item) or '' }}
   {% endif %}
{% endfor %}

{% set inherited_methods = [] %}
{% set own_methods = [] %}
{% for item in methods %}
   {% if item != '__init__' %}
      {% if item in inherited_members %}
         {{ inherited_methods.append(item) or '' }}
      {% else %}
         {{ own_methods.append(item) or '' }}
      {% endif %}
   {% endif %}
{% endfor %}

{% set inherited_special_methods = [] %}
{% set own_special_methods = [] %}
{% for item in special_methods %}
   {% if item in inherited_members %}
      {{ inherited_special_methods.append(item) or '' }}
   {% else %}
      {{ own_special_methods.append(item) or '' }}
   {% endif %}
{% endfor %}


{# Class name, signatures, superclass (if any) and docstring. #}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}


{# Table of contents #}

   {% block attributes %}
   {% if own_attributes %}
   .. rubric:: Attributes

   .. autosummary::
   {% for item in own_attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}

   {% if inherited_attributes %}
   .. rubric:: Attributes (inherited)

   .. autoinherit::
   {% for item in inherited_attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% if own_methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree:
   {% for item in own_methods %}
      ~{{ name }}.{{ item }}
   {% endfor %}
   {% endif %}

   {% if inherited_methods %}
   .. rubric:: Methods (inherited)

   .. autoinherit::
   {%- for item in inherited_methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block special_methods %}
   {% if own_special_methods %}
   .. rubric:: Special methods

   .. autosummary::
   {% for item in own_special_methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}

   {% if inherited_special_methods %}
   .. rubric:: Special methods (inherited)

   .. autoinherit::
   {% for item in inherited_special_methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}


{# Detailed documentation #}

   {% block details %}
   .. rubric:: Details

   {% for item in own_attributes %}
   {% if item not in inherited_members %}
   .. autoattribute:: {{ name }}.{{ item }}
   {% endif %}
   {% endfor %}

   {% for item in own_special_methods %}
   {% if item not in inherited_members %}
   .. automethod:: {{ name }}.{{ item }}
   {% endif %}
   {% endfor %}
   {% endblock %}
