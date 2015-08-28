{% extends "!autosummary/class.rst" %}

{# Taken from scipy's sphinx documentation setup (https://github.com/scipy/scipy/blob/master/doc/source/_templates/autosummary/class.rst). #}

{% block methods %}
{% if methods %}
   .. HACK -- the point here is that we don't want this to appear in the output, but the autosummary should still generate the pages.
      .. autosummary::
         :toctree:
      {% for item in all_methods %}
         {# We want to build dunder methods if they exist, but not every kind of dunder. These are the dunders provided by default on `object` #}
         {%- if not item.startswith('_') or (item not in ['__class__',
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
         {{ name }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
{% endif %}
{% endblock %}

{% block attributes %}
{% if attributes %}
   .. HACK -- the point here is that we don't want this to appear in the output, but the autosummary should still generate the pages.
      .. autosummary::
         :toctree:
      {% for item in all_attributes %}
         {%- if not item.startswith('_') %}
         {{ name }}.{{ item }}
         {%- endif -%}
      {%- endfor %}
{% endif %}
{% endblock %}
