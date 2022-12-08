..
  The extra blocks seem to be required for autosummary to populate our
  docstrings correctly. In the original template, these would recusively call
  autosummary again, but we already do that in the narratively convenient
  places.

.. automodule:: {{ fullname }}

   {% block attributes %}
   {% endblock %}

   {% block functions %}
   {% endblock %}

   {% block classes %}
   {% endblock %}

   {% block exceptions %}
   {% endblock %}
