API Reference
=============

Since ``effluent`` is intended for research applications, users are encouraged
to fork the code and implement modifications to suit their needs (still
:ref:`citing <citation>` the original source, of course).
A basic API reference is provided for users who want to familiarize themselves
with the source code.

Unlike the :doc:`configuration file <../config>`, the structure of the source
code documented in the API reference is not considered a part of the public
software interface. It may therefore change with new minor releases of the
package. In contrast, the structure of the configuration file is more or less
fixed. Breaking changes are not introduced unless there are very good reasons
to do so, in which case the major version number is updated.

The API reference is auto-generated from docstrings using
`sphinx-autoapi <https://github.com/readthedocs/sphinx-autoapi>`_. Each
submodule is listed below:

.. toctree::
   :titlesonly:
   :maxdepth: 1

   /autoapi/effluent/index
   {% for page in pages %}
   {% if page.type == "module" and page.display %}
   {{ page.include_path }}
   {% endif %}
   {% endfor %}

