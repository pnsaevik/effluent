from datetime import datetime
from sphinx.application import Sphinx
from sphinx.util.docfields import Field
import os
import sys


def setup(app: Sphinx):
    app.add_object_type(
        'confval',
        'confval',
        objname='configuration value',
        indextemplate='pair: %s; configuration value',
        doc_field_types=[
            Field('type', label='Type', has_arg=False, names=('type',)),
            Field('default', label='Default', has_arg=False, names=('default',)),
            Field('units', label='Units', has_arg=False, names=('units',)),
        ]
    )


# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
sys.path.insert(0, os.path.abspath('../../src'))


# -- Project information -----------------------------------------------------

project = 'Effluent'
# noinspection PyShadowingBuiltins
copyright = f'{datetime.now().year}, Pål Næverlid Sævik'
author = 'Pål Næverlid Sævik'
source_suffix = '.rst'


# The full version, including alpha/beta/rc tags
def getversion():
    import effluent
    return effluent.__version__


release = getversion()


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'autoapi.extension',
    'matplotlib.sphinxext.plot_directive',
]

nitpicky = True
html_css_files = [
    'css/custom.css',
]

# Matplotlib extension options
plot_html_show_source_link = False
plot_formats = ['png']
plot_html_show_formats = False
plot_pre_code = ""

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for Sphinx AutoAPI -----------------------------------------------

autoapi_dirs = ['../../src']
autoapi_add_toctree_entry = True
autoapi_member_order = 'alphabetical'
autoapi_template_dir = '_templates/autoapi'
autoapi_options = [
    'members',
    'show-inheritance',
    'show-module-summary',
    'imported-members',
]
