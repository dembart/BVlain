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
import os
import sys

# Define the canonical URL if you are using a custom domain on Read the Docs
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")

# Tell Jinja2 templates the build is running on Read the Docs
if os.environ.get("READTHEDOCS", "") == "True":
    if "html_context" not in globals():
        html_context = {}
    html_context["READTHEDOCS"] = True
    

sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(1, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------


project = 'bvlain'

copyright = '2024, Artem Dembitskiy'
author = 'Artem Dembitskiy'

# The full version, including alpha/beta/rc tags
release = '0.24.8'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
	"myst_parser",
	"sphinx.ext.duration",
	"sphinx.ext.autosectionlabel",
	"sphinx.ext.autodoc",
	"sphinx.ext.napoleon",
	"sphinx.ext.mathjax",
	"nbsphinx"
]

#add_module_names = False
autodoc_mock_imports = ["ase", "networkx", "pymatgen", "joblib", "numpy", "scipy"]

mathjax_config = {
    'TeX': {'equationNumbers': {'autoNumber': 'AMS', 'useLabelIds': True}},
}

source_suffix = {
    '.md': 'markdown',
}







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
html_theme = 'furo'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
