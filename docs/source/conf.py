# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

# sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../src/'))


# -- Project information -----------------------------------------------------

project = 'FunFAM'
copyright = '2020, Ramiro Checa-Garcia'
author = 'Ramiro Checa-Garcia'

# The full version, including alpha/beta/rc tags
release = '1.0'





# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [ 'sphinx.ext.graphviz',
               'sphinx.ext.autodoc',
               'sphinx.ext.doctest',
               'sphinx.ext.todo',
               'sphinx.ext.coverage',
               'sphinx.ext.imgmath',
               'sphinx.ext.viewcode',
               'sphinx.ext.napoleon',
               'sphinxcontrib.bibtex'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for LaTeX output ------------------------------------------------


#master_doc = "index"

#LATEX_STYLING = "_static/latex-styles.tex"
#try:
#    with open(LATEX_STYLING, mode="r+") as latex_styling:
#        PREAMBLE = latex_styling.read()
#except FileNotFoundError:
#    print("Could not read {0}".format(LATEX_STYLING), file=sys.stderr)
#    PREAMBLE = ""

#latex_elements = {"papersize": "a4paper", "pointsize": "10pt", "preamble": PREAMBLE}
#latex_engine = "xelatex"

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
#latex_documents = [(master_doc, filename + ".tex", project, author, "manual")]
