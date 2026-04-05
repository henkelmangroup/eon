# -*- coding: utf-8 -*-
import sys, os
from datetime import date
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.append(os.path.abspath('.'))
###############################################################
####                PROJECT INFO
###############################################################
project = u"Eon"                #name
copyright = f"2010-{date.today().year},  "       # YYYY-YYYY, author
version = "1.1"                 # short verion X.Y
release = "1.1.2549"                 # full version, alpha/beta/rc
###############################################################
####                BUILD SUPPORT
###############################################################
#suffix of source files
source_suffix = {".txt": 'restructuredtext'}
#master Table of Contents tree document (root_doc alias)
master_doc = "index" # EMDW: changed from "contents" to "index"; doesn't seem to break anything
#directories to ignore
exclude_trees = ["_build"]
#mathjax support
extensions = ["sphinx.ext.autodoc","sphinx.ext.todo","sphinx.ext.coverage","sphinx.ext.mathjax"]
# Output file base name for HTML help builder.
htmlhelp_basename = 'eOndoc'
#custom templates
templates_path = ['_templates']


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"
#location of custom static and theme files
html_static_path = ["_static"]
html_theme_path = ["."]
#name of theme to use, either included in Sphinx or custom
html_theme = "henk_theme"
# The name for this set of Sphinx documents.  I
html_title = "EON: Long Timescale Dynamics"
# The name of an image file (relative to this directory) 
html_logo = "_static/eon_logo.png"
# The name of a Windows icon file (within the static path) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/favicon.ico"
html_css_files = ["henkv2.css"] # location relative to _static folder
html_show_copyright = True
html_show_sphinx = True             # have Created using Sphinx in footer
html_show_search_summary = True
html_show_sourcelink = True
html_sidebars = { # NOTE: order determines order will appear in sidebar
    "**" : [
        "localtoc.html",
        "sourcelink.html"
    ]
}
###############################################################
####                MATHJAX SUPPORT
###############################################################
html_math_renderer = "mathjax"          #default value
###############################################################
####                LATEX SUPPORT
###############################################################
latex_documents = [
  ('index', 'eOn.tex', u'eOn Documentation',
   u' ', 'manual'),
]
