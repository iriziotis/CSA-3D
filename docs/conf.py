"""Sphinx configuration for CSA-3D documentation."""

import sys
from pathlib import Path

# Add parent directory to path for autodoc
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Project information
project = 'CSA-3D'
copyright = '2024, Ioannis Riziotis'
author = 'Ioannis Riziotis'
release = '0.3.1'

# Sphinx extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'myst_parser',
]

# MyST (Markdown) configuration
myst_enable_extensions = ["colon_fence"]
myst_heading_anchors = 3

# Napoleon settings for Google-style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_method = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_rtype = False

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

# Autodoc settings
autodoc_typehints = 'description'
autodoc_class_signature = 'separated'
autodoc_inherit_docstrings = True

# Mock imports for optional dependencies that might not be available on RTD
autodoc_mock_imports = [
    'rdkit',
    'rdkit.Chem',
    'rdkit.DataStructs',
]

# Intersphinx mappings
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'biopython': ('https://biopython.org/wiki/Documentation', None),
}

# Theme configuration
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': 'view',
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

# Source and build paths
templates_path = ['_templates']
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# HTML output settings
html_static_path = []
html_logo = None
html_favicon = None

# Misc settings
master_doc = 'index'
language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False

# Autodoc type hints
autodoc_typehints = 'description'

# Suppress warnings for missing references
suppress_warnings = ['ref.python']


# Setup function for autodoc
def setup(app):
    """Configure Sphinx extensions."""
    app.setup_extension('sphinx.ext.autodoc')
    # Ensure autodoc can find the package
    return {'version': release}
