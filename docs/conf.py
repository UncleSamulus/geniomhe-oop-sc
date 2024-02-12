import sys
from pathlib import Path

from datetime import datetime
from packaging.version import Version
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sphinx.application import Sphinx


import oopsc

HERE = Path(__file__).parent

nitpicky = True  # Warn about broken links. This is here for a reason: Do not change.
needs_sphinx = "4.0"  # Nicer param docs
suppress_warnings = [
    "myst.header",  # https://github.com/executablebooks/MyST-Parser/issues/262
]

# General information
project = "oopsc - Scanpy Plotly"
author = "Océane Saibou & Samuel Ortion"
repository_url = "https://github.com/UncleSamulus/geniomhe-oop-sc"
copyright = f"{datetime.now():%Y}, the oopsc development team."
version = oopsc.__version__.replace(".dirty", "")

# Bumping the version updates all docs, so don't do that
if Version(version).is_devrelease:
    parsed = Version(version)
    version = f"{parsed.major}.{parsed.minor}.{parsed.micro}.dev"

release = version

# default settings
templates_path = ["_templates"]
master_doc = "index"
default_role = "literal"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**.ipynb_checkpoints"]

extensions = [
    "myst_nb",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    # "matplotlib.sphinxext.plot_directive",
    "sphinx_autodoc_typehints",  # needs to be after napoleon
    # "sphinx.ext.linkcode",
    "sphinx_design",
    "sphinx_search.extension",
    "sphinxext.opengraph",
    *[p.stem for p in (HERE / "extensions").glob("*.py") if p.stem not in {"git_ref"}],
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
todo_include_todos = False
api_dir = HERE / "api"  # function_images
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
    "html_admonition",
]
myst_url_schemes = ("http", "https", "mailto")
nb_output_stderr = "remove"
nb_execution_mode = "off"
nb_merge_streams = True


ogp_site_url = "https://unclesamulus.github.io/geniomhe-oop-sc/en/stable/"
# ogp_image = "https://unclesamulus.github.io/en/stable/_static/Scanpy_Logo_BrightFG.svg"

typehints_defaults = "braces"

pygments_style = "default"
pygments_dark_style = "native"

# intersphinx_mapping = dict(
#     anndata=("https://anndata.readthedocs.io/en/stable/", None),
#     bbknn=("https://bbknn.readthedocs.io/en/latest/", None),
#     cuml=("https://docs.rapids.ai/api/cuml/stable/", None),
#     cycler=("https://matplotlib.org/cycler/", None),
#     dask=("https://docs.dask.org/en/stable/", None),
#     dask_ml=("https://ml.dask.org/", None),
#     h5py=("https://docs.h5py.org/en/stable/", None),
#     ipython=("https://ipython.readthedocs.io/en/stable/", None),
#     igraph=("https://python.igraph.org/en/stable/api/", None),
#     leidenalg=("https://leidenalg.readthedocs.io/en/latest/", None),
#     louvain=("https://louvain-igraph.readthedocs.io/en/latest/", None),
#     matplotlib=("https://matplotlib.org/stable/", None),
#     networkx=("https://networkx.org/documentation/stable/", None),
#     numpy=("https://numpy.org/doc/stable/", None),
#     pandas=("https://pandas.pydata.org/pandas-docs/stable/", None),
#     pynndescent=("https://pynndescent.readthedocs.io/en/latest/", None),
#     pytest=("https://docs.pytest.org/en/latest/", None),
#     python=("https://docs.python.org/3", None),
#     rapids_singlecell=("https://rapids-singlecell.readthedocs.io/en/latest/", None),
#     scipy=("https://docs.scipy.org/doc/scipy/", None),
#     seaborn=("https://seaborn.pydata.org/", None),
#     sklearn=("https://scikit-learn.org/stable/", None),
#     tutorials=("https://scanpy-tutorials.readthedocs.io/en/latest/", None),
# )


# -- Options for HTML output ----------------------------------------------

# The theme is sphinx-book-theme, with patches for readthedocs-sphinx-search
html_theme = "sphinx_book_theme"
html_theme_options = {
    "repository_url": repository_url,
    "use_repository_button": True,
}
html_static_path = ["_static"]
html_show_sphinx = False
# html_logo = "_static/img/Scanpy_Logo_BrightFG.svg"
html_title = "oopsc - scanpy plotly"


def setup(app: Sphinx):
    """App setup hook."""
    app.add_config_value(
        "recommonmark_config",
        {
            "auto_toc_tree_section": "Contents",
            "enable_auto_toc_tree": True,
            "enable_math": True,
            "enable_inline_math": False,
            "enable_eval_rst": True,
        },
        True,
    )


# -- Options for other output formats ------------------------------------------

htmlhelp_basename = f"{project}doc"
doc_title = f"{project} Documentation"
latex_documents = [(master_doc, f"{project}.tex", doc_title, author, "manual")]
man_pages = [(master_doc, project, doc_title, [author], 1)]
texinfo_documents = [
    (
        master_doc,
        project,
        doc_title,
        author,
        project,
        "One line description of project.",
        "Miscellaneous",
    )
]


# -- Suppress link warnings ----------------------------------------------------

# Options for plot examples

plot_include_source = True
plot_formats = [("png", 90)]
plot_html_show_formats = False
plot_html_show_source_link = False
plot_working_directory = HERE.parent  # Project root

# extlinks config
extlinks = {
    "issue": ("https://github.com/UncleSamulus/geniomhe-oop-sc/issues/%s", "issue%s"),
    "pr": ("https://github.com/UncleSamulus/geniomhe-oop-sc/pull/%s", "pr%s"),
}
