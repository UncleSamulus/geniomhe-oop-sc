from __future__ import annotations

import numpy as np
import pandas as pd

from typing import TYPE_CHECKING

import scanpy as sc
import plotly.graph_objects as go
import plotly.express as px

# if TYPE_CHECKING:
from anndata import AnnData


def highest_expr_genes(
    adata: AnnData,
    n_top: int = 30,
    show: bool | None = None,
    save: str | bool | None = None,
    gene_symbols: str | None = None,
    log: bool = False,
    **kwds,
):
    """\
    Fraction of counts assigned to each gene over all cells.
    
    Parameters
    ----------
    adata
        Annotated data matrix
    n_top
        Number of top
    gene_symbols
        Key for field in .var that stores gene symbols
    log
        Plot x-axis in log scale
    **kwds
        Are passed to :func:`~plotly.express.box`
        Warning: this is probably not compatible with the :func:`~seaborn.boxplot` keyword arguments.
    
    Returns
    -------
    A plotly Fig    
    """
    # Scanpy code
    from scipy.sparse import issparse

    # compute the percentage of each gene per cell
    norm_dict = sc.pp.normalize_total(adata, target_sum=100, inplace=False)

    # identify the genes with the highest mean
    if issparse(norm_dict["X"]):
        mean_percent = norm_dict["X"].mean(axis=0).A1
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict["X"][:, top_idx].A
    else:
        mean_percent = norm_dict["X"].mean(axis=0)
        top_idx = np.argsort(mean_percent)[::-1][:n_top]
        counts_top_genes = norm_dict["X"][:, top_idx]
    columns = (
        adata.var_names[top_idx]
        if gene_symbols is None
        else adata.var[gene_symbols][top_idx]
    )
    counts_top_genes = pd.DataFrame(
        counts_top_genes, index=adata.obs_names, columns=columns
    )

    # oop-sc contribution
    fig = px.box(counts_top_genes, log_x=log, orientation="h", color="variable", **kwds)
    fig.update_layout(
        title="Highest expressed genes",
        xaxis_title="% of total counts",
        yaxis_title="Genes",
        showlegend=False
    )
    return fig