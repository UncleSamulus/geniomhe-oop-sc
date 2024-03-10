from __future__ import annotations

import numpy as np
import pandas as pd

from typing import TYPE_CHECKING

import scanpy as sc
import plotly.graph_objects as go
import plotly.express as px
from typing import Sequence, Iterable


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

def rank_genes_groups_violin(
    adata: AnnData,
    groups: Sequence[str] | None = None,
    *,
    n_genes: int = 20,
    gene_names: Iterable[str] | None = None,
    gene_symbols: str | None = None,
    use_raw: bool | None = None,
    key: str | None = None,
    split: bool = True,
    scale: str = "width",
    strip: bool = True,
    jitter: int | float | bool = True,
    size: int = 1,
    show: bool | None = None,
    save: bool | None = None,
 ):
    """\
    Plot ranking of genes for all tested comparisons.

    Parameters
    ----------
    adata
        Annotated data matrix.
    groups
        List of group names.
    n_genes
        Number of genes to show. Is ignored if `gene_names` is passed.
    gene_names
        List of genes to plot. Is only useful if interested in a custom gene list,
        which is not the result of :func:`scanpy.tl.rank_genes_groups`.
    gene_symbols
        Key for field in `.var` that stores gene symbols if you do not want to
        use `.var_names` displayed in the plot.
    use_raw
        Use `raw` attribute of `adata` if present. Defaults to the value that
        was used in :func:`~scanpy.tl.rank_genes_groups`.
    split
        Whether to split the violins or not.
    scale
        See :func:`~seaborn.violinplot`.
    strip
        Show a strip plot on top of the violin plot.
    jitter
        If set to 0, no points are drawn. See :func:`~seaborn.stripplot`.
    size
        Size of the jitter points.
    {show_save_ax}
    """
    # Scanpy code
    if key is None:
        key = "rank_genes_groups"
    groups_key = str(adata.uns[key]["params"]["groupby"])
    if use_raw is None:
        use_raw = bool(adata.uns[key]["params"]["use_raw"])
    reference = str(adata.uns[key]["params"]["reference"])
    groups_names = adata.uns[key]["names"].dtype.names if groups is None else groups
    if isinstance(groups_names, str):
        groups_names = [groups_names]
    fig = go.Figure()
    for group_name in groups_names:
        if gene_names is None:
            _gene_names = adata.uns[key]["names"][group_name][:n_genes]
        else:
            _gene_names = gene_names
        if isinstance(_gene_names, np.ndarray):
            _gene_names = _gene_names.tolist()
        df = sc.get.obs_df(adata, _gene_names, use_raw=use_raw, gene_symbols=gene_symbols)
        new_gene_names = df.columns
        df["hue"] = adata.obs[groups_key].astype(str).values
        if reference == "rest":
            df.loc[df["hue"] != group_name, "hue"] = "rest"
        else:
            df.loc[~df["hue"].isin([group_name, reference]), "hue"] = np.nan
        df["hue"] = df["hue"].astype("category")
        df_tidy = pd.melt(df, id_vars="hue", value_vars=new_gene_names)

        #oopsc contribution
        df_tidy2 = df_tidy[(df_tidy["hue"]=='0')|(df_tidy["hue"]=='1')] # Get rid of NaN 
        
        hue_order = [group_name, reference]
        for index, value in enumerate(df_tidy2['hue'].unique()):
            subset_df=df_tidy2[df_tidy2["hue"]==value]
            x= subset_df["variable"]
            y= subset_df["value"]
            side = 'negative' if index % 2 == 0 else 'positive'  # Alternate between positive and negative sides
            line_color = 'blue' if side == 'negative' else 'red'
            trace = go.Violin(x=x, 
                                    y=y, 
                                    legendgroup=gene_names, 
                                    scalegroup=gene_names, 
                                    name=gene_names,
                                    side=side,
                                    line_color=line_color)
            fig.add_trace(trace)
        
    fig.update_traces(meanline_visible=True)
    fig.update_layout(violingap=0, violinmode='overlay')
    return fig