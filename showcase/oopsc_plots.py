"""
Showcase for the plot we reimplemented in oopsc
"""

import scanpy as sc
import oopsc.scanpy

from dash import Dash, html, dcc, callback, Output, Input

from joblib import Memory
memory = Memory("cachedir")

@memory.cache
def load_data():
    adata = sc.read_10x_mtx(
    'docs/tutorials/data/filtered_gene_bc_matrices/hg19',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)      
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    adata.var_names_make_unique()
    return adata


adata = load_data()

app = Dash(__name__)


app.layout = html.Div(
    [
        html.H1("oopsc showcase"),
        html.P("This dashboard shows the plot we reproduced in oopsc."),
        html.H2("Gene expression"),
        dcc.Graph(figure=sc.pli.umap(adata, color=['CST3', 'NKG7', 'PPBP'])),
        dcc.Graph(figure=sc.pli.highest_expr_genes(adata, n_top=20, )),
    ])

if __name__ == "__main__":
    app.run_server(debug=True)