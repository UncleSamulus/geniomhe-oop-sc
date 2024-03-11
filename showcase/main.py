"""A dashboard for single-cell quality control

ref. https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html"""

from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation
import oopsc.scanpy

sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=80,
    facecolor="white",
    frameon=False,
)

from joblib import Memory
memory = Memory("cachedir")
@memory.cache
def prepare_data():
    adata = sc.read_10x_h5(
        filename="filtered_feature_bc_matrix.h5",
        backup_url="https://figshare.com/ndownloader/files/39546196",
    )
    adata.var_names_make_unique()

    # Load the data
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    return adata

adata = prepare_data()

p1 = px.histogram(adata.obs, x="total_counts", nbins=100)
p2 = px.violin(adata.obs, y="pct_counts_mt")
p3 = px.scatter(adata.obs, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt")

# Create the app
app = Dash(__name__)

app.layout = html.Div(
    [
        html.H1("Single-cell quality control"),
        html.Div(
            children="This dashboard shows the quality control metrics for a single-cell RNA-seq dataset.",
        ),
        html.Div(
            children=[
                html.Div(
                    children=[
                        html.H2("Total counts"),
                        dcc.Graph(figure=p1),
                    ],
                    style={"width": "33%", "display": "inline-block"},
                ),
                html.Div(
                    children=[
                        html.H2("Percent mitochondrial counts"),
                        dcc.Graph(figure=p2),
                    ],
                    style={"width": "33%", "display": "inline-block"},
                ),
                html.Div(
                    children=[
                        html.H2("Total counts vs. number of genes"),
                        dcc.Graph(figure=p3),
                    ],
                    style={"width": "33%", "display": "inline-block"},
                ),
            ]
        ),
        # Add a footer with source
        html.Footer(
            children=[
                "Source: ",
                html.A(
                    "Single-cell Best Practices",
                    href="https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html",
                ),
            ]
        ),
    ]
)

if __name__ == "__main__":
    app.run_server(debug=True)
