from __future__ import annotations

import collections.abc as cabc
import inspect
import sys
from collections.abc import Mapping, Sequence  # noqa: TCH003
from copy import copy
from functools import partial
from itertools import combinations, product
from numbers import Integral
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import pandas as pd
from anndata import AnnData  # noqa: TCH002
from cycler import Cycler  # noqa: TCH002
from matplotlib import colormaps, colors, patheffects, rcParams
from matplotlib import pyplot as plt
from matplotlib.axes import Axes  # noqa: TCH002
from matplotlib.colors import Colormap, Normalize
from matplotlib.figure import Figure  # noqa: TCH002
from numpy.typing import NDArray  # noqa: TCH002

import scanpy as sc

# For drop-in parameters compatibility with scanpy version
from scanpy.plotting._utils import (
    ColorLike, 
    VBound,
    _FontSize, 
    _FontWeight,
    _IGraphLayout,
    check_colornorm,
    check_projection,
    circles,
)

from scanpy._utils import (
    sanitize_anndata
)

from scanpy.plotting._tools.scatterplots import (
    _get_basis,
    _components_to_dimensions,
    _broadcast_args,
    _color_vector,
    _get_vboundnorm,
)


import plotly.express as px
import plotly.graph_objects as go


def embedding(
    adata: AnnData,
    basis: str,
    *,
    color: str | Sequence[str] | None = None,
    mask: NDArray[np.bool_] | str | None = None,
    gene_symbols: str | None = None,
    use_raw: bool | None = None,
    sort_order: bool = True,
    edges: bool = False,
    edges_width: float = 0.1,
    edges_color: str | Sequence[float] | Sequence[str] = "grey",
    neighbors_key: str | None = None,
    arrows: bool = False,
    arrows_kwds: Mapping[str, Any] | None = None,
    groups: str | Sequence[str] | None = None,
    components: str | Sequence[str] | None = None,
    dimensions: tuple[int, int] | Sequence[tuple[int, int]] | None = None,
    layer: str | None = None,
    projection: Literal["2d", "3d"] = "2d",
    scale_factor: float | None = None,
    color_map: Colormap | str | None = None,
    cmap: Colormap | str | None = None,
    palette: str | Sequence[str] | Cycler | None = None,
    na_color: ColorLike = "lightgray",
    na_in_legend: bool = True,
    size: float | Sequence[float] | None = None,
    frameon: bool | None = None,
    legend_fontsize: int | float | _FontSize | None = None,
    legend_fontweight: int | _FontWeight = "bold",
    legend_loc: str = "right margin",
    legend_fontoutline: int | None = None,
    colorbar_loc: str | None = "right",
    vmax: VBound | Sequence[VBound] | None = None,
    vmin: VBound | Sequence[VBound] | None = None,
    vcenter: VBound | Sequence[VBound] | None = None,
    norm: Normalize | Sequence[Normalize] | None = None,
    add_outline: bool | None = False,
    outline_width: tuple[float, float] = (0.3, 0.05),
    outline_color: tuple[str, str] = ("black", "white"),
    ncols: int = 4,
    hspace: float = 0.25,
    wspace: float | None = None,
    title: str | Sequence[str] | None = None,
    show: bool | None = None,
    save: bool | str | None = None,
    ax: Axes | None = None,
    return_fig: bool | None = None,
    marker: str | Sequence[str] = ".",
    **kwargs,
) -> go.Figure | None:
    """\
    Scatter plot for user specified embedding basis (e.g. umap, pca, etc)

    Parameters
    ----------
    basis
        Name of the `obsm` basis to use.
    # {adata_color_etc}
    # {edges_arrows}
    # {scatter_bulk}
    # {show_save_ax}

    Returns
    -------
    a :class:`~plotly.graph_objects.Figure`.
    """
    # scanpy
    #####################
    # Argument handling #
    #####################

    check_projection(projection)
    sanitize_anndata(adata)

    basis_values = _get_basis(adata, basis)
    dimensions = _components_to_dimensions(
        components, dimensions, projection=projection, total_dims=basis_values.shape[1]
    )
    args_3d = dict(projection="3d") if projection == "3d" else {}

    # Checking the mask format and if used together with groups
    if groups is not None and mask is not None:
        raise ValueError("Groups and mask arguments are incompatible.")
    if mask is not None:
        mask = sc.get._check_mask(adata, mask, "obs")

    # Figure out if we're using raw
    if use_raw is None:
        # check if adata.raw is set
        use_raw = layer is None and adata.raw is not None
    if use_raw and layer is not None:
        raise ValueError(
            "Cannot use both a layer and the raw representation. Was passed:"
            f"use_raw={use_raw}, layer={layer}."
        )
    if use_raw and adata.raw is None:
        raise ValueError(
            "`use_raw` is set to True but AnnData object does not have raw. "
            "Please check."
        )

    if isinstance(groups, str):
        groups = [groups]

    # Color map
    if color_map is not None:
        if cmap is not None:
            raise ValueError("Cannot specify both `color_map` and `cmap`.")
        else:
            cmap = color_map
    cmap = copy(colormaps.get_cmap(cmap))
    cmap.set_bad(na_color)
    kwargs["cmap"] = cmap
    # Prevents warnings during legend creation
    na_color = colors.to_hex(na_color, keep_alpha=True)

    if "edgecolor" not in kwargs:
        # by default turn off edge color. Otherwise, for
        # very small sizes the edge will not reduce its size
        # (https://github.com/scverse/scanpy/issues/293)
        kwargs["edgecolor"] = "none"

    # Vectorized arguments

    # turn color into a python list
    color = [color] if isinstance(color, str) or color is None else list(color)

    # turn marker into a python list
    marker = [marker] if isinstance(marker, str) else list(marker)

    if title is not None:
        # turn title into a python list if not None
        title = [title] if isinstance(title, str) else list(title)

    # turn vmax and vmin into a sequence
    if isinstance(vmax, str) or not isinstance(vmax, cabc.Sequence):
        vmax = [vmax]
    if isinstance(vmin, str) or not isinstance(vmin, cabc.Sequence):
        vmin = [vmin]
    if isinstance(vcenter, str) or not isinstance(vcenter, cabc.Sequence):
        vcenter = [vcenter]
    if isinstance(norm, Normalize) or not isinstance(norm, cabc.Sequence):
        norm = [norm]

    # Size
    if "s" in kwargs and size is None:
        size = kwargs.pop("s")
    if size is not None:
        # check if size is any type of sequence, and if so
        # set as ndarray
        if (
            size is not None
            and isinstance(size, (cabc.Sequence, pd.Series, np.ndarray))
            and len(size) == adata.shape[0]
        ):
            size = np.array(size, dtype=float)
    else:
        size = 120000 / adata.shape[0]

    ##########
    # Layout #
    ##########
    # Most of the code is for the case when multiple plots are required

    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams["figure.figsize"][0] + 0.02

    if components is not None:
        color, dimensions = list(zip(*product(color, dimensions)))

    color, dimensions, marker = _broadcast_args(color, dimensions, marker)

    # 'color' is a list of names that want to be plotted.
    # Eg. ['Gene1', 'louvain', 'Gene2'].
    # component_list is a list of components [[0,1], [1,2]]
    if (
        not isinstance(color, str)
        and isinstance(color, cabc.Sequence)
        and len(color) > 1
    ) or len(dimensions) > 1:
        if ax is not None:
            raise ValueError(
                "Cannot specify `ax` when plotting multiple panels "
                "(each for a given value of 'color')."
            )

    #     # each plot needs to be its own panel
    #     fig, grid = _panel_grid(hspace, wspace, ncols, len(color))
    # else:
    #     grid = None
    #     if ax is None:
    #         fig = plt.figure()
    #         ax = fig.add_subplot(111, **args_3d)

    ############
    # Plotting #
    ############
    # axs = []

    # use itertools.product to make a plot for each color and for each component
    # For example if color=[gene1, gene2] and components=['1,2, '2,3'].
    # The plots are: [
    #     color=gene1, components=[1,2], color=gene1, components=[2,3],
    #     color=gene2, components = [1, 2], color=gene2, components=[2,3],
    # ]
    for count, (value_to_plot, dims) in enumerate(zip(color, dimensions)):
        color_source_vector = _get_color_source_vector(
            adata,
            value_to_plot,
            layer=layer,
            mask=mask,
            use_raw=use_raw,
            gene_symbols=gene_symbols,
            groups=groups,
        )
        color_vector, categorical = _color_vector(
            adata,
            value_to_plot,
            values=color_source_vector,
            palette=palette,
            na_color=na_color,
        )

        # Order points
        order = slice(None)
        if sort_order is True and value_to_plot is not None and categorical is False:
            # Higher values plotted on top, null values on bottom
            order = np.argsort(-color_vector, kind="stable")[::-1]
        elif sort_order and categorical:
            # Null points go on bottom
            order = np.argsort(~pd.isnull(color_source_vector), kind="stable")
        # Set orders
        if isinstance(size, np.ndarray):
            size = np.array(size)[order]
        color_source_vector = color_source_vector[order]
        color_vector = color_vector[order]
        coords = basis_values[:, dims][order, :]

        # # if plotting multiple panels, get the ax from the grid spec
        # # else use the ax value (either user given or created previously)
        # if grid:
        #     ax = plt.subplot(grid[count], **args_3d)
        #     axs.append(ax)
        # if not (settings._frameon if frameon is None else frameon):
        #     ax.axis("off")
        # if title is None:
        #     if value_to_plot is not None:
        #         ax.set_title(value_to_plot)
        #     else:
        #         ax.set_title("")
        # else:
        #     try:
        #         ax.set_title(title[count])
        #     except IndexError:
        #         logg.warning(
        #             "The title list is shorter than the number of panels. "
        #             "Using 'color' value instead for some plots."
        #         )
        #         ax.set_title(value_to_plot)

        # if not categorical:
        #     vmin_float, vmax_float, vcenter_float, norm_obj = _get_vboundnorm(
        #         vmin, vmax, vcenter, norm=norm, index=count, colors=color_vector
        #     )
        #     normalize = check_colornorm(
        #         vmin_float,
        #         vmax_float,
        #         vcenter_float,
        #         norm_obj,
        #     )
        # else:
        #     normalize = None

        # make the scatter plot
        if projection == "3d":
            fig = px.scatter_3d(
                x=coords[:, 0],
                y=coords[:, 1],
                z=coords[:, 2],
                color=color_source_vector,
            )
        else:
            fig = px.scatter(
                x=coords[:, 0],
                y=coords[:, 1],
                color=color_source_vector,
            )
        # Remove color from legend
        fig.update_traces(showlegend=False)
        # Add annotations to hover data based on color_source_vector
        # else:
        #     scatter = (
        #         partial(ax.scatter, s=size, plotnonfinite=True)
        #         if scale_factor is None
        #         else partial(
        #             circles, s=size, ax=ax, scale_factor=scale_factor
        #         )  # size in circles is radius
        #     )

    #         if add_outline:
    #             # the default outline is a black edge followed by a
    #             # thin white edged added around connected clusters.
    #             # To add an outline
    #             # three overlapping scatter plots are drawn:
    #             # First black dots with slightly larger size,
    #             # then, white dots a bit smaller, but still larger
    #             # than the final dots. Then the final dots are drawn
    #             # with some transparency.

    #             bg_width, gap_width = outline_width
    #             point = np.sqrt(size)
    #             gap_size = (point + (point * gap_width) * 2) ** 2
    #             bg_size = (np.sqrt(gap_size) + (point * bg_width) * 2) ** 2
    #             # the default black and white colors can be changes using
    #             # the contour_config parameter
    #             bg_color, gap_color = outline_color

    #             # remove edge from kwargs if present
    #             # because edge needs to be set to None
    #             kwargs["edgecolor"] = "none"

    #             # remove alpha for outline
    #             alpha = kwargs.pop("alpha") if "alpha" in kwargs else None

    #             ax.scatter(
    #                 coords[:, 0],
    #                 coords[:, 1],
    #                 s=bg_size,
    #                 c=bg_color,
    #                 rasterized=settings._vector_friendly,
    #                 norm=normalize,
    #                 marker=marker[count],
    #                 **kwargs,
    #             )
    #             ax.scatter(
    #                 coords[:, 0],
    #                 coords[:, 1],
    #                 s=gap_size,
    #                 c=gap_color,
    #                 rasterized=settings._vector_friendly,
    #                 norm=normalize,
    #                 marker=marker[count],
    #                 **kwargs,
    #             )
    #             # if user did not set alpha, set alpha to 0.7
    #             kwargs["alpha"] = 0.7 if alpha is None else alpha

    #         cax = scatter(
    #             coords[:, 0],
    #             coords[:, 1],
    #             c=color_vector,
    #             rasterized=settings._vector_friendly,
    #             norm=normalize,
    #             marker=marker[count],
    #             **kwargs,
    #         )

    #     # remove y and x ticks
    #     ax.set_yticks([])
    #     ax.set_xticks([])
    #     if projection == "3d":
    #         ax.set_zticks([])

    #     # set default axis_labels
    #     name = _basis2name(basis)
    #     axis_labels = [name + str(d + 1) for d in dims]

    #     ax.set_xlabel(axis_labels[0])
    #     ax.set_ylabel(axis_labels[1])
    #     if projection == "3d":
    #         # shift the label closer to the axis
    #         ax.set_zlabel(axis_labels[2], labelpad=-7)
    #     ax.autoscale_view()

    #     if edges:
    #         _utils.plot_edges(
    #             ax, adata, basis, edges_width, edges_color, neighbors_key=neighbors_key
    #         )
    #     if arrows:
    #         _utils.plot_arrows(ax, adata, basis, arrows_kwds)

    #     if value_to_plot is None:
    #         # if only dots were plotted without an associated value
    #         # there is not need to plot a legend or a colorbar
    #         continue

    #     if legend_fontoutline is not None:
    #         path_effect = [
    #             patheffects.withStroke(linewidth=legend_fontoutline, foreground="w")
    #         ]
    #     else:
    #         path_effect = None

    #     # Adding legends
    #     if categorical or color_vector.dtype == bool:
    #         _add_categorical_legend(
    #             ax,
    #             color_source_vector,
    #             palette=_get_palette(adata, value_to_plot),
    #             scatter_array=coords,
    #             legend_loc=legend_loc,
    #             legend_fontweight=legend_fontweight,
    #             legend_fontsize=legend_fontsize,
    #             legend_fontoutline=path_effect,
    #             na_color=na_color,
    #             na_in_legend=na_in_legend,
    #             multi_panel=bool(grid),
    #         )
    #     elif colorbar_loc is not None:
    #         plt.colorbar(
    #             cax, ax=ax, pad=0.01, fraction=0.08, aspect=30, location=colorbar_loc
    #         )

    return fig

# API


def umap(adata: AnnData, **kwargs) -> Figure | Axes | list[Axes] | None:
    """\
    Scatter plot in UMAP basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        import oopsc.scanpy
        adata = sc.datasets.pbmc68k_reduced()
        sc.pli.umap(adata)

    Colour points by discrete variable (Louvain clusters).

    .. plot::
        :context: close-figs

        sc.pli.umap(adata, color="louvain")

    Colour points by gene expression.

    .. plot::
        :context: close-figs

        sc.pli.umap(adata, color="HES4")

    Plot muliple umaps for different gene expressions.

    .. plot::
        :context: close-figs

        sc.pli.umap(adata, color=["HES4", "TNFRSF4"])

    .. currentmodule:: oopsc

    See also
    --------
    """
    return embedding(adata, "umap", **kwargs)


def tsne(adata: AnnData, **kwargs) -> Figure | Axes | list[Axes] | None:
    """\
    Scatter plot in tSNE basis.

    Parameters
    ----------
    {adata_color_etc}
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.tsne(adata)
        sc.pl.tsne(adata, color='bulk_labels')

    .. currentmodule:: scanpy

    See also
    --------
    tl.tsne
    """
    return embedding(adata, "tsne", **kwargs)


def diffmap(adata: AnnData, **kwargs) -> Figure | Axes | list[Axes] | None:
    """\
    Scatter plot in Diffusion Map basis.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.diffmap(adata)
        sc.pl.diffmap(adata, color='bulk_labels')

    .. currentmodule:: scanpy

    See also
    --------
    tl.diffmap
    """
    return embedding(adata, "diffmap", **kwargs)


def draw_graph(
    adata: AnnData, *, layout: _IGraphLayout | None = None, **kwargs
) -> Figure | Axes | list[Axes] | None:
    """\
    Scatter plot in graph-drawing basis.

    Parameters
    ----------
    {adata_color_etc}
    layout
        One of the :func:`~scanpy.tl.draw_graph` layouts.
        By default, the last computed layout is used.
    {edges_arrows}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------
    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.draw_graph(adata)
        sc.pl.draw_graph(adata, color=['phase', 'bulk_labels'])

    .. currentmodule:: scanpy

    See also
    --------
    tl.draw_graph
    """
    if layout is None:
        layout = str(adata.uns["draw_graph"]["params"]["layout"])
    basis = f"draw_graph_{layout}"
    if f"X_{basis}" not in adata.obsm_keys():
        raise ValueError(
            f"Did not find {basis} in adata.obs. Did you compute layout {layout}?"
        )

    return embedding(adata, basis, **kwargs)

def pca(
    adata: AnnData,
    *,
    annotate_var_explained: bool = False,
    show: bool | None = None,
    return_fig: bool | None = None,
    save: bool | str | None = None,
    **kwargs,
) -> Figure | Axes | list[Axes] | None:
    """\
    Scatter plot in PCA coordinates.

    Use the parameter `annotate_var_explained` to annotate the explained variance.

    Parameters
    ----------
    {adata_color_etc}
    annotate_var_explained
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc3k_processed()
        sc.pli.pca(adata)

    Colour points by discrete variable (Louvain clusters).

    .. plot::
        :context: close-figs

        sc.pli.pca(adata, color="louvain")

    Colour points by gene expression.

    .. plot::
        :context: close-figs

        sc.pli.pca(adata, color="CST3")

    .. currentmodule:: oopsc

    See also
    --------
    pp.pca
    """
    if not annotate_var_explained:
        return embedding(
            adata, "pca", show=show, return_fig=return_fig, save=save, **kwargs
        )
    if "pca" not in adata.obsm.keys() and "X_pca" not in adata.obsm.keys():
        raise KeyError(
            f"Could not find entry in `obsm` for 'pca'.\n"
            f"Available keys are: {list(adata.obsm.keys())}."
        )

    label_dict = {
        f"PC{i + 1}": f"PC{i + 1} ({round(v * 100, 2)}%)"
        for i, v in enumerate(adata.uns["pca"]["variance_ratio"])
    }

    if return_fig is True:
        # edit axis labels in returned figure
        fig = embedding(adata, "pca", return_fig=return_fig, **kwargs)
        for ax in fig.axes:
            if xlabel := label_dict.get(ax.xaxis.get_label().get_text()):
                ax.set_xlabel(xlabel)
            if ylabel := label_dict.get(ax.yaxis.get_label().get_text()):
                ax.set_ylabel(ylabel)
        return fig

    # get the axs, edit the labels and apply show and save from user
    fig = embedding(adata, "pca", show=False, save=False, **kwargs)
    return fig

'''
def spatial(
    adata: AnnData,
    *,
    basis: str = "spatial",
    img: np.ndarray | None = None,
    img_key: str | None | Empty = _empty,
    library_id: str | None | Empty = _empty,
    crop_coord: tuple[int, int, int, int] | None = None,
    alpha_img: float = 1.0,
    bw: bool | None = False,
    size: float = 1.0,
    scale_factor: float | None = None,
    spot_size: float | None = None,
    na_color: ColorLike | None = None,
    show: bool | None = None,
    return_fig: bool | None = None,
    save: bool | str | None = None,
    **kwargs,
) -> go.Figure :
    """\
    Scatter plot in spatial coordinates.

    This function allows overlaying data on top of images.
    Use the parameter `img_key` to see the image in the background
    And the parameter `library_id` to select the image.
    By default, `'hires'` and `'lowres'` are attempted.

    Use `crop_coord`, `alpha_img`, and `bw` to control how it is displayed.
    Use `size` to scale the size of the Visium spots plotted on top.

    As this function is designed to for imaging data, there are two key assumptions
    about how coordinates are handled:

    1. The origin (e.g `(0, 0)`) is at the top left – as is common convention
    with image data.

    2. Coordinates are in the pixel space of the source image, so an equal
    aspect ratio is assumed.

    If your anndata object has a `"spatial"` entry in `.uns`, the `img_key`
    and `library_id` parameters to find values for `img`, `scale_factor`,
    and `spot_size` arguments. Alternatively, these values be passed directly.

    Parameters
    ----------
    {adata_color_etc}
    {scatter_spatial}
    {scatter_bulk}
    {show_save_ax}

    Returns
    -------
    If `show==False` a :class:`~matplotlib.axes.Axes` or a list of it.

    Examples
    --------
    This function behaves very similarly to other embedding plots like
    :func:`~scanpy.pl.umap`

    >>> import scanpy as sc
    >>> adata = sc.datasets.visium_sge("Targeted_Visium_Human_Glioblastoma_Pan_Cancer")
    >>> sc.pp.calculate_qc_metrics(adata, inplace=True)
    >>> sc.pl.spatial(adata, color="log1p_n_genes_by_counts")

    See Also
    --------
    :func:`scanpy.datasets.visium_sge`
        Example visium data.
    :doc:`tutorials:spatial/basic-analysis`
        Tutorial on spatial analysis.
    """
    # get default image params if available
    library_id, spatial_data = sc.pl._scatter._check_spatial_data(adata.uns, library_id)
    img, img_key = sc.plotting._scatter._check_img(spatial_data, img, img_key, bw=bw)
    spot_size = sc.plottng._scatter._check_spot_size(spatial_data, spot_size)
    scale_factor = _check_scale_factor(
        spatial_data, img_key=img_key, scale_factor=scale_factor
    )
    crop_coord = _check_crop_coord(crop_coord, scale_factor)
    na_color = _check_na_color(na_color, img=img)

    if bw:
        cmap_img = "gray"
    else:
        cmap_img = None
    circle_radius = size * scale_factor * spot_size * 0.5

    axs = embedding(
        adata,
        basis=basis,
        scale_factor=scale_factor,
        size=circle_radius,
        na_color=na_color,
        show=False,
        save=False,
        **kwargs,
    )
    if not isinstance(axs, list):
        axs = [axs]
    for ax in axs:
        cur_coords = np.concatenate([ax.get_xlim(), ax.get_ylim()])
        if img is not None:
            ax.imshow(img, cmap=cmap_img, alpha=alpha_img)
        else:
            ax.set_aspect("equal")
            ax.invert_yaxis()
        if crop_coord is not None:
            ax.set_xlim(crop_coord[0], crop_coord[1])
            ax.set_ylim(crop_coord[3], crop_coord[2])
        else:
            ax.set_xlim(cur_coords[0], cur_coords[1])
            ax.set_ylim(cur_coords[3], cur_coords[2])
    _utils.savefig_or_show("show", show=show, save=save)
    if return_fig:
        return axs[0].figure
    show = settings.autoshow if show is None else show
    if show:
        return None
    return axs
'''


"""Ugly fix, present in later version of scanpy"""


def _get_color_source_vector(
    adata: AnnData,
    value_to_plot: str,
    *,
    mask: NDArray[np.bool_] | None = None,
    use_raw: bool = False,
    gene_symbols: str | None = None,
    layer: str | None = None,
    groups: Sequence[str] | None = None,
) -> np.ndarray | pd.api.extensions.ExtensionArray:
    """
    Get array from adata that colors will be based on.
    """
    if value_to_plot is None:
        # Points will be plotted with `na_color`. Ideally this would work
        # with the "bad color" in a color map but that throws a warning. Instead
        # _color_vector handles this.
        # https://github.com/matplotlib/matplotlib/issues/18294
        return np.broadcast_to(np.nan, adata.n_obs)
    if (
        gene_symbols is not None
        and value_to_plot not in adata.obs.columns
        and value_to_plot not in adata.var_names
    ):
        # We should probably just make an index for this, and share it over runs
        # TODO: Throw helpful error if this doesn't work
        value_to_plot = adata.var.index[adata.var[gene_symbols] == value_to_plot][0]
    if use_raw and value_to_plot not in adata.obs.columns:
        values = adata.raw.obs_vector(value_to_plot)
    else:
        values = adata.obs_vector(value_to_plot, layer=layer)
    if mask is not None:
        values[~mask] = np.nan
    if groups and isinstance(values, pd.Categorical):
        values = values.remove_categories(values.categories.difference(groups))
    return values

