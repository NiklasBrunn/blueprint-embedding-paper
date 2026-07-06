"""Explore tab: layout + callbacks for the client-side recolor UMAP.

The plot is a plain ``<div id="exp-graph">`` driven by ``assets/explore.js``
(window.initExplore / applyGeneVector, exposed as the ``explore`` clientside
namespace). The server only ships the bundle once and answers gene-vector
requests; recolor/temporal scrubbing happen entirely in-browser.
"""
from __future__ import annotations

import dash
import dash_bootstrap_components as dbc
from dash import dcc, html, Input, Output, State

from webapp import coloring, explore_bundle


def panel() -> html.Div:
    # --- left column: data + blueprint + downloads (reused markup) --------
    data_card = dbc.Card([
        dbc.CardHeader("Data"),
        dbc.CardBody([
            dcc.Upload(id="exp-upload-h5ad", accept=".h5ad,.h5ad.gz",
                       children=html.Div(["Drag & drop or ",
                                          html.A("select an h5ad")]),
                       className="ukf-drop"),
            html.Div(id="exp-h5ad-name", className="fst-italic mt-1"),
            html.Div("— or —", className="text-muted text-center my-1 small"),
            dbc.InputGroup([
                dbc.Button("Browse…", id="exp-h5ad-browse", n_clicks=0,
                           color="primary", outline=True),
                dbc.Input(id="exp-h5ad-path", type="text", debounce=True,
                          placeholder="select or paste a local .h5ad path"),
            ]),
            dcc.Loading(html.Div(id="exp-h5ad-info",
                                 className="small mt-1",
                                 style={"minHeight": "20px"}),
                        type="dot", color="#344A9A"),
            html.Div("Selecting only picks the file — click Map to load and "
                     "embed it (a progress bar shows while it runs). A local "
                     "file is read straight from disk, so there is no size limit.",
                     className="text-muted small mt-1"),
        ]),
    ])
    bp_card = dbc.Card([
        dbc.CardHeader("Blueprint"),
        dbc.CardBody([
            html.Div(id="exp-bp-ready"),
            dbc.RadioItems(id="exp-bp-source", value="session",
                           options=[{"label": "Use in-session blueprint",
                                     "value": "session"},
                                    {"label": "Use uploaded blueprint",
                                     "value": "upload"}]),
            html.Div(id="exp-bp-upload-wrap", style={"display": "none"},
                     children=[
                         dcc.Upload(id="exp-upload-bp", accept=".py",
                                    children=html.Div(
                                        ["Drag & drop or ",
                                         html.A("select a blueprint .py")]),
                                    className="ukf-drop"),
                         html.Div(id="exp-bp-name", className="fst-italic mt-1"),
                         html.Div("— or —",
                                  className="text-muted text-center my-1 small"),
                         dbc.InputGroup([
                             dbc.Button("Browse…", id="exp-bp-browse",
                                        n_clicks=0, color="primary",
                                        outline=True),
                             dbc.Input(id="exp-bp-path", type="text",
                                       debounce=True,
                                       placeholder="select or paste a local "
                                                   "blueprint .py path"),
                         ]),
                         html.Div(id="exp-bp-info", className="small mt-1",
                                  style={"minHeight": "20px"}),
                     ]),
            dbc.Button("Map", id="exp-map-btn", n_clicks=0,
                       color="success", class_name="mt-3 w-100"),
            dbc.Button("Clear file paths", id="exp-clear-paths", n_clicks=0,
                       color="secondary", outline=True, size="sm",
                       class_name="mt-2 w-100"),
            html.Div("Clears the typed/Browsed paths so a dragged-in file is "
                     "used instead (a path otherwise takes precedence).",
                     className="text-muted small mt-1"),
            dcc.Loading(
                html.Div(id="exp-coverage", className="mt-2",
                         style={"minHeight": "24px"}),
                custom_spinner=html.Div(className="mt-2", children=[
                    html.Div("Mapping… running the embedder + UMAP "
                             "(can take a while for large data)",
                             className="text-muted small mb-1"),
                    dbc.Progress(value=100, striped=True, animated=True,
                                 color="success", style={"height": "12px"}),
                ]),
            ),
        ]),
    ])
    downloads = dbc.Card(class_name="mt-3", children=[
        dbc.CardHeader("Downloads"),
        dbc.CardBody([
            dbc.Button("Download augmented h5ad", id="exp-dl-h5ad-btn",
                       n_clicks=0, color="primary", size="sm"),
            dcc.Download(id="exp-dl-h5ad"),
            html.Div(className="mt-2", children=[
                dbc.RadioItems(id="exp-fig-fmt", inline=True, value="png",
                               options=[{"label": "PNG", "value": "png"},
                                        {"label": "SVG", "value": "svg"}]),
                dbc.Button("Download plot", id="exp-dl-fig-btn", n_clicks=0,
                           color="primary", size="sm", class_name="mt-1"),
            ]),
            html.Div(className="mt-2", children=[
                dbc.Label("Tree mode", class_name="mb-1 small text-muted"),
                dbc.RadioItems(id="exp-tree-kind", inline=True, value="celltype",
                               options=[{"label": "Cell-type", "value": "celltype"},
                                        {"label": "Time", "value": "time"}]),
                dbc.RadioItems(id="exp-tree-fmt", inline=True, value="png",
                               class_name="mt-1",
                               options=[{"label": "PNG", "value": "png"},
                                        {"label": "PDF", "value": "pdf"},
                                        {"label": "SVG", "value": "svg"}]),
                dbc.Button("Download tree plot", id="exp-dl-tree-btn", n_clicks=0,
                           color="primary", size="sm", class_name="mt-1"),
                dcc.Download(id="exp-dl-tree"),
            ]),
        ]),
    ])

    # --- right column: the two-pane client-side explore area --------------
    controls = html.Div(className="controls", children=[
        html.Label(["Color by: ", html.Select(id="exp-mode")]),
        html.Div(id="exp-branch-group", className="control-group",
                 children=[html.Label(["Branch: ", html.Select(id="exp-branch")])]),
        html.Div(id="exp-palette-group", className="control-group", children=[
            html.Label(["Colors: ", html.Select(id="exp-palette", children=[
                html.Option("Viridis (unsigned)", value="viridis"),
                html.Option("Diverging (signed)", value="divbo"),
            ])])]),
        # The scrubber is a plain native <input type="range"> created by
        # explore.js inside this host — NOT a React-controlled dcc.Input, whose
        # static min/max props would fight the per-branch range we set in JS.
        html.Div(id="exp-time-group", className="control-group", children=[
            html.Label("t: "), html.Span(id="exp-time-host")]),
        html.Span(id="exp-stateline"),
    ])
    explore_area = html.Div(className="exp-app", children=[
        html.Div(id="exp-sidebar"),
        html.Div(className="exp-main", children=[
            html.Div("No mapping yet — pick a dataset + blueprint and click Map.",
                     id="exp-current", className="exp-current"),
            controls,
            html.Div(id="exp-ticks"),
            html.Div(id="exp-cattitle"),
            html.Div(id="exp-catlegend"),
            html.Div(id="exp-graph"),
            html.P("Drag the time slider to scrub; the sidebar snaps to the "
                   "dominant axis.", id="exp-footnote", className="note"),
        ]),
    ])
    canvas = dbc.Card([dbc.CardBody([explore_area])])

    return html.Div([
        dbc.Row([
            dbc.Col([data_card, html.Div(className="mt-3"), bp_card, downloads],
                    md=4),
            dbc.Col([canvas], md=8),
        ]),
        dcc.Store(id="explore-bundle"),     # client JS reads this
        dcc.Store(id="gene-vector"),        # server -> JS gene vector
        dcc.Store(id="exp-paths"),          # server-side paths for dl + gene fetch
        dcc.Store(id="exp-init-dummy"),     # clientside init callback output sink
        dcc.Store(id="exp-gene-dummy"),     # clientside applyGene output sink
        dcc.Store(id="exp-fig-dummy"),      # clientside figure-download output sink
        # JS requests a gene by setting this input's value (via a native input
        # event so the change flows through dcc.Input's normal callback path —
        # clientside set_props updates props WITHOUT triggering callbacks).
        # Kept in the DOM but visually hidden (type="text" so React binds onChange).
        dcc.Input(id="exp-gene-select", type="text", value="",
                  debounce=False, style={"display": "none"}),
    ])


def gene_vector(h5ad_path, gene, cell_ids=None):
    """log1p of one gene (case-insensitive) from an h5ad, aligned to the plot.

    ``cell_ids`` is the displayed cell order (from
    ``explore_bundle.displayed_cell_ids``); the returned values are reindexed to
    that order so they line up one-to-one with the bundle's coords (h5ad obs
    order may differ from the embedding order, and the plot may be downsampled).
    Cells absent from the h5ad get 0.0.

    Returns ``{"gene": gene, "values": [float, ...]}`` or ``None`` when inputs
    are missing or the gene isn't present."""
    if not h5ad_path or not gene:
        return None
    try:
        import anndata as ad
        import numpy as np
        import pandas as pd

        a = ad.read_h5ad(h5ad_path)
        X = a.X
        var_names = list(map(str, a.var_names))
        lower = [v.lower() for v in var_names]
        key = str(gene).lower()
        if key not in lower:
            return None
        col = X[:, lower.index(key)]
        col = np.asarray(col.todense()).ravel() if hasattr(col, "todense") \
            else np.asarray(col).ravel()
        vals = np.log1p(col)
        if cell_ids is not None:
            series = pd.Series(vals, index=a.obs_names.astype(str))
            vals = series.reindex([str(c) for c in cell_ids]).fillna(0.0).to_numpy()
    except Exception:
        return None
    return {"gene": gene, "values": [float(v) for v in np.asarray(vals).ravel()]}


def register(app) -> None:
    """Wire the explore tab's clientside + server callbacks."""
    # 1. Ship the bundle to the browser once it's mapped.
    app.clientside_callback(
        dash.dependencies.ClientsideFunction(namespace="explore",
                                             function_name="init"),
        Output("exp-init-dummy", "data"),
        Input("explore-bundle", "data"),
        prevent_initial_call=True,
    )

    # 2. Apply a fetched gene vector in-browser.
    app.clientside_callback(
        dash.dependencies.ClientsideFunction(namespace="explore",
                                             function_name="applyGene"),
        Output("exp-gene-dummy", "data"),
        Input("gene-vector", "data"),
        prevent_initial_call=True,
    )

    # 3. Server gene fetch: JS sets exp-gene-select, we return the vector.
    @app.callback(
        Output("gene-vector", "data"),
        Input("exp-gene-select", "value"),
        State("exp-paths", "data"),
        prevent_initial_call=True,
    )
    def _gene(gene, paths):
        if not gene or not paths:
            return None
        cell_ids = None
        parquet = paths.get("parquet")
        if parquet:
            try:
                cell_ids = explore_bundle.displayed_cell_ids(parquet)
            except Exception:
                cell_ids = None
        return gene_vector(paths.get("h5ad"), gene, cell_ids)

    # 4. Client-side figure export (PNG/SVG via Plotly.downloadImage).
    app.clientside_callback(
        """
        function(n, fmt){
          if(!n) return window.dash_clientside.no_update;
          var gd = document.getElementById("exp-graph");
          if(gd && window.Plotly){
            window.Plotly.downloadImage(gd, {format: fmt||"png", filename: "umap",
              width: (gd.clientWidth||1000), height: (gd.clientHeight||700)});
          }
          return window.dash_clientside.no_update;
        }
        """,
        Output("exp-fig-dummy", "data"),
        Input("exp-dl-fig-btn", "n_clicks"),
        State("exp-fig-fmt", "value"),
        prevent_initial_call=True,
    )

    # 5. Clear the typed/Browsed file paths so a dragged-in file wins (a path
    #    otherwise takes precedence). allow_duplicate: the Browse callbacks in
    #    server.py also output these path values.
    @app.callback(
        Output("exp-h5ad-path", "value", allow_duplicate=True),
        Output("exp-bp-path", "value", allow_duplicate=True),
        Input("exp-clear-paths", "n_clicks"),
        prevent_initial_call=True,
    )
    def _clear_paths(n):
        return "", ""
