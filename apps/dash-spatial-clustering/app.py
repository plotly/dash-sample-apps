import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.graph_objects as go
import plotly.express as px

import pandas as pd

from helpers import (
    austin_listings,
    zc_link,
    apply_clustering,
    rating_clustering,
    review_columns,
)

app = dash.Dash(__name__)
app.title = "Real Estate Spatial Clustering"
server = app.server
app.config.suppress_callback_exceptions = True

# CONSTANTS
grouped = austin_listings.groupby("zipcode").size()

mapbox_token = "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNrOWJqb2F4djBnMjEzbG50amg0dnJieG4ifQ.Zme1-Uzoi75IaFbieBDl3A"

geo_colors = [
    "#8dd3c7",
    "#ffd15f",
    "#bebada",
    "#fb8072",
    "#80b1d3",
    "#fdb462",
    "#b3de69",
    "#fccde5",
    "#d9d9d9",
    "#bc80bd",
    "#ccebc5",
]

bar_coloway = [
    "#fa4f56",
    "#8dd3c7",
    "#ffffb3",
    "#bebada",
    "#80b1d3",
    "#fdb462",
    "#b3de69",
    "#fccde5",
    "#d9d9d9",
    "#bc80bd",
    "#ccebc5",
    "#ffed6f",
]

intro_text = """
**About this app**

This app applies spatial clustering and regionalization analysis to discover the [dataset of AirBnb listings in the
city of Austin](http://insideairbnb.com/get-the-data.html). Models are created using [pysal](https://pysal.readthedocs.io/en/latest/)
and scikit-learn.

Select the type of model from radioitem, click on the button to run clustering and visualize output regions geographically on the map, computing may take seconds to finish. Click
on regions on the map to update the number of airbnb listings from your highlighted group.

"""

listing_txt = """
out of *{}* total listings
""".format(
    grouped.sum()
)

INIT_THRESHOLD_VAL = 5


def header_section():
    return html.Div(
        [
            html.Img(src=app.get_asset_url("dash-logo.png"), className="logo"),
            html.H4("Spatial Clustering"),
        ],
        className="header__title",
    )


def make_base_map():
    # Scattermapbox with geojson layer, plot all listings on mapbox
    customdata = list(
        zip(
            austin_listings["host_name"],
            austin_listings["name"],
            austin_listings["host_since"],
            austin_listings["price"],
            austin_listings["accommodates"],
            austin_listings["availability_365"],
            round(austin_listings["availability_365"] / 365 * 100, 1),
        )
    )
    mapbox_figure = dict(
        type="scattermapbox",
        lat=austin_listings["latitude"],
        lon=austin_listings["longitude"],
        marker=dict(size=7, opacity=0.7, color="#550100"),
        customdata=customdata,
        name="Listings",
        hovertemplate="<b>Host: %{customdata[0]}</b><br><br>"
        "<b>%{customdata[1]}</b><br>"
        "<b>Host Since: </b>%{customdata[2]}<br>"
        "<b>Price: </b>%{customdata[3]} / night<br>"
        "<b>Person to accommodate: </b>%{customdata[4]}<br>"
        "<b>Yearly Availability: </b>%{customdata[5]} days/year (%{customdata[6]} %)",
    )

    layout = dict(
        mapbox=dict(
            style="streets",
            uirevision=True,
            accesstoken=mapbox_token,
            zoom=9,
            center=dict(
                lon=austin_listings["longitude"].mean(),
                lat=austin_listings["latitude"].mean(),
            ),
        ),
        shapes=[
            {
                "type": "rect",
                "xref": "paper",
                "yref": "paper",
                "x0": 0,
                "y0": 0,
                "x1": 1,
                "y1": 1,
                "line": {"width": 1, "color": "#B0BEC5"},
            }
        ],
        margin=dict(l=10, t=10, b=10, r=10),
        height=900,
        showlegend=True,
        hovermode="closest",
    )

    figure = {"data": [mapbox_figure], "layout": layout}
    return figure


def make_map_with_clustering(sel_ind, c_type, stored_data):
    """
    Update layers on map from clustering regions.
    :param sel_ind: lasso-select index from map.selectedData.
    :param c_type: cluster type.
    :param stored_data: datastore from computing.
    :return: Plotly figure object.
    """
    # Group based on zipcode
    figure = make_base_map()

    # Decrease opacity of scatter
    figure["data"][0]["marker"].update(opacity=0.02)
    figure["layout"].update(
        dragmode="lasso"
    )  # clickmode doesn't work but drag will select scatters

    db = pd.DataFrame()
    if c_type == "ht-cluster":
        db = pd.read_json(stored_data["ht"]["data"])
    elif c_type == "rating-cluster":
        db, p_val = pd.read_json(stored_data["rt"]["data"]), stored_data["rt"]["p_val"]

    for ind, i in enumerate(db["cl"].unique()):
        # Choro cluster by zipcode, split into different colored choro layer after clustering or regionalization.
        figure["data"].append(
            dict(
                type="choroplethmapbox",
                showlegend=True,
                geojson=zc_link,
                locations=db[db["cl"] == i]["zipcode"],
                z=list(1 for _ in range(len(db[db["cl"] == i]["zipcode"]))),
                hoverinfo="location",
                name="Group {}".format(ind + 1),
                customdata=list(ind for _ in range(len(db[db["cl"] == i]["zipcode"]))),
                selected=dict(marker=dict(opacity=1)),
                unselected=dict(marker=dict(opacity=0.2)),
                selectedpoints="" if ind == sel_ind or sel_ind is None else [],
                marker=dict(opacity=0.8, line=dict(width=1)),
                colorscale=[[0, geo_colors[ind]], [1, geo_colors[ind]]],
                showscale=False,
            )
        )

    # move scatter trace at the end and disable its hover effect
    figure["data"].append(figure["data"].pop(0))
    figure["data"][-1].update(dict(hoverinfo="skip", hovertemplate=""))

    return figure


def make_original_property_graph():
    # Property type without any grouping
    types = pd.get_dummies(austin_listings["property_type"])
    prop_types = types.join(austin_listings["zipcode"]).groupby("zipcode").sum()
    prop_types_pct = (prop_types * 100.0).div(prop_types.sum(axis=1), axis=0)

    # Plot horizontal bars
    bar_fig = []

    for prop in list(prop_types_pct.columns):
        bar_fig.append(
            go.Bar(
                name=prop,
                orientation="h",
                y=[str(i) for i in prop_types_pct.index],
                x=prop_types_pct[prop],
            )
        )
    figure = go.Figure(
        data=bar_fig,
        layout=dict(
            barmode="stack",
            colorway=bar_coloway,
            yaxis_type="category",
            margin=dict(l=10, r=10, t=10, b=10),
            showlegend=False,
            height=1000,
            yaxis=dict(title="ZipCode"),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        ),
    )
    return figure


def make_property_graph(pct: pd.DataFrame) -> dict:
    # Demographic explore patterns clustering results
    bar_fig = []
    for prop in pct.columns:
        bar_fig.append(
            go.Bar(
                name=prop,
                x=["Group {}".format(int(i) + 1) for i in pct.index],
                y=pct[prop],
                marker=dict(opacity=0.8, line=dict(color="#ddd")),
                orientation="v",
            )
        )

    fig = go.Figure(
        data=bar_fig,
        layout=dict(
            barmode="stack",
            colorway=bar_coloway,
            margin=dict(l=10, r=10, t=10, b=10),
            legend=dict(traceorder="reversed", font=dict(size=9)),
            height=400,
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        ),
    )

    return fig


def make_review_chart(pct_d, original=False):
    fig = []
    for ind, rows in pct_d.iterrows():
        fig.append(
            go.Scatter(
                x=[x.split("_")[-1].capitalize() for x in review_columns],
                y=rows,
                name="Group {}".format(ind + 1) if not original else "All listings",
                mode="markers+lines",
                hovertemplate="<b> Review Score %{x}: </b> %{y}",
                marker=dict(
                    size=12,
                    opacity=0.8,
                    color=geo_colors[ind],
                    line=dict(width=1, color="#ffffff"),
                ),
            )
        )
    return {
        "data": fig,
        "layout": dict(
            yaxis=dict(title="Scores"),
            margin=dict(l=40, r=10, t=10, b=70),
            hovermode="closest",
        ),
    }


def populate_init_data():
    # Initialize data store with default computing
    ht_res = apply_clustering()
    rt_res = rating_clustering(INIT_THRESHOLD_VAL)
    return {
        "ht": {"data": ht_res[0].to_json(), "pct": ht_res[1].to_json()},
        "rt": {
            "data": rt_res[0].to_json(),
            "p_val": rt_res[1],
            "pct_d": rt_res[2].to_json(),
        },
    }


# Dash App Layout
app.layout = html.Div(
    children=[
        dcc.Store(id="cluster-data-store", data=populate_init_data()),
        # Banner
        header_section(),
        html.Div(
            [
                html.Div(
                    children=[
                        html.Div(id="intro-text", children=dcc.Markdown(intro_text)),
                        html.P("Austin Airbnb listings"),
                        html.Hr(),
                        dcc.Graph(id="map", config={"responsive": True}),
                    ],
                    className="eight columns named-card",
                ),
                # Categorical properties by cluster e.g (property type stacked bar)
                html.Div(
                    children=[
                        html.P("Map Options"),
                        html.Hr(),
                        dcc.RadioItems(
                            id="cluster-ctl",
                            options=[
                                {"label": "All Listings", "value": "no-cluster"},
                                {
                                    "label": "Cluster based on house type",
                                    "value": "ht-cluster",
                                },
                                {
                                    "label": "Regionalization based on Ratings",
                                    "value": "rating-cluster",
                                },
                            ],
                            value="no-cluster",
                        ),
                        html.Div(
                            id="threshold-div",
                            children=[
                                html.P("Minimum Threshold (% of reviews / region)"),
                                dcc.Slider(
                                    id="regionalization-threshold",
                                    min=5,
                                    max=30,
                                    marks={
                                        5: "5%",
                                        10: "10%",
                                        15: "15%",
                                        20: "20%",
                                        25: "25%",
                                        30: "30%",
                                    },
                                    value=INIT_THRESHOLD_VAL,
                                    tooltip={"placement": "bottom"},
                                ),
                            ],
                        ),
                        html.Div(
                            html.Button(
                                "Run Clustering & Update Map",
                                id="btn-updt-map",
                                title="Click to run spatial clustering, computing could take seconds to complete.",
                                n_clicks=0,
                            ),
                            className="btn-outer",
                        ),
                        html.Hr(),
                        html.Div(
                            id="listing-div",
                            children=[
                                html.H6("1000", id="listing-ind"),
                                dcc.Markdown(listing_txt, id="total-listings"),
                            ],
                        ),
                        html.P("Property Types"),
                        html.Hr(),
                        dcc.Graph(id="geodemo-chart"),
                        html.P("User satisfactions"),
                        html.Hr(),
                        dcc.Graph(id="user-satisfaction"),
                    ],
                    className="four columns named-card",
                ),
            ],
            className="twelve columns",
        ),
    ],
    className="container twelve columns",
)


# =====Callbacks=====
@app.callback(
    Output("map", "figure"),
    [Input("map", "clickData"), Input("cluster-data-store", "data")],
    [State("cluster-ctl", "value")],
)
def update_map(region_select, ds, clustering_type):
    # Update map based on selectedData and stored calculation
    ctx = dash.callback_context
    print("cluster type", clustering_type)
    if ds is not None:
        if "ht" in clustering_type or "rating" in clustering_type:
            if ctx.triggered and "data-store" in ctx.triggered[0]["prop_id"]:
                figure = make_map_with_clustering(None, clustering_type, ds)
                return figure
            if region_select is not None:
                # Empty select will return the entire non-selected figure
                if not len(region_select["points"]):
                    return make_map_with_clustering(None, clustering_type, ds)
                else:
                    for point in region_select["points"]:
                        if len(str(point["customdata"])) == 1:
                            print(
                                "selected region: Group {}".format(
                                    point["customdata"] + 1
                                )
                            )
                            return make_map_with_clustering(
                                point["customdata"], clustering_type, ds
                            )
            return make_map_with_clustering(None, clustering_type, ds)
        else:
            return make_base_map()
    return make_base_map()


@app.callback(
    Output("cluster-data-store", "data"),
    [Input("btn-updt-map", "n_clicks")],
    [
        State("cluster-ctl", "value"),
        State("cluster-data-store", "data"),
        State("regionalization-threshold", "value"),
    ],
)
def update_ds(n_clicks, clustering_type, cur_ds, thr):
    # Apply algorithm and only apply computing once upon button=click,save it for figure loading
    if n_clicks:
        if clustering_type == "ht-cluster":
            # Apply KMeans and update datastore
            ht_res = apply_clustering()
            cur_ds.update(ht={"data": ht_res[0].to_json(), "pct": ht_res[1].to_json()})
        elif clustering_type == "rating-cluster":
            rt_res = rating_clustering(thr)
            cur_ds.update(
                rt={
                    "data": rt_res[0].to_json(),
                    "p_val": rt_res[1],
                    "pct_d": rt_res[2].to_json(),
                }
            )
        else:
            return cur_ds
        return cur_ds
    return cur_ds


@app.callback(
    [Output("listing-ind", "children"), Output("listing-ind", "style")],
    [Input("map", "clickData"), Input("btn-updt-map", "n_clicks")],
    [State("cluster-ctl", "value"), State("cluster-data-store", "data")],
)
def update_indicator(map_select, n_click, cluster_type, ds):
    ctx = dash.callback_context
    if ctx.triggered and "clickData" in ctx.triggered[0]["prop_id"]:
        ht = pd.read_json(ds["ht"]["data"])
        rt = pd.read_json(ds["rt"]["data"])
        if cluster_type == "ht-cluster":
            dff = ht
        elif cluster_type == "rating-cluster":
            dff = rt
        if map_select is None:
            return str(len(austin_listings)), {"color": "#550100"}

        else:
            for point in map_select["points"]:
                if len(str(point["customdata"])) == 1:
                    print("selected region: Group {}".format(point["customdata"] + 1))
                    sel_ind = point["customdata"]
                    zips = dff[dff["cl"] == sel_ind]["zipcode"]
                    count = 0
                    for zip in zips:
                        count += grouped[str(zip)]
                    return str(count), {"color": geo_colors[sel_ind]}

    return str(len(austin_listings)), {"color": "#550100"}


@app.callback(
    Output("geodemo-chart", "figure"),
    [Input("cluster-data-store", "data")],
    [State("cluster-ctl", "value")],
)
def update_geodemo_chart(ds, clustering_type):
    # Update upon clustering updates ds
    if ds:
        if clustering_type == "ht-cluster":
            pct = pd.read_json(ds["ht"]["pct"]).drop(review_columns, axis=1)
            return make_property_graph(pct)
        if clustering_type == "rating-cluster":
            pct_d = pd.read_json(ds["rt"]["pct_d"]).drop(review_columns, axis=1)
            return make_property_graph(pct_d)
        elif clustering_type == "no-cluster":
            return make_original_property_graph()
    return make_original_property_graph()


@app.callback(
    Output("user-satisfaction", "figure"),
    [Input("cluster-data-store", "data")],
    [State("cluster-ctl", "value")],
)
def update_review_chart(ds, clustering_type):
    # y: Average score, x: category, color: group
    empty_fig = {
        "data": [],
        "layout": dict(
            yaxis=dict(title="Scores"),
            margin=dict(l=30, r=10, t=10, b=70),
            hovermode="closest",
        ),
    }

    if clustering_type == "rating-cluster":
        pct_d = pd.read_json(ds["rt"]["pct_d"])[review_columns]
        pct_d[review_columns[1:]] = pct_d[review_columns[1:]] * 10
        return make_review_chart(pct_d)
    elif clustering_type == "ht-cluster":
        ht = pd.read_json(ds["ht"]["pct"])[review_columns]
        ht[review_columns[1:]] = ht[review_columns[1:]] * 10
        return make_review_chart(ht)
    elif clustering_type == "no-cluster":
        # Avg review scores for all data
        dff = austin_listings[review_columns].mean(axis=0, skipna=True)
        dff_t = pd.DataFrame(dff).transpose()
        dff_t[review_columns[1:]] = dff_t[review_columns[1:]] * 10
        return make_review_chart(dff_t, original=True)
    return empty_fig


@app.callback(Output("threshold-div", "style"), [Input("cluster-ctl", "value")])
def update_div(clustering_type):
    if clustering_type == "rating-cluster":
        return {"display": "block"}
    return {"display": "none"}


if __name__ == "__main__":
    app.run_server(
        debug=True, port=8051, dev_tools_hot_reload=False, use_reloader=False
    )
