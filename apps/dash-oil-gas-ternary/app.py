import pathlib
import os

import pandas as pd
import numpy as np

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State

import constants

# app initialize
app = dash.Dash(
    __name__,
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, initial-scale=1.0"}
    ],
)
server = app.server
app.config["suppress_callback_exceptions"] = True

# mapbox
mapbox_access_token = "pk.eyJ1IjoicGxvdGx5bWFwYm94IiwiYSI6ImNrOWJqb2F4djBnMjEzbG50amg0dnJieG4ifQ.Zme1-Uzoi75IaFbieBDl3A"

# Load data
APP_PATH = str(pathlib.Path(__file__).parent.resolve())

df = pd.read_csv(os.path.join(APP_PATH, os.path.join("data", "test_composition.csv")))
df_prod = pd.read_csv(
    os.path.join(APP_PATH, os.path.join("data", "YearlyProduction_table_1.csv"))
)

# Assign color to legend
colormap = {}
for ind, formation_name in enumerate(df["fm_name"].unique().tolist()):
    colormap[formation_name] = constants.colors[ind]


def build_banner():
    return html.Div(
        id="banner",
        className="banner",
        children=[
            html.Img(src=app.get_asset_url("dash-logo.png")),
            html.H6("Oil and gas ternary map"),
        ],
    )


def build_graph_title(title):
    return html.P(className="graph-title", children=title)


def generate_production_plot(processed_data):
    """
    :param processed_data: List containing two lists, one containing well ID information, and the second containing
    rock formation type associated with the well
    :return: Figure object
    """
    layout = dict(
        xaxis=dict(title="Year"), yaxis=dict(title="GAS Production (mcf)", type="log")
    )

    data = []
    for well_id, formation in list(
        zip(processed_data["well_id"], processed_data["formation"])
    ):
        well_prod = df_prod[df_prod["RecordNumber"] == well_id]
        new_trace = dict(
            x=well_prod["Year"],
            y=well_prod["VolumeMCF"],
            name=str(well_id),
            mode="lines+markers",
            hoverinfo="x+y+name",
            marker=dict(
                symbol="hexagram-open", line={"width": "0.5"}, color=colormap[formation]
            ),
            line=dict(shape="spline"),
            showlegend=True,
        )
        data.append(new_trace)
    return {"data": data, "layout": layout}


def generate_well_map(dff, selected_data, style):
    """
    Generate well map based on selected data.

    :param dff: dataframe for generate plot.
    :param selected_data: Processed dictionary for plot generation with defined selected points.
    :param style: mapbox visual style.
    :return: Plotly figure object.
    """

    layout = go.Layout(
        clickmode="event+select",
        dragmode="lasso",
        showlegend=True,
        autosize=True,
        hovermode="closest",
        margin=dict(l=0, r=0, t=0, b=0),
        mapbox=go.layout.Mapbox(
            accesstoken=mapbox_access_token,
            bearing=0,
            center=go.layout.mapbox.Center(lat=37.497562, lon=-82.755728),
            pitch=0,
            zoom=8,
            style=style,
        ),
        legend=dict(
            bgcolor="#1f2c56",
            orientation="h",
            font=dict(color="white"),
            x=0,
            y=0,
            yanchor="bottom",
        ),
    )

    formations = dff["fm_name"].unique().tolist()

    data = []

    for formation in formations:
        selected_index = None
        if formation in selected_data:
            selected_index = selected_data[formation]

        text_list = list(
            map(
                lambda item: "Well ID:" + str(int(item)),
                dff[dff["fm_name"] == formation]["RecordNumber"],
            )
        )
        op_list = dff[dff["fm_name"] == formation]["op"].tolist()

        text_list = [op_list[i] + "<br>" + text_list[i] for i in range(len(text_list))]

        new_trace = go.Scattermapbox(
            lat=dff[dff["fm_name"] == formation]["nlat83"],
            lon=dff[dff["fm_name"] == formation]["wlon83"],
            mode="markers",
            marker={"color": colormap[formation], "size": 9},
            text=text_list,
            name=formation,
            selectedpoints=selected_index,
            customdata=dff[dff["fm_name"] == formation]["RecordNumber"],
        )
        data.append(new_trace)

    return {"data": data, "layout": layout}


def generate_ternary_map(dff, selected_data, contour_visible, marker_visible):
    """
    Generate ternary plot based on selected data.

    :param dff: dataframe for generate plot.
    :param selected_data: Processed dictionary for plot generation with defined selected points.
    :param contour_visible: Contour trace visibility.
    :param marker_visible: Marker trace visibility.
    :return: ternary map figure object.
    """

    # Generate contour

    contour_traces = []
    for ind, key in enumerate(constants.ternary_contour.keys()):
        trace = dict(
            name=key,
            type="scatterternary",
            a=[k["Quartz"] for k in constants.ternary_contour[key]],
            b=[k["Carbonate"] for k in constants.ternary_contour[key]],
            c=[k["Clay"] for k in constants.ternary_contour[key]],
            mode="lines",
            line=dict(color="#444", width=0.5),
            fill="toself",
            fillcolor=constants.ternary_color[ind],
            opacity=0.8,
            hoverinfo="none",
            showlegend=False,
            visible=contour_visible,
        )
        contour_traces.append(trace)

    contour_text = generate_contour_text_layer(contour_visible)

    # Layout
    layout = {
        "dragmode": "lasso",
        "ternary": {
            "sum": 100,
            "aaxis": {
                "title": {
                    "text": "Quartz",
                    "font": {"family": "Open Sans", "size": 15, "color": "white"},
                },
                "min": -2,
                "linewidth": 1.5,
                "ticks": "outside",
            },
            "baxis": {
                "title": {
                    "text": "Carbonate",
                    "font": {"family": "Open Sans", "size": 15, "color": "white"},
                },
                "min": -2,
                "linewidth": 1.5,
                "ticks": "outside",
            },
            "caxis": {
                "title": {
                    "text": "Clay",
                    "font": {"family": "Open Sans", "size": 15, "color": "white"},
                },
                "min": -2,
                "linewidth": 1.5,
                "ticks": "outside",
            },
        },
        "margin": dict(l=110, r=50, t=50, b=50),
        "paper_bgcolor": "#192444",
        "plot_bgcolor": "#192444",
        "showLegend": False,
        "font": {"color": "white"},
        "annotations": {"visible": False},
        "autosize": True,
    }

    hovertemplate = "<b> %{text}</b><br><br> Quartz: %{a:.0f}<br>Carbonate : %{b:.0f}<br> Clay: %{c:.0f}<extra></extra>"

    formations = dff["fm_name"].unique().tolist()

    data_traces = []
    for key in formations:
        if selected_data:
            select_indices = selected_data[key]
        else:
            select_indices = None

        new_data_trace = dict(
            text=list(
                map(
                    lambda item: "Well ID:" + str(int(item)),
                    dff[dff["fm_name"] == key]["RecordNumber"],
                )
            ),
            name=key,
            customdata=dff[dff["fm_name"] == key]["RecordNumber"],
            type="scatterternary",
            a=dff[dff["fm_name"] == key]["Quartz"],
            b=dff[dff["fm_name"] == key]["Carbonate"],
            c=dff[dff["fm_name"] == key]["Clay"],
            mode="markers",
            hovertemplate=hovertemplate,
            showlegend=False,
            marker={
                "color": colormap[key],
                "size": 8,
                "line": {"color": "#000000", "width": 0.2},
            },
            selectedpoints=select_indices,
            visible=marker_visible,
        )
        data_traces.append(new_data_trace)

    return {"data": contour_traces + contour_text + data_traces, "layout": layout}


def generate_contour_text_layer(contour_visible):
    layer = []
    for key, value in constants.ternary_contour.items():
        a = np.mean([i["Quartz"] for i in value])
        b = np.mean([i["Carbonate"] for i in value])
        c = np.mean([i["Clay"] for i in value])

        key_br = key.replace(" ", "<br>")

        new_trace = generate_contour_text(a, b, c, key, key_br, contour_visible)
        layer.append(new_trace)

    return layer


def generate_contour_text(a, b, c, name, text, visible):
    return dict(
        type="scatterternary",
        a=[a],
        b=[b],
        c=[c],
        name=name,
        text=text,
        mode="text",
        hoverinfo="none",
        textposition="middle center",
        textfont={"size": 11, "color": "#000000", "family": "sans-serif"},
        showlegend=False,
        legendgroup="Rock type",
        visible=visible,
    )


def generate_formation_bar(dff, selected_data):
    """
    Generate bar plot based on selected data.

        :param dff: dataframe for generate plot.
        :param selected_data: Processed dictionary for plot generation with defined selected points.
        :return: ternary map figure object.

    """

    layout = go.Layout(
        showlegend=False,
        hovermode="closest",
        xaxis=dict(tickangle=-45, title="Formations"),
        yaxis=dict(title="Well Counts"),
        clickmode="event+select",
    )

    formations = dff["fm_name"].unique().tolist()

    if selected_data:
        data = []
        for i in formations:
            selected_points = []
            select_indices = selected_data[i]
            if select_indices is not None and len(select_indices) > 0:
                selected_points = [0]
            new_trace = go.Bar(
                x=[i],
                y=[len(dff[dff["fm_name"] == i])],
                name=i,
                hoverinfo="x+y",
                marker={"color": colormap[i]},
                selectedpoints=selected_points,
            )
            data.append(new_trace)

    else:
        data = []
        for i in formations:
            new_trace = go.Bar(
                x=[i],
                y=[len(dff[dff["fm_name"] == i])],
                name=i,
                marker={"color": colormap[i]},
                selectedpoints=None,
            )
            data.append(new_trace)

    return {"data": data, "layout": layout}


# Helper for extracting select index from mapbox and tern selectData
def get_selection(data, formation, selection_data, starting_index):
    ind = []
    current_curve = data["fm_name"].unique().tolist().index(formation)
    for point in selection_data["points"]:
        if point["curveNumber"] - starting_index == current_curve:
            ind.append(point["pointNumber"])
    return ind


# Helper for extracting select index from bar
def get_selection_by_bar(bar_selected_data):
    dict = {}
    if bar_selected_data is not None:
        for point in bar_selected_data["points"]:
            if point["x"] is not None:
                dict[(point["x"])] = list(range(0, point["y"]))
    return dict


app.layout = html.Div(
    children=[
        html.Div(
            id="top-row",
            children=[
                html.Div(
                    className="row",
                    id="top-row-header",
                    children=[
                        html.Div(
                            id="header-container",
                            children=[
                                build_banner(),
                                html.P(
                                    id="instructions",
                                    children="Select data points from the well map, ternary map or bar graph to "
                                    "visualize cross-filtering to other plots. Selection could be done by "
                                    "clicking on individual data points or using the lasso tool to capture "
                                    "multiple data points or bars. With the box tool from modebar, multiple "
                                    "regions can be selected by holding the SHIFT key while clicking and "
                                    "dragging.",
                                ),
                                build_graph_title("Select Operator"),
                                dcc.Dropdown(
                                    id="operator-select",
                                    options=[
                                        {"label": i, "value": i}
                                        for i in df["op"].unique().tolist()
                                    ],
                                    multi=True,
                                    value=[
                                        df["op"].unique().tolist()[0],
                                        df["op"].unique().tolist()[1],
                                    ],
                                ),
                            ],
                        )
                    ],
                ),
                html.Div(
                    className="row",
                    id="top-row-graphs",
                    children=[
                        # Well map
                        html.Div(
                            id="well-map-container",
                            children=[
                                build_graph_title("Well Map"),
                                dcc.RadioItems(
                                    id="mapbox-view-selector",
                                    options=[
                                        {"label": "basic", "value": "basic"},
                                        {"label": "satellite", "value": "satellite"},
                                        {"label": "outdoors", "value": "outdoors"},
                                        {
                                            "label": "satellite-street",
                                            "value": "mapbox://styles/mapbox/satellite-streets-v9",
                                        },
                                    ],
                                    value="basic",
                                ),
                                dcc.Graph(
                                    id="well-map",
                                    figure={
                                        "layout": {
                                            "paper_bgcolor": "#192444",
                                            "plot_bgcolor": "#192444",
                                        }
                                    },
                                    config={"scrollZoom": True, "displayModeBar": True},
                                ),
                            ],
                        ),
                        # Ternary map
                        html.Div(
                            id="ternary-map-container",
                            children=[
                                html.Div(
                                    id="ternary-header",
                                    children=[
                                        build_graph_title(
                                            "Shale Mineralogy Composition"
                                        ),
                                        dcc.Checklist(
                                            id="ternary-layer-select",
                                            options=[
                                                {
                                                    "label": "Well Data",
                                                    "value": "Well Data",
                                                },
                                                {
                                                    "label": "Rock Type",
                                                    "value": "Rock Type",
                                                },
                                            ],
                                            value=["Well Data", "Rock Type"],
                                        ),
                                    ],
                                ),
                                dcc.Graph(
                                    id="ternary-map",
                                    figure={
                                        "layout": {
                                            "paper_bgcolor": "#192444",
                                            "plot_bgcolor": "#192444",
                                        }
                                    },
                                    config={
                                        "scrollZoom": True,
                                        "displayModeBar": False,
                                    },
                                ),
                            ],
                        ),
                    ],
                ),
            ],
        ),
        html.Div(
            className="row",
            id="bottom-row",
            children=[
                # Formation bar plots
                html.Div(
                    id="form-bar-container",
                    className="six columns",
                    children=[
                        build_graph_title("Well count by formations"),
                        dcc.Graph(id="form-by-bar"),
                    ],
                ),
                html.Div(
                    # Selected well productions
                    id="well-production-container",
                    className="six columns",
                    children=[
                        build_graph_title("Individual well annual production"),
                        dcc.Graph(id="production-fig"),
                    ],
                ),
            ],
        ),
    ]
)


# Update bar plot
@app.callback(
    Output("form-by-bar", "figure"),
    [
        Input("well-map", "selectedData"),
        Input("ternary-map", "selectedData"),
        Input("operator-select", "value"),
    ],
)
def update_bar(map_selected_data, tern_selected_data, op_select):
    dff = df[df["op"].isin(op_select)]

    formations = dff["fm_name"].unique().tolist()
    # Find which one has been triggered
    ctx = dash.callback_context

    prop_id = ""
    prop_type = ""
    if ctx.triggered:
        splitted = ctx.triggered[0]["prop_id"].split(".")
        prop_id = splitted[0]
        prop_type = splitted[1]

    processed_data = {}
    if prop_id == "well-map" and prop_type == "selectedData":
        for formation in formations:
            if map_selected_data is None:
                processed_data[formation] = [
                    0
                ]  # [0] is the default value to select current bar
            else:
                processed_data[formation] = get_selection(
                    dff, formation, map_selected_data, 0
                )

    elif prop_id == "ternary-map" and prop_type == "selectedData":

        for formation in formations:
            if tern_selected_data is None:
                processed_data[formation] = [0]
            else:
                processed_data[formation] = get_selection(
                    dff, formation, tern_selected_data, 32
                )

    else:
        for formation in formations:
            processed_data[formation] = [0]

    return generate_formation_bar(dff, processed_data)


# Update ternary map
@app.callback(
    Output("ternary-map", "figure"),
    [
        Input("well-map", "selectedData"),
        Input("form-by-bar", "selectedData"),
        Input("form-by-bar", "clickData"),
        Input("operator-select", "value"),
        Input("ternary-layer-select", "value"),
    ],
    state=[State("ternary-map", "figure")],
)
def update_ternary_map(
    map_selected_data,
    bar_selected_data,
    bar_click_data,
    op_select,
    layer_select,
    curr_fig,
):
    marker_visible = contour_visible = True

    dff = df[df["op"].isin(op_select)]
    formations = dff["fm_name"].unique().tolist()

    # Find which one has been triggered
    ctx = dash.callback_context

    if ctx.triggered:
        splitted = ctx.triggered[0]["prop_id"].split(".")
        prop_id = splitted[0]
        prop_type = splitted[1]
    else:
        return curr_fig

    processed_data = {}

    if prop_id != "ternary-layer-select":
        if prop_id == "well-map" and prop_type == "selectedData":

            for formation in formations:
                if map_selected_data is None:
                    processed_data[formation] = None
                else:
                    processed_data[formation] = get_selection(
                        dff, formation, map_selected_data, 0
                    )

        elif prop_id == "form-by-bar" and prop_type == "selectedData":

            processed_data = get_selection_by_bar(bar_selected_data)

            for formation in formations:
                if bar_selected_data is None:
                    processed_data[formation] = None
                elif formation not in processed_data:
                    processed_data[formation] = []

        elif prop_id == "form-by-bar" and prop_type == "clickData":

            processed_data = get_selection_by_bar(bar_click_data)
            for formation in formations:
                if bar_click_data is None:
                    processed_data[formation] = None
                elif formation not in processed_data:
                    processed_data[formation] = []

        else:

            for formation in formations:
                processed_data[formation] = None

        return generate_ternary_map(
            dff, processed_data, contour_visible, marker_visible
        )

    if prop_id == "ternary-layer-select":
        if curr_fig is not None:
            if "Well Data" not in layer_select:
                marker_visible = "legendonly"
            if "Rock Type" not in layer_select:
                contour_visible = "legendonly"

            for contour_dict in curr_fig["data"][:32]:
                contour_dict["visible"] = contour_visible

            for marker_dict in curr_fig["data"][32:]:
                marker_dict["visible"] = marker_visible
            return curr_fig
        else:
            return curr_fig


# Update well map
@app.callback(
    Output("well-map", "figure"),
    [
        Input("ternary-map", "selectedData"),
        Input("form-by-bar", "selectedData"),
        Input("form-by-bar", "clickData"),
        Input("operator-select", "value"),
        Input("mapbox-view-selector", "value"),
    ],
)
def update_well_map(
    tern_selected_data, bar_selected_data, bar_click_data, op_select, mapbox_view
):
    dff = df[df["op"].isin(op_select)]
    formations = dff["fm_name"].unique().tolist()

    # Find which one has been triggered
    ctx = dash.callback_context

    prop_id = ""
    prop_type = ""
    if ctx.triggered:
        splitted = ctx.triggered[0]["prop_id"].split(".")
        prop_id = splitted[0]
        prop_type = splitted[1]

    processed_data = {}

    if prop_id == "ternary-map":
        for formation in formations:
            if tern_selected_data is None:
                processed_data[formation] = None
            else:
                processed_data[formation] = get_selection(
                    dff, formation, tern_selected_data, 32
                )

    elif prop_id == "form-by-bar":

        bar_data = ""
        if prop_type == "selectedData":
            bar_data = bar_selected_data
        elif prop_type == "clickData":
            bar_data = bar_click_data

        processed_data = get_selection_by_bar(bar_data)

        for formation in formations:
            if bar_data is None:
                processed_data[formation] = None
            elif formation not in processed_data:
                processed_data[formation] = []

    else:
        for formation in formations:
            processed_data[formation] = None

    return generate_well_map(dff, processed_data, mapbox_view)


# Update production plot
@app.callback(
    Output("production-fig", "figure"),
    [
        Input("well-map", "selectedData"),
        Input("ternary-map", "selectedData"),
        Input("form-by-bar", "selectedData"),
        Input("operator-select", "value"),
    ],
)
def update_production(map_select, tern_select, bar_select, op_select):
    dff = df[df["op"].isin(op_select)]

    # Find which one has been triggered
    ctx = dash.callback_context

    prop_id = ""
    prop_type = ""
    if ctx.triggered:
        splitted = ctx.triggered[0]["prop_id"].split(".")
        prop_id = splitted[0]
        prop_type = splitted[1]

    processed_data_init = {}
    processed_data_init["well_id"] = dff["RecordNumber"].tolist()
    processed_data_init["formation"] = dff["fm_name"].tolist()

    if prop_id == "well-map" and prop_type == "selectedData":
        if map_select is not None:
            processed_data = {"well_id": [], "formation": []}
            for point in map_select["points"]:
                processed_data["well_id"].append(point["customdata"])
                processed_data["formation"].append(
                    dff[dff["RecordNumber"] == point["customdata"]]["fm_name"].tolist()[
                        0
                    ]
                )
        else:
            processed_data = processed_data_init

    elif prop_id == "ternary-map" and prop_type == "selectedData":
        if tern_select is not None:
            processed_data = {"well_id": [], "formation": []}
            for point in tern_select["points"]:
                if "customdata" in point:
                    processed_data["well_id"].append(point["customdata"])
                    processed_data["formation"].append(
                        dff[dff["RecordNumber"] == point["customdata"]][
                            "fm_name"
                        ].tolist()[0]
                    )

        else:
            processed_data = processed_data_init

    elif prop_id == "form-by-bar" and prop_type == "selectedData":
        if bar_select is not None:
            processed_data = {"well_id": [], "formation": []}

            # Find all wells according to selected formation category
            for point in bar_select["points"]:
                selected_form = point["x"]
                selected_well = dff[dff["fm_name"] == point["x"]][
                    "RecordNumber"
                ].tolist()
                for well in selected_well:
                    processed_data["well_id"].append(int(well))
                    processed_data["formation"].append(selected_form)

        else:
            processed_data = processed_data_init
    else:
        processed_data = processed_data_init

    return generate_production_plot(processed_data)


# Running the server
if __name__ == "__main__":
    app.run_server(debug=True)
