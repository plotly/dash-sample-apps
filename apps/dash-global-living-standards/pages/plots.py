import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from ipywidgets import widgets

##################################################################################################

# datasets needed for plots
data = pd.read_csv("df_final.csv")

city_options = data["City"].tolist()

app = dash.Dash()

app.layout = html.Div(
    [
        html.H2("Bubble"),
        html.Div(
            [
                dcc.Dropdown(
                    id="City",
                    options=[{"label": i, "value": i} for i in city_options],
                    value="All cities"
                    # multi=True
                )
            ],
            style={"width": "25%", "display": "inline-block"},
        ),
        # dcc.Graph(id='funnel-graph'),
        dcc.Graph(id="radar"),
        dcc.Graph(id="bubble"),
    ]
)

# plotting age groups for each country
# @app.callback(
# dash.dependencies.Output('funnel-graph', 'figure'),
# [dash.dependencies.Input('Country', 'value')])


def update_demo(Country):
    if Country == "All Countries":
        df_plot = data.copy()
    else:
        df_plot = data[data["Country"] == Country]

    trace1 = go.Bar(x=df_plot.Country, y=df_plot["Under 5"], name="Under 5")
    trace2 = go.Bar(x=df_plot.Country, y=df_plot["Aged 5-14"], name="Aged 5-14")
    trace3 = go.Bar(x=df_plot.Country, y=df_plot["Aged 15-24"], name="Aged 15-24")
    trace4 = go.Bar(x=df_plot.Country, y=df_plot["Aged 25-64"], name="Aged 25-64")
    trace5 = go.Bar(x=df_plot.Country, y=df_plot["Over 65"], name="Over 65")

    return {
        "data": [trace1, trace2, trace3, trace4, trace5],
        "layout": go.Layout(
            title="Age demographics for {}".format(Country),
            orientation=180,
            barmode="stack",
        ),
    }


@app.callback(
    dash.dependencies.Output("radar", "figure"),
    [dash.dependencies.Input("City", "value")],
)

# radar plot to compare index values
def update_radar(city):
    # creating a subset dataframe
    df = data[
        [
            "City",
            "Cost of Living Index",
            "Purchasing Power Index",
            "Safety Index",
            "Health Care Index",
            "Pollution Index",
        ]
    ]

    # categories
    cat = df.columns[1:].tolist()

    select_df = df[df["City"] == city]

    Row_list = []
    r = []
    # Iterate over each row
    for index, rows in select_df.iterrows():
        for i in range(len(cat)):
            # Create list for the current
            r.append(rows[cat[i]])

            # append the list to the final list
        Row_list.append(r)
        Row_list = list(np.concatenate(Row_list).flat)

    fig = go.Figure()

    fig.add_trace(
        go.Barpolar(
            r=Row_list,
            theta=cat,
            name="",
            marker_color=["rgb(243,203,70)"] * 6,
            marker_line_color="white",
            hoverinfo=["theta"] * 9,
            opacity=0.7,
            base=0,
        )
    )

    fig.add_trace(
        go.Barpolar(
            r=df.mean(axis=0).tolist(),
            theta=cat,
            name="",
            marker_color=["#d1d1cf"] * 6,
            marker_line_color="white",
            hoverinfo=["theta"] * 9,
            opacity=0.7,
            base=0,
        )
    )

    fig.update_layout(
        title="",
        font_size=12,
        polar=dict(
            bgcolor="rgba(0,0,0,0)",
            angularaxis=dict(linewidth=3, showline=False, showticklabels=True),
            radialaxis=dict(
                showline=False,
                showticklabels=False,
                linewidth=2,
                gridcolor="white",
                gridwidth=2,
            ),
        ),
    )

    return fig


md = data[
    [
        "City",
        "Employment",
        "Startup",
        "Tourism",
        "Housing",
        "Transport",
        "Health",
        "Food",
        "Internet Speed",
        "Access to Contraception",
        "Gender Equality",
        "Immigration Tolerance",
        "LGBT Friendly",
        "Nightscene",
        "Beer",
        "Festival",
    ]
]

for column in md.columns.tolist()[1:]:
    md["{column}_q".format(column=column)] = pd.qcut(
        md[column].rank(method="first"), 4, labels=False
    )

# Build parcats dimensions
quartiles = [
    "Startup_q",
    "Internet Speed_q",
    "Gender Equality_q",
    "Immigration Tolerance_q",
    "LGBT Friendly_q",
    "Nightscene_q",
]


dimensions = [dict(values=md[label], label=label) for label in quartiles]

# Build colorscale
color = np.zeros(len(md), dtype="uint8")
colorscale = [[0, "gray"], [1, "rgb(243,203,70)"]]

size = md["Employment"] * 30

# bubble plot for happiness and related indicators
def build_figure():

    # Build figure as FigureWidget
    fig = go.Figure(
        data=[
            go.Scatter(
                x=md["Food"],
                y=md["Health"],
                marker={"color": "gray", "size": size},
                mode="markers",
                selected={"marker": {"color": "rgb(243,203,70)"}},
                unselected={"marker": {"opacity": 0.3}},
            ),
            go.Parcats(
                domain={"y": [0, 0.4]},
                dimensions=dimensions,
                line={
                    "colorscale": colorscale,
                    "cmin": 0,
                    "cmax": 1,
                    "color": color,
                    "shape": "hspline",
                },
            ),
        ]
    )

    fig.update_layout(
        height=800,
        xaxis={"title": "Employment"},
        yaxis={"title": "Health", "domain": [0.6, 1]},
        dragmode="lasso",
        hovermode="closest",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        autosize=False,
        bargap=0.35,
    )
    return fig


app = dash.Dash(prevent_initial_callbacks=True)
app.layout = html.Div([dcc.Graph(figure=build_figure(), id="graph")])


@app.callback(
    dash.dependencies.Output("graph", "figure"),
    [
        dash.dependencies.Input("graph", "selectedData"),
        dash.dependencies.Input("graph", "clickData"),
    ],
    [State("graph", "figure")],
)

# Update color callback
def update_color(selectedData, clickData, fig):

    # trace, points, state):
    # Update scatter selection
    # fig.data[0].selectedpoints = points.point_inds

    selection = None

    # Update selection based on which event triggered the update.
    trigger = dash.callback_context.triggered[0]["prop_id"]
    if trigger == "graph.clickData":
        selection = [point["pointNumber"] for point in clickData["points"]]
    if trigger == "graph.selectedData":
        selection = [point["pointIndex"] for point in selectedData["points"]]
    # Update scatter selection
    fig["data"][0]["selectedpoints"] = selection
    # Update parcats colors
    new_color = np.zeros(len(md), dtype="uint8")
    new_color[selection] = 1
    fig["data"][1]["line"]["color"] = new_color
    return fig
    # Register callback on scatter selection...
    # fig.data[0].on_selection(update_color)
    # and parcats click
    # fig.data[1].on_click(update_color)


if __name__ == "__main__":
    app.run_server(debug=True)
