import plotly.express as px
import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import aerosandbox as asb
import aerosandbox.numpy as np
import copy
import plotly.figure_factory as ff
import pandas as pd

from app_components import *

### Build the app
app = dash.Dash(
    __name__, external_stylesheets=[dbc.themes.MINTY], title="Airfoil Analysis"
)
server = app.server

btn_style = {"margin-top": "13px"}

app.layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
                        dcc.Markdown(
                            """
                # Airfoil Analysis with [AeroSandbox](https://github.com/peterdsharpe/AeroSandbox) and [Dash](https://plotly.com/dash/)
                
                By [Peter Sharpe](https://peterdsharpe.github.io/). Uses potential flow theory (viscous effects neglected, for now). [Original repository here](https://github.com/peterdsharpe/Automotive-Airfoil-Design).
                """
                        )
                    ],
                    width=8,
                ),
                dbc.Col(
                    [
                        html.Div(
                            [
                                html.A(
                                    dbc.Button(
                                        "Enterprise Demo",
                                        color="primary",
                                        size="md",
                                        className="mr-1",
                                    ),
                                    href="https://plotly.com/get-demo/",
                                    target="_blank",
                                    style=btn_style,
                                ),
                                html.A(
                                    dbc.Button(
                                        "Source Code",
                                        size="md",
                                        className="mr-1",
                                        color="secondary",
                                    ),
                                    href="https://github.com/plotly/dash-sample-apps/tree/main/apps/dash-airfoil-design",
                                ),
                                html.Img(
                                    src="assets/MIT-logo-red-gray-72x38.svg",
                                    alt="MIT Logo",
                                    height="100%",
                                    style={"margin-left": "15px"},
                                ),
                            ],
                            style={
                                "float": "right",
                                "height": "80px",
                                "padding-bottom": "20px",
                            },
                        ),
                    ],
                    width=4,
                ),
            ],
            align="end",
        ),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Button(
                            "Modify Operating Conditions", id="operating_button"
                        ),
                        dbc.Collapse(
                            dbc.Card(dbc.CardBody(operating_slider_components,)),
                            id="operating_collapse",
                            is_open=False,
                        ),
                        html.Hr(),
                        dbc.Button(
                            "Modify Shape Parameters (Kulfan)", id="shape_button"
                        ),
                        dbc.Collapse(
                            dbc.Card(dbc.CardBody(kulfan_slider_components,)),
                            id="shape_collapse",
                            is_open=False,
                        ),
                        html.Hr(),
                        dbc.Button(
                            "Show Raw Coordinates (*.dat format)",
                            id="coordinates_button",
                        ),
                        dbc.Collapse(
                            dbc.Card(
                                dbc.CardBody(dcc.Markdown(id="coordinates_output"))
                            ),
                            id="coordinates_collapse",
                            is_open=False,
                        ),
                        html.Hr(),
                        dcc.Markdown("##### Commands"),
                        dbc.Button(
                            "Analyze",
                            id="analyze",
                            color="primary",
                            style={"margin": "5px"},
                        ),
                        html.Hr(),
                        dcc.Markdown("##### Aerodynamic Performance"),
                        dbc.Spinner(html.P(id="text_output"), color="primary",),
                    ],
                    width=3,
                ),
                dbc.Col(
                    [dcc.Graph(id="display", style={"height": "90vh"}),],
                    width=9,
                    align="start",
                ),
            ]
        ),
        html.Hr(),
        dcc.Markdown(
            """
        Aircraft design tools powered by [AeroSandbox](https://github.com/peterdsharpe/AeroSandbox). Build beautiful UIs for your scientific computing apps with [Plot.ly](https://plotly.com/) and [Dash](https://plotly.com/dash/)!
        """
        ),
    ],
    fluid=True,
)


### Callback to make shape parameters menu expand
@app.callback(
    Output("shape_collapse", "is_open"),
    [Input("shape_button", "n_clicks")],
    [State("shape_collapse", "is_open")],
)
def toggle_shape_collapse(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open


### Callback to make operating parameters menu expand
@app.callback(
    Output("operating_collapse", "is_open"),
    [Input("operating_button", "n_clicks")],
    [State("operating_collapse", "is_open")],
)
def toggle_shape_collapse(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open


### Callback to make coordinates menu expand
@app.callback(
    Output("coordinates_collapse", "is_open"),
    [Input("coordinates_button", "n_clicks")],
    [State("coordinates_collapse", "is_open")],
)
def toggle_shape_collapse(n_clicks, is_open):
    if n_clicks:
        return not is_open
    return is_open


### Callback to make operating sliders display proper values
@app.callback(
    Output("alpha_slider_output", "children"),
    [Input("alpha_slider_input", "drag_value")],
)
def display_alpha_slider(drag_value):
    return f"Angle of Attack: {drag_value}"


@app.callback(
    Output("height_slider_output", "children"),
    [Input("height_slider_input", "drag_value")],
)
def display_alpha_slider(drag_value):
    return f"Height: {drag_value}"


@app.callback(
    Output("streamline_density_slider_output", "children"),
    [Input("streamline_density_slider_input", "drag_value")],
)
def display_streamline_density_slider(drag_value):
    return f"Streamline Density: {drag_value}"


### The callback to make the kulfan sliders display proper values
for side in sides:
    for i in range(n_kulfan_inputs_per_side):

        @app.callback(
            Output(f"kulfan_{side.lower()}_{i}_output", "children"),
            [Input(f"kulfan_{side.lower()}_{i}_input", "drag_value")],
        )
        def display_slider_value(drag_value):
            return f"Parameter: {drag_value}"


def make_table(dataframe):
    return dbc.Table.from_dataframe(
        dataframe, bordered=True, hover=True, responsive=True, striped=True, style={}
    )


last_analyze_timestamp = None

n_clicks_last = 0

### The callback to draw the airfoil on the graph
@app.callback(
    Output("display", "figure"),
    Output("text_output", "children"),
    Output("coordinates_output", "children"),
    [
        Input("analyze", "n_clicks"),
        Input("alpha_slider_input", "value"),
        Input("height_slider_input", "value"),
        Input("streamline_density_slider_input", "value"),
        Input("operating_checklist", "value"),
    ]
    + [
        Input(f"kulfan_{side.lower()}_{i}_input", "value")
        for side in sides
        for i in range(n_kulfan_inputs_per_side)
    ],
)
def display_graph(
    n_clicks, alpha, height, streamline_density, operating_checklist, *kulfan_inputs
):
    ### Figure out if a button was pressed
    global n_clicks_last
    if n_clicks is None:
        n_clicks = 0

    analyze_button_pressed = n_clicks > n_clicks_last
    n_clicks_last = n_clicks

    ### Parse the checklist
    ground_effect = "ground_effect" in operating_checklist

    ### Start constructing the figure
    airfoil = asb.Airfoil(
        coordinates=asb.get_kulfan_coordinates(
            lower_weights=np.array(kulfan_inputs[n_kulfan_inputs_per_side:]),
            upper_weights=np.array(kulfan_inputs[:n_kulfan_inputs_per_side]),
            TE_thickness=0,
            enforce_continuous_LE_radius=False,
            n_points_per_side=200,
        )
    )

    ### Do coordinates output
    coordinates_output = "\n".join(
        ["```"]
        + ["AeroSandbox Airfoil"]
        + ["\t%f\t%f" % tuple(coordinate) for coordinate in airfoil.coordinates]
        + ["```"]
    )

    ### Continue doing the airfoil things
    airfoil = airfoil.rotate(angle=-np.radians(alpha))
    airfoil = airfoil.translate(0, height + 0.5 * np.sind(alpha))
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=airfoil.x(),
            y=airfoil.y(),
            mode="lines",
            name="Airfoil",
            fill="toself",
            line=dict(color="blue"),
        )
    )

    ### Default text output
    text_output = 'Click "Analyze" to compute aerodynamics!'

    xrng = (-0.5, 1.5)
    yrng = (-0.6, 0.6) if not ground_effect else (0, 1.2)

    if analyze_button_pressed:

        analysis = asb.AirfoilInviscid(
            airfoil=airfoil.repanel(50),
            op_point=asb.OperatingPoint(velocity=1, alpha=0,),
            ground_effect=ground_effect,
        )

        x = np.linspace(*xrng, 100)
        y = np.linspace(*yrng, 100)
        X, Y = np.meshgrid(x, y)
        u, v = analysis.calculate_velocity(x_field=X.flatten(), y_field=Y.flatten())
        U = u.reshape(X.shape)
        V = v.reshape(Y.shape)

        streamline_fig = ff.create_streamline(
            x,
            y,
            U,
            V,
            arrow_scale=1e-16,
            density=streamline_density,
            line=dict(color="#ff82a3"),
            name="Streamlines",
        )

        fig = go.Figure(data=streamline_fig.data + fig.data)

        text_output = make_table(
            pd.DataFrame(
                {"Engineering Quantity": ["C_L"], "Value": [f"{analysis.Cl:.3f}"]}
            )
        )

    fig.update_layout(
        xaxis_title="x/c",
        yaxis_title="y/c",
        showlegend=False,
        yaxis=dict(scaleanchor="x", scaleratio=1),
        margin={"t": 0},
        title=None,
    )

    fig.update_xaxes(range=xrng)
    fig.update_yaxes(range=yrng)

    return fig, text_output, [coordinates_output]


if __name__ == "__main__":
    app.run_server(debug=True)
