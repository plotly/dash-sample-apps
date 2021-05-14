import plotly.express as px
import plotly.graph_objects as go
import plotly.subplots as sub
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import aerosandbox as asb
import casadi as cas
from airplane import make_airplane
import numpy as np
import pandas as pd

app = dash.Dash(external_stylesheets=[dbc.themes.MINTY])
app.title = "Aircraft CFD"
server = app.server

app.layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H2("Solar Aircraft Design with AeroSandbox and Dash"),
                        html.H5("Peter Sharpe"),
                    ],
                    width=True,
                ),
                # dbc.Col([
                #     html.Img(src="assets/MIT-logo-red-gray-72x38.svg", alt="MIT Logo", height="30px"),
                # ], width=1)
            ],
            align="end",
        ),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Div(
                            [
                                html.H5("Key Parameters"),
                                html.P("Number of booms:"),
                                dcc.Slider(
                                    id="n_booms",
                                    min=1,
                                    max=3,
                                    step=1,
                                    value=3,
                                    marks={1: "1", 2: "2", 3: "3",},
                                ),
                                html.P("Wing Span [m]:"),
                                dcc.Input(id="wing_span", value=43, type="number"),
                                html.P("Angle of Attack [deg]:"),
                                dcc.Input(id="alpha", value=7.0, type="number"),
                            ]
                        ),
                        html.Hr(),
                        html.Div(
                            [
                                html.H5("Commands"),
                                dbc.Button(
                                    "Display (1s)",
                                    id="display_geometry",
                                    color="primary",
                                    style={"margin": "5px"},
                                    n_clicks_timestamp="0",
                                ),
                                dbc.Button(
                                    "LL Analysis (3s)",
                                    id="run_ll_analysis",
                                    color="secondary",
                                    style={"margin": "5px"},
                                    n_clicks_timestamp="0",
                                ),
                                dbc.Button(
                                    "VLM Analysis (15s)",
                                    id="run_vlm_analysis",
                                    color="secondary",
                                    style={"margin": "5px"},
                                    n_clicks_timestamp="0",
                                ),
                            ]
                        ),
                        html.Hr(),
                        html.Div(
                            [
                                html.H5("Aerodynamic Performance"),
                                dbc.Spinner(html.P(id="output"), color="primary",),
                            ]
                        ),
                    ],
                    width=3,
                ),
                dbc.Col(
                    [
                        # html.Div(id='display')
                        dbc.Spinner(
                            dcc.Graph(id="display", style={"height": "80vh"}),
                            color="primary",
                        )
                    ],
                    width=True,
                ),
            ]
        ),
        html.Hr(),
        html.P(
            [
                html.A(
                    "Source code",
                    href="https://github.com/peterdsharpe/AeroSandbox-Interactive-Demo",
                ),
                ". Aircraft design tools powered by ",
                html.A(
                    "AeroSandbox", href="https://peterdsharpe.github.com/AeroSandbox"
                ),
                ". Build beautiful UIs for your scientific computing apps with ",
                html.A("Plot.ly ", href="https://plotly.com/"),
                "and ",
                html.A("Dash", href="https://plotly.com/dash/"),
                "!",
            ]
        ),
    ],
    fluid=True,
)


def make_table(dataframe):
    return dbc.Table.from_dataframe(
        dataframe, bordered=True, hover=True, responsive=True, striped=True, style={}
    )


@app.callback(
    [Output("display", "figure"), Output("output", "children")],
    [
        Input("display_geometry", "n_clicks_timestamp"),
        Input("run_ll_analysis", "n_clicks_timestamp"),
        Input("run_vlm_analysis", "n_clicks_timestamp"),
    ],
    [State("n_booms", "value"), State("wing_span", "value"), State("alpha", "value"),],
)
def display_geometry(
    display_geometry, run_ll_analysis, run_vlm_analysis, n_booms, wing_span, alpha,
):
    ### Figure out which button was clicked
    try:
        button_pressed = np.argmax(
            np.array(
                [
                    float(display_geometry),
                    float(run_ll_analysis),
                    float(run_vlm_analysis),
                ]
            )
        )
        assert button_pressed is not None
    except:
        button_pressed = 0

    ### Make the airplane
    airplane = make_airplane(n_booms=n_booms, wing_span=wing_span,)
    op_point = asb.OperatingPoint(density=0.10, velocity=20, alpha=alpha,)
    if button_pressed == 0:
        # Display the geometry
        figure = airplane.draw(show=False, colorbar_title=None)
        output = "Please run an analysis to display the data."
    elif button_pressed == 1:
        # Run an analysis
        opti = cas.Opti()  # Initialize an analysis/optimization environment
        ap = asb.Casll1(
            airplane=airplane, op_point=op_point, opti=opti, run_setup=False
        )
        ap.setup(verbose=False)
        # Solver options
        p_opts = {}
        s_opts = {}
        # s_opts["mu_strategy"] = "adaptive"
        opti.solver("ipopt", p_opts, s_opts)
        # Solve
        try:
            sol = opti.solve()
        except RuntimeError:
            sol = opti.debug
            raise Exception("An error occurred!")

        figure = ap.draw(show=False)  # Generates figure

        output = make_table(
            pd.DataFrame(
                {
                    "Figure": ["CL", "CD", "CDi", "CDp", "L/D"],
                    "Value": [
                        sol.value(ap.CL),
                        sol.value(ap.CD),
                        sol.value(ap.CDi),
                        sol.value(ap.CDp),
                        sol.value(ap.CL / ap.CD),
                    ],
                }
            )
        )

    elif button_pressed == 2:
        # Run an analysis
        opti = cas.Opti()  # Initialize an analysis/optimization environment
        ap = asb.Casvlm1(
            airplane=airplane, op_point=op_point, opti=opti, run_setup=False
        )
        ap.setup(verbose=False)
        # Solver options
        p_opts = {}
        s_opts = {}
        # s_opts["mu_strategy"] = "adaptive"
        opti.solver("ipopt", p_opts, s_opts)
        # Solve
        try:
            sol = opti.solve()
        except RuntimeError:
            sol = opti.debug
            raise Exception("An error occurred!")

        figure = ap.draw(show=False)  # Generates figure

        output = make_table(
            pd.DataFrame(
                {
                    "Figure": ["CL", "CDi", "L/Di"],
                    "Value": [
                        sol.value(ap.CL),
                        sol.value(ap.CDi),
                        sol.value(ap.CL / ap.CDi),
                    ],
                }
            )
        )

    figure.update_layout(
        autosize=True,
        # width=1000,
        # height=700,
        margin=dict(l=0, r=0, b=0, t=0,),
    )

    return (figure, output)


try:  # wrapping this, since a forum post said it may be deprecated at some point.
    app.title = "Aircraft Design with Dash"
except:
    print("Could not set the page title!")


if __name__ == "__main__":
    app.run_server(debug=False)
