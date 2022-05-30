import dash
from dash import Dash, html, dcc, Input, Output, State, callback, callback_context
import dash_bootstrap_components as dbc


from utils.components import (
    header,
    key_parameters_card,
    commands_card,
    aerodynamic_performance_card,
    figure_card,
)
from utils.figures import analysis_and_display

app = dash.Dash(external_stylesheets=[dbc.themes.MINTY])
app.title = "Aircraft CFD"
server = app.server


app.layout = dbc.Container(
    [
        dbc.Row(
            header(
                app,
                "white",
                "Solar Aircraft Design with AeroSandbox and Dash",
                "Peter Sharpe",
            ),
        ),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(
                    [
                        key_parameters_card("n_booms", "wing_span", "alpha"),
                        html.Hr(),
                        commands_card(
                            "display_geometry", "run_ll_analysis", "run_vlm_analysis"
                        ),
                        html.Hr(),
                        aerodynamic_performance_card("output"),
                    ],
                    width=3,
                ),
                figure_card("display"),
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


@callback(
    Output("display", "figure"),
    Output("output", "children"),
    Input("display_geometry", "n_clicks_timestamp"),
    Input("run_ll_analysis", "n_clicks_timestamp"),
    Input("run_vlm_analysis", "n_clicks_timestamp"),
    State("n_booms", "value"),
    State("wing_span", "value"),
    State("alpha", "value"),
)
def run_display_geometry(
    display_geometry,
    run_ll_analysis,
    run_vlm_analysis,
    n_booms,
    wing_span,
    alpha,
):
    return analysis_and_display(
        display_geometry, run_ll_analysis, run_vlm_analysis, n_booms, wing_span, alpha
    )


try:  # wrapping this, since a forum post said it may be deprecated at some point.
    app.title = "Aircraft Design with Dash"
except:
    print("Could not set the page title!")


if __name__ == "__main__":
    app.run_server(debug=False)
