import base64
from io import BytesIO

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import floris.tools as wfct
import matplotlib.pyplot as plt

import reusable_components as rc  # see reusable_components.py


# ############ Create helper functions ############
def mpl_to_b64(fig, format="png", dpi=300, **kwargs):
    b_io = BytesIO()
    fig.savefig(b_io, format=format, bbox_inches="tight", dpi=dpi, **kwargs)
    b64_enc = base64.b64encode(b_io.getvalue()).decode("utf-8")
    return f"data:image/{format};base64," + b64_enc


def build_visualizations(x_loc, y_loc, yaw_1, wd, gch, minSpeed=4, maxSpeed=8.0):
    fi = wfct.floris_interface.FlorisInterface("./data/example_input.json")
    fi.set_gch(gch)
    fi.reinitialize_flow_field(
        wind_direction=wd, layout_array=((0, 126 * 7, 126 * 14), (0, 0, 0))
    )
    fi.calculate_wake(yaw_angles=[yaw_1, 0, 0])

    # Horizontal plane
    fig, ax = plt.subplots()
    wfct.visualization.visualize_cut_plane(
        fi.get_hor_plane(), ax=ax, minSpeed=minSpeed, maxSpeed=maxSpeed
    )
    ax.axhline(y_loc, color="w", ls="--", lw=1)
    ax.axvline(x_loc, color="w", ls="--", lw=1)
    horiz_b64 = mpl_to_b64(fig)
    plt.close(fig)

    # Cross (x-normal) plane
    fig, ax = plt.subplots()
    wfct.visualization.visualize_cut_plane(
        fi.get_cross_plane(x_loc=x_loc), ax=ax, minSpeed=minSpeed, maxSpeed=maxSpeed
    )
    wfct.visualization.reverse_cut_plane_x_axis_in_plot(ax)
    x_plane_b64 = mpl_to_b64(fig)
    plt.close(fig)

    # Cross (y-normal) plane
    fig, ax = plt.subplots()
    wfct.visualization.visualize_cut_plane(
        fi.get_y_plane(y_loc=y_loc), ax=ax, minSpeed=minSpeed, maxSpeed=maxSpeed
    )
    wfct.visualization.reverse_cut_plane_x_axis_in_plot(ax)
    y_plane_b64 = mpl_to_b64(fig)
    plt.close(fig)

    return horiz_b64, x_plane_b64, y_plane_b64


# ############ Initialize app ############
app = dash.Dash(__name__, external_stylesheets=[rc.MATERALIZE_CSS])
server = app.server


# ############ Build components and layouts ############
navbar = html.Nav(
    html.Div(
        className="nav-wrapper teal",
        children=[
            html.Img(
                src=app.get_asset_url("dash-logo.png"),
                style={"float": "right", "height": "100%", "padding-right": "15px"},
            ),
            html.A(
                "GCH and Cut Plane Visualization in FLORIS",
                className="brand-logo",
                href="https://plotly.com/dash/",
                style={"padding-left": "15px"},
            ),
        ],
    )
)

controls = [
    rc.CustomSlider(id="wind-direction", min=250, max=290, label="Wind Direction"),
    rc.CustomSlider(id="yaw-angle", min=-30, max=30, label="Yaw angle T1"),
    rc.CustomSlider(
        id="x-loc", min=0, max=3000, value=500, label="X Normal Plane Intercept"
    ),
    rc.CustomSlider(id="y-loc", min=-100, max=100, label="Y Normal Plane Intercept"),
]

left_section = rc.Card(
    rc.CardContent(
        [
            rc.CardTitle("Horizontal Cut Plane"),
            html.Img(id="gch-horizontal", style={"width": "100%"}),
            rc.CardTitle("Cross (X-Normal) Cut Plane"),
            html.Img(id="gch-x-normal", style={"width": "100%"}),
            rc.CardTitle("Cross (Y-Normal) Cut Plane"),
            html.Img(id="gch-y-normal", style={"width": "100%"}),
        ]
    )
)

right_section = rc.Card(
    rc.CardContent(
        [
            rc.CardTitle("Horizontal Cut Plane"),
            html.Img(id="no-gch-horizontal", style={"width": "100%"}),
            rc.CardTitle("Cross (X-Normal) Cut Plane"),
            html.Img(id="no-gch-x-normal", style={"width": "100%"}),
            rc.CardTitle("Cross (Y-Normal) Cut Plane"),
            html.Img(id="no-gch-y-normal", style={"width": "100%"}),
        ]
    )
)

app.layout = html.Div(
    style={"--slider_active": "teal"},
    # className="container",
    children=[
        navbar,
        html.Br(),
        rc.Row(
            rc.Col(
                rc.Card(rc.CardContent(rc.Row([rc.Col(c, width=3) for c in controls]))),
                width=12,
            )
        ),
        rc.Row(
            [
                rc.Col([html.H4("Results with GCH"), left_section], width=6),
                rc.Col([html.H4("Results without GCH"), right_section], width=6),
            ]
        ),
    ],
)


@app.callback(
    Output("gch-horizontal", "src"),
    Output("gch-x-normal", "src"),
    Output("gch-y-normal", "src"),
    Input("x-loc", "value"),
    Input("y-loc", "value"),
    Input("yaw-angle", "value"),
    Input("wind-direction", "value"),
)
def gch_update(x_loc, y_loc, yaw_1, wd):
    return build_visualizations(x_loc, y_loc, yaw_1, wd, gch=True)


@app.callback(
    Output("no-gch-horizontal", "src"),
    Output("no-gch-x-normal", "src"),
    Output("no-gch-y-normal", "src"),
    Input("x-loc", "value"),
    Input("y-loc", "value"),
    Input("yaw-angle", "value"),
    Input("wind-direction", "value"),
)
def no_gch_update(x_loc, y_loc, yaw_1, wd):
    return build_visualizations(x_loc, y_loc, yaw_1, wd, gch=False)


if __name__ == "__main__":
    app.run_server(debug=True, threaded=False, processes=2)
