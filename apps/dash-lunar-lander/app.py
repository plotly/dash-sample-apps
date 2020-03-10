import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import numpy as np
import casadi
from casadi import SX, DM
from math import cos, sin
import plotly.graph_objects as go

app = dash.Dash(__name__)
server = app.server

app.layout = html.Div(
    children=[
        html.Div(
            children=[
                html.Img(src=app.get_asset_url("dash-logo.png"), id="logo"),
                html.H4("Dash Lunar Lander"),
                html.A(
                    id="github-link",
                    children=["View source code on Github"],
                    href="https://github.com/plotly/dash-sample-apps/blob/master/apps/dash-lunar-lander/app.py",
                ),
            ],
            className="app-header",
        ),
        html.Div(
            [
                html.H3("Lunar Lander Trade Study App"),
                html.P(
                    "This app allows for exploring the trade space of a lunar lander. You can adjust different "
                    "parameters on the lander, like max angular velocity and inital mass, using the sliders and a new "
                    "optimal trajectory will be recomputed in near real time using Casadi."
                ),
                html.A(
                    "The lunar lander model was based off of one in this paper,",
                    href="https://arxiv.org/pdf/1610.08668.pdf",
                    target="_blank",
                ),
                html.P(
                    "This model does not allow for the lander to gimbal its engine, instead it must turn the entire "
                    "spacecraft and then thrust to cancel out any horizontal velocity, giving the mass optimal control "
                    "case a distinctive hooked shape. Switching the optimizer to time-optimal control results in a "
                    "smoother and more expected shape."
                ),
            ]
        ),
        # Display the Trajectory
        dcc.Loading(
            [dcc.Graph(id="indicator-graphic")],
            style={"height": "450px", "verticalAlign": "middle"},
            type="dot",
        ),
        # Adjust the Spacecraft Parameters
        dcc.Slider(
            id="m0Slider",
            min=5000,
            max=15000,
            value=10000,
            step=10,
            marks={5000: {"label": "Inital Mass", "style": {"transform": "none"}}},
        ),
        dcc.Slider(
            id="c1Slider",
            min=35000,
            max=44000 * 1.1,
            value=44000,
            step=10,
            marks={35000: {"label": "Max Thrust", "style": {"transform": "none"}}},
        ),
        dcc.Slider(
            id="isp",
            min=300,
            max=330,
            value=311,
            step=1,
            marks={300: {"label": "ISP", "style": {"transform": "none"}}},
        ),
        dcc.Slider(
            id="c3Slider",
            min=0.01,
            max=0.1,
            value=0.0698,
            step=0.0001,
            marks={
                0.01: {"label": "Max Angular Velocity", "style": {"transform": "none"}}
            },
        ),
        dcc.Slider(
            id="gSlider",
            min=0.5,
            max=2,
            value=1.6229,
            step=0.0001,
            marks={0.5: {"label": "Gravity", "style": {"transform": "none"}}},
        ),
        html.Div(
            [
                # Display Current Spacecraft Parameters
                html.Div(
                    [
                        html.H5("Parameters"),
                        html.P("Inital Mass (kg):", id="m0Out"),
                        html.P("Max Thrust (N):", id="maThrustOut"),
                        html.P("ISP (s):", id="ispOut"),
                        html.P("C3 (rad/s):", id="c3Out"),
                        html.P("Gravity (N):", id="gravOut"),
                    ],
                ),
                html.Div(
                    [
                        # Choose Between Different Cost Functions
                        html.H5("Outputs"),
                        dcc.RadioItems(
                            options=[
                                {"label": "Mass Optimal Control", "value": 1},
                                {"label": "Time Optimal Control", "value": 2},
                            ],
                            id="costFun",
                            value=1,
                        ),
                        # Display Final Cost Functions
                        html.P("Final  Mass (kg):", id="mfOut"),
                        html.P("TOF (s):", id="tof"),
                    ],
                ),
                html.Div(
                    [
                        html.H5("Adjust Initial Position"),
                        html.Div(
                            [
                                # DPad for Adjusting the Spacecrafts Initial Position
                                html.Div(
                                    [html.Button("Left", id="Left",),],
                                    className="direction-button",
                                ),
                                html.Div(
                                    [
                                        html.Button("Up", id="Up",),
                                        html.Button("Down", id="Down",),
                                    ],
                                    style={
                                        "display": "flex",
                                        "flex-direction": "column",
                                    },
                                    className="direction-button",
                                ),
                                html.Div(
                                    [html.Button("Right", id="Right",),],
                                    className="direction-button",
                                ),
                            ],
                            id="button-grid",
                        ),
                    ],
                    id="adjust-init-pos",
                ),
            ],
            id="bottom-cards",
        ),
    ],
    className="entire-app",
)


@app.callback(
    [
        Output("indicator-graphic", "figure"),
        Output("m0Out", "children"),
        Output("mfOut", "children"),
        Output("maThrustOut", "children"),
        Output("ispOut", "children"),
        Output("c3Out", "children"),
        Output("gravOut", "children"),
        Output("tof", "children"),
    ],
    [
        Input(component_id="m0Slider", component_property="value"),
        Input(component_id="c1Slider", component_property="value"),
        Input(component_id="isp", component_property="value"),
        Input(component_id="c3Slider", component_property="value"),
        Input(component_id="gSlider", component_property="value"),
        Input(component_id="Left", component_property="n_clicks"),
        Input(component_id="Right", component_property="n_clicks"),
        Input(component_id="Up", component_property="n_clicks"),
        Input(component_id="Down", component_property="n_clicks"),
        Input(component_id="costFun", component_property="value"),
    ],
)
def update_output_div(
    m0, c1, isp, c3, g, n_clicksL, n_clicksR, n_clicksU, n_clicksD, optimal
):

    # Adjust the Spacecraft's Intal Position
    if n_clicksL is None:
        n_clicksL = 0
    if n_clicksR is None:
        n_clicksR = 0
    if n_clicksU is None:
        n_clicksU = 0
    if n_clicksD is None:
        n_clicksD = 0

    motionPerClickX = 10
    motionPerClickY = 100

    clicksX = n_clicksR - n_clicksL
    clicksY = n_clicksU - n_clicksD
    x0 = -30 + (motionPerClickX * clicksX)
    y0 = 1000 + (motionPerClickY * clicksY)

    # Setup and Solve the Optimal Control Problem
    IC = np.array([x0, y0, 10, 20, m0, 0])
    g0 = 9.81  # For Calculating Spacecraft Engine Efficiency
    argsRWSC = [c1, isp, g0, c3, g]
    resultsRWSC = wrapperRWSC(IC, argsRWSC, optimal)

    # Unpack the Solution
    states = np.array(resultsRWSC["states"])
    dt = resultsRWSC["dt"]
    cont = np.array(resultsRWSC["controls"])
    tof = 30 * dt
    x = states[:, 0]
    y = states[:, 1]

    # Create the Graphics Object
    trace0 = go.Scatter(y=y, x=x, mode="lines", name="Trajectory")
    trace1 = go.Scatter(
        y=np.array(y[0]), x=np.array(x[0]), mode="markers", name="Inital Location"
    )
    trace2 = go.Scatter(
        y=np.array(y[-1]), x=np.array(x[-1]), mode="markers", name="Target"
    )
    data = [trace1, trace2, trace0]

    if optimal == 1:
        titleName = "Lunar Lander - Mass Optimal Trajectory"
    else:
        titleName = "Lunar Lander - Time Optimal Trajectory"

    layout = go.Layout(
        title=titleName,
        yaxis={"title": "Height (meters)", "scaleanchor": "x", "scaleratio": 1},
        xaxis={"title": "Distance uprange from target (meters)"},
    )
    return (
        {"data": data, "layout": layout},
        "Inital Mass: " + str(m0) + " (kg)",
        "Final  Mass: " + str(np.round(states[-1, 4], 2)) + " (kg)",
        "Max Thrust: " + str(c1) + " (N)",
        "ISP: " + str(isp) + " (S)",
        "Max c3: " + str(c3) + "(rad/s)",
        "Gravity: " + str(g) + "(N)",
        "TOF: " + str(tof) + "(s)",
    )


def RWSC(states, controls, args):
    # ODE for the Reaction Wheel Spacecraft
    # See https://arxiv.org/pdf/1610.08668.pdf Equation 13 for EOM

    # Unpacking the Input Variables
    x = states[:, 0]
    y = states[:, 1]
    xDot = states[:, 2]
    yDot = states[:, 3]
    m = states[:, 4]
    scTheta = states[:, 5]

    beta = controls[:, 0]
    u = controls[:, 1]

    c1 = args[0]
    isp = args[1]
    g0 = args[2]
    c2 = isp * g0
    c3 = args[3]
    g = args[4]

    # Equations of Motion
    xDotDot = c1 * u * np.sin(scTheta) / m
    yDotDot = (c1 * u * np.cos(scTheta) / m) - g
    scThetaDot = c3 * beta
    mDot = -c1 * u / c2

    # Repacking Time Derivative of the State Vector
    statesDot = 0 * states
    statesDot[:, 0] = xDot
    statesDot[:, 1] = yDot
    statesDot[:, 2] = xDotDot
    statesDot[:, 3] = yDotDot
    statesDot[:, 4] = mDot
    statesDot[:, 5] = scThetaDot
    return statesDot


def wrapperRWSC(IC, args, optimal):
    # Converting the Optimal Control Problem into a Non-Linear Programming  Problem

    numStates = 6
    numInputs = 2
    nodes = 30  # Keep this Number Small to Reduce Runtime

    dt = SX.sym("dt")
    states = SX.sym("state", nodes, numStates)
    controls = SX.sym("controls", nodes, numInputs)

    variables_list = [dt, states, controls]
    variables_name = ["dt", "states", "controls"]
    variables_flat = casadi.vertcat(*[casadi.reshape(e, -1, 1) for e in variables_list])
    pack_variables_fn = casadi.Function(
        "pack_variables_fn", variables_list, [variables_flat], variables_name, ["flat"]
    )
    unpack_variables_fn = casadi.Function(
        "unpack_variables_fn",
        [variables_flat],
        variables_list,
        ["flat"],
        variables_name,
    )

    # Bounds
    bds = [
        [np.sqrt(np.finfo(float).eps), np.inf],
        [-100, 300],
        [0, np.inf],
        [-np.inf, np.inf],
        [-np.inf, np.inf],
        [np.sqrt(np.finfo(float).eps), np.inf],
        [-1, 1],
        [np.sqrt(np.finfo(float).eps), 1],
    ]

    lower_bounds = unpack_variables_fn(flat=-float("inf"))
    lower_bounds["dt"][:, :] = bds[0][0]

    lower_bounds["states"][:, 0] = bds[1][0]
    lower_bounds["states"][:, 1] = bds[2][0]
    lower_bounds["states"][:, 4] = bds[5][0]

    lower_bounds["controls"][:, 0] = bds[6][0]
    lower_bounds["controls"][:, 1] = bds[7][0]

    upper_bounds = unpack_variables_fn(flat=float("inf"))
    upper_bounds["dt"][:, :] = bds[0][1]

    upper_bounds["states"][:, 0] = bds[1][1]

    upper_bounds["controls"][:, 0] = bds[6][1]
    upper_bounds["controls"][:, 1] = bds[7][1]

    # Set Initial Conditions
    #   Casadi does not accept equality constraints, so boundary constraints are
    #   set as box constraints with 0 area.
    lower_bounds["states"][0, 0] = IC[0]
    lower_bounds["states"][0, 1] = IC[1]
    lower_bounds["states"][0, 2] = IC[2]
    lower_bounds["states"][0, 3] = IC[3]
    lower_bounds["states"][0, 4] = IC[4]
    lower_bounds["states"][0, 5] = IC[5]
    upper_bounds["states"][0, 0] = IC[0]
    upper_bounds["states"][0, 1] = IC[1]
    upper_bounds["states"][0, 2] = IC[2]
    upper_bounds["states"][0, 3] = IC[3]
    upper_bounds["states"][0, 4] = IC[4]
    upper_bounds["states"][0, 5] = IC[5]

    # Set Final Conditions
    #   Currently set for a soft touchdown at the origin
    lower_bounds["states"][-1, 0] = 0
    lower_bounds["states"][-1, 1] = 0
    lower_bounds["states"][-1, 2] = 0
    lower_bounds["states"][-1, 3] = 0
    lower_bounds["states"][-1, 5] = 0
    upper_bounds["states"][-1, 0] = 0
    upper_bounds["states"][-1, 1] = 0
    upper_bounds["states"][-1, 2] = 0
    upper_bounds["states"][-1, 3] = 0
    upper_bounds["states"][-1, 5] = 0

    # Initial Guess Generation
    #   Generate the initial guess as a line between initial and final conditions
    xIG = np.array(
        [
            np.linspace(IC[0], 0, nodes),
            np.linspace(IC[1], 0, nodes),
            np.linspace(IC[2], 0, nodes),
            np.linspace(IC[3], 0, nodes),
            np.linspace(IC[4], IC[4] * 0.5, nodes),
            np.linspace(IC[5], 0, nodes),
        ]
    ).T
    uIG = np.array([np.linspace(0, 1, nodes), np.linspace(1, 1, nodes)]).T

    ig_list = [60 / nodes, xIG, uIG]
    ig_flat = casadi.vertcat(*[casadi.reshape(e, -1, 1) for e in ig_list])

    # Generating Defect Vector
    xLow = states[0 : (nodes - 1), :]
    xHigh = states[1:nodes, :]

    contLow = controls[0 : (nodes - 1), :]
    contHi = controls[1:nodes, :]
    contMid = (contLow + contHi) / 2

    # Use a RK4 Method for Generating the Defects
    k1 = RWSC(xLow, contLow, args)
    k2 = RWSC(xLow + (0.5 * dt * k1), contMid, args)
    k3 = RWSC(xLow + (0.5 * dt * k2), contMid, args)
    k4 = RWSC(xLow + k3, contHi, args)
    xNew = xLow + ((dt / 6) * (k1 + (2 * k2) + (2 * k3) + k4))
    defect = casadi.reshape(xNew - xHigh, -1, 1)

    # Choose the Cost Function
    if optimal == 1:
        J = -states[-1, 4]  # Mass Optimal Cost Function
    elif optimal == 2:
        J = dt[0]  # Time Optimal Cost Function

    # Put the OCP into a form that Casadi can solve
    nlp = {"x": variables_flat, "f": J, "g": defect}
    solver = casadi.nlpsol(
        "solver", "bonmin", nlp
    )  # Use bonmin instead of ipopt due to speed
    result = solver(
        x0=ig_flat,
        lbg=0.0,
        ubg=0.0,
        lbx=pack_variables_fn(**lower_bounds)["flat"],
        ubx=pack_variables_fn(**upper_bounds)["flat"],
    )

    results = unpack_variables_fn(flat=result["x"])

    return results


if __name__ == "__main__":
    app.run_server(debug=True)
