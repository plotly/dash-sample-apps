from dash import html, dcc
from textwrap import dedent

def Header(app):
    app_name = html.Span("Support Vector Machine (SVM) Explorer")
    logo = html.Img(src=app.get_asset_url("images/plotly-logo.png"))
    link = html.A(logo, href="https://plotly.com/dash/", target="_blank") 
    demo_link = html.A("ENTERPRISE DEMO", href="https://plotly.com/get-demo/", target="_blank", className="demo-button")
    return html.Div([html.Div(app_name, className="header-title"), html.Div([demo_link, link], className="header-logos")], className="header")


def _omit(omitted_keys, d):
    return {k: v for k, v in d.items() if k not in omitted_keys}

def FormattedSlider(**kwargs):
    return html.Div(
        style=kwargs.get("style", {}), children=dcc.Slider(**_omit(["style"], kwargs))
    )

def NamedSlider(name, **kwargs):
    return html.Div(
        style={"padding": "20px 10px 25px 4px"},
        children=[
            html.P(f"{name}:"),
            html.Div(style={"margin-left": "6px"}, children=dcc.Slider(**kwargs)),
        ],
    )

def NamedDropdown(name, **kwargs):
    return html.Div(
        style={"margin": "10px 0px"},
        children=[
            html.P(children=f"{name}:", style={"margin-left": "3px"}),
            dcc.Dropdown(**kwargs),
        ],
    )


def NamedRadioItems(name, **kwargs):
    return html.Div(
        style={"padding": "20px 10px 25px 4px"},
        children=[html.P(children=f"{name}:"), dcc.RadioItems(**kwargs)],
    )


controls_first = children=[
    NamedDropdown(
        name="Select Dataset",
        id="dropdown-select-dataset",
        options=[
            {"label": "Moons", "value": "moons"},
            {
                "label": "Linearly Separable",
                "value": "linear",
            },
            {
                "label": "Circles",
                "value": "circles",
            },
        ],
        clearable=False,
        searchable=False,
        value="moons",
    ),
    NamedSlider(
        name="Sample Size",
        id="slider-dataset-sample-size",
        min=100,
        max=500,
        step=100,
        marks={
            str(i): str(i)
            for i in [100, 200, 300, 400, 500]
        },
        value=300,
    ),
    NamedSlider(
        name="Noise Level",
        id="slider-dataset-noise-level",
        min=0,
        max=1,
        marks={
            i / 10: str(i / 10)
            for i in range(0, 11, 2)
        },
        step=0.1,
        value=0.2,
    ),
]

controls_second = [
    NamedSlider(
        name="Threshold",
        id="slider-threshold",
        min=0,
        max=1,
        value=0.4,
        step=0.01,
        marks={ i / 10: str(i / 10) for i in range(0, 11, 2) },
    ),
    html.Button(
        "Reset Threshold",
        id="button-zero-threshold",
    ),
]

controls_third = [
    NamedDropdown(
        name="Kernel",
        id="dropdown-svm-parameter-kernel",
        options=[
            {
                "label": "Radial basis function (RBF)",
                "value": "rbf",
            },
            {"label": "Linear", "value": "linear"},
            {
                "label": "Polynomial",
                "value": "poly",
            },
            {
                "label": "Sigmoid",
                "value": "sigmoid",
            },
        ],
        value="rbf",
        clearable=False,
        searchable=False,
    ),
    NamedSlider(
        name="Cost (C)",
        id="slider-svm-parameter-C-power",
        min=-2,
        max=4,
        value=0,
        marks={ i: str(10 ** i) for i in range(-2, 5, 2) },
    ),
    FormattedSlider(
        id="slider-svm-parameter-C-coef",
        min=1,
        max=9,
        value=1,
    ),
    NamedSlider(
        name="Degree",
        id="slider-svm-parameter-degree",
        min=2,
        max=10,
        value=3,
        step=1,
        marks={
            str(i): str(i) for i in range(2, 11, 2)
        },
    ),
    NamedSlider(
        name="Gamma",
        id="slider-svm-parameter-gamma-power",
        min=-5,
        max=0,
        value=-1,
        marks={
            i: "{}".format(10 ** i)
            for i in range(-5, 1, 2)
        },
    ),
    FormattedSlider(
        id="slider-svm-parameter-gamma-coef",
        min=1,
        max=9,
        value=5,
    ),
    NamedRadioItems(
        name="Shrinking",
        id="radio-svm-parameter-shrinking",
        labelStyle={
            "margin-right": "7px",
            "display": "inline-block",
        },
        options=[
            {
                "label": " Enabled",
                "value": "True",
            },
            {
                "label": " Disabled",
                "value": "False",
            },
        ],
        value="True",
    ),
]