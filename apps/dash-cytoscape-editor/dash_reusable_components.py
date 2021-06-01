import dash_core_components as dcc
import dash_html_components as html


# Display utility functions
def _merge(a, b):
    return dict(a, **b)


def _omit(omitted_keys, d):
    return {k: v for k, v in d.items() if k not in omitted_keys}


# Custom Display Components
def Card(children, **kwargs):
    return html.Section(
        children,
        style=_merge(
            {
                "padding": 20,
                "margin": 5,
                "borderRadius": 5,
                "border": "thin lightgrey solid",
                "background-color": "white",
                # Remove possibility to select the text for better UX
                "user-select": "none",
                "-moz-user-select": "none",
                "-webkit-user-select": "none",
                "-ms-user-select": "none",
            },
            kwargs.get("style", {}),
        ),
        **_omit(["style"], kwargs),
    )


def SectionTitle(title, size, align="center", color="#222"):
    return html.Div(
        style={"text-align": align, "color": color},
        children=dcc.Markdown("#" * size + " " + title),
    )


def NamedCard(title, size, children, **kwargs):
    size = min(size, 6)
    size = max(size, 1)

    return html.Div(
        [Card([SectionTitle(title, size, align="left")] + children, **kwargs)]
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


def NamedInput(name, **kwargs):
    return html.Div(children=[html.P(children=f"{name}:"), dcc.Input(**kwargs)])


# Utils
def DropdownOptionsList(*args):
    return [{"label": val.capitalize(), "value": val} for val in args]
