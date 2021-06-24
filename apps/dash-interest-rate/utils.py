import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc

def Header(name, app):
    title = html.H2(name, style={"margin-top": 7})
    logo = html.Img(
        src=app.get_asset_url("dash-logo.png"), style={"float": "right", "height": 60}
    )
    link = html.A(logo, href="https://plotly.com/dash/")
    btn_style = {'margin-top': '13px', 'float': 'right', 'margin-right': '10px'}
    demo_btn = html.A(dbc.Button("Enterprise Demo", style=btn_style, color="primary"), href="https://plotly.com/get-demo/")
    code_btn = html.A(dbc.Button("Source Code", style=btn_style, color="secondary"), href="https://github.com/plotly/dash-sample-apps/tree/main/apps/dash-interest-rate")

    return dbc.Row([dbc.Col(title, md=7), dbc.Col([link, demo_btn, code_btn], md=5)])

def OptionMenu(values, label, **kwargs):
    options = [{"label": s.replace("_", " ").capitalize(), "value": s} for s in values]
    kwargs["value"] = kwargs.get("value", values[0])

    if len(options) <= 4:
        component = dbc.RadioItems
        kwargs["inline"] = True
    else:
        component = dbc.Select

    return dbc.FormGroup([dbc.Label(label), component(options=options, **kwargs)])


def CustomRangeSlider(values, label, **kwargs):
    values = sorted(values)
    marks = {i: f"{i//1000}k" for i in values}

    return dbc.FormGroup(
        [
            dbc.Label(label),
            dcc.RangeSlider(
                min=values[0],
                max=values[-1],
                step=1000,
                value=[values[1], values[-2]],
                marks=marks,
                **kwargs,
            ),
        ]
    )


def get_unique(connection, db, col):
    query = f"""
    SELECT DISTINCT {col}
    FROM {db}.PUBLIC.LOAN_CLEAN;
    """
    return [x[0] for x in connection.execute(query).fetchall()]


def get_range(connection, db, col):
    query = f"""
    SELECT MIN({col}), MAX({col})
    FROM {db}.PUBLIC.LOAN_CLEAN;
    """
    return connection.execute(query).fetchall()[0]
