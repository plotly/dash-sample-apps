import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc


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
