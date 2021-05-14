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


def get_unique(connection, db, table, col):
    query = f"""
    SELECT DISTINCT {col}
    FROM {db}.dbo.{table};
    """
    return [x[0] for x in connection.execute(query).fetchall()]


def get_range(connection, db, table, col):
    query = f"""
    SELECT MIN({col}), MAX({col})
    FROM {db}.dbo.{table};
    """
    return connection.execute(query).fetchall()[0]


def get_column_strings(df):
    # Load the actual csv file

    # Create SQL columns based on the columns of that dataframe
    types = (
        df.dtypes.copy()
        .replace("float64", "FLOAT")
        .replace("int64", "INT")
        .replace("object", "VARCHAR(100) COLLATE Latin1_General_BIN2")
    )

    ls = [f"{ix.lower()} {t}" for ix, t in zip(types.index, types.values)]
    return ",\n".join(ls)
