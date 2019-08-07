import math
import dash
import dash_html_components as html

from sfManager import sf_Manager

app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)

app.config.suppress_callback_exceptions = True

sf_manager = sf_Manager()

millnames = ["", " K", " M", " B", " T"]  # used to convert numbers


# return html Table with dataframe values
def df_to_table(df):
    return html.Table(
        [html.Tr([html.Th(col) for col in df.columns])]
        + [
            html.Tr([html.Td(df.iloc[i][col]) for col in df.columns])
            for i in range(len(df))
        ]
    )


# returns most significant part of a number
def millify(n):
    n = float(n)
    millidx = max(
        0,
        min(
            len(millnames) - 1, int(math.floor(0 if n == 0 else math.log10(abs(n)) / 3))
        ),
    )

    return "{:.0f}{}".format(n / 10 ** (3 * millidx), millnames[millidx])


# returns top indicator div
def indicator(color, text, id_value):
    return html.Div(
        [
            html.P(id=id_value, className="indicator_value"),
            html.P(text, className="twelve columns indicator_text"),
        ],
        className="four columns indicator pretty_container",
    )
