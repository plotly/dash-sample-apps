import dash
import pandas as pd
import pathlib
import dash_html_components as html
import dash_core_components as dcc

from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
from helpers import make_dash_table, create_plot


app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

server = app.server

DATA_PATH = pathlib.Path(__file__).parent.joinpath("data").resolve()

# read from datasheet
df = pd.read_csv(DATA_PATH.joinpath("small_molecule_drugbank.csv")).drop(
    ["Unnamed: 0"], axis=1
)

STARTING_DRUG = "Levobupivacaine"
DRUG_DESCRIPTION = df.loc[df["NAME"] == STARTING_DRUG]["DESC"].iloc[0]
DRUG_IMG = df.loc[df["NAME"] == STARTING_DRUG]["IMG_URL"].iloc[0]
FIGURE = create_plot(
    x=df["PKA"],
    y=df["LOGP"],
    z=df["SOL"],
    size=df["MW"],
    color=df["MW"],
    name=df["NAME"],
)

app.layout = html.Div(
    [
        html.Div(
            [html.Img(src=app.get_asset_url("dash-logo.png"))], className="app__banner"
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.H3(
                                    "dash for drug discovery",
                                    className="uppercase title",
                                ),
                                html.Span("Hover ", className="uppercase bold"),
                                html.Span(
                                    "over a drug in the graph to see its structure."
                                ),
                                html.Br(),
                                html.Span("Select ", className="uppercase bold"),
                                html.Span(
                                    "a drug in the dropdown to add it to the drug candidates at the bottom."
                                ),
                            ]
                        )
                    ],
                    className="app__header",
                ),
                html.Div(
                    [
                        dcc.Dropdown(
                            id="chem_dropdown",
                            multi=True,
                            value=[STARTING_DRUG],
                            options=[{"label": i, "value": i} for i in df["NAME"]],
                        )
                    ],
                    className="app__dropdown",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                dcc.RadioItems(
                                    id="charts_radio",
                                    options=[
                                        {"label": "3D Scatter", "value": "scatter3d"},
                                        {"label": "2D Scatter", "value": "scatter"},
                                        {
                                            "label": "2D Histogram",
                                            "value": "histogram2d",
                                        },
                                    ],
                                    labelClassName="radio__labels",
                                    inputClassName="radio__input",
                                    value="scatter3d",
                                    className="radio__group",
                                ),
                                dcc.Graph(
                                    id="clickable-graph",
                                    hoverData={"points": [{"pointNumber": 0}]},
                                    figure=FIGURE,
                                ),
                            ],
                            className="two-thirds column",
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Img(
                                            id="chem_img",
                                            src=DRUG_IMG,
                                            className="chem__img",
                                        )
                                    ],
                                    className="chem__img__container",
                                ),
                                html.Div(
                                    [
                                        html.A(
                                            STARTING_DRUG,
                                            id="chem_name",
                                            href="https://www.drugbank.ca/drugs/DB01002",
                                            target="_blank",
                                        ),
                                        html.P(DRUG_DESCRIPTION, id="chem_desc"),
                                    ],
                                    className="chem__desc__container",
                                ),
                            ],
                            className="one-third column",
                        ),
                    ],
                    className="container card app__content bg-white",
                ),
                html.Div(
                    [
                        html.Table(
                            make_dash_table([STARTING_DRUG], df),
                            id="table-element",
                            className="table__container",
                        )
                    ],
                    className="container bg-white p-0",
                ),
            ],
            className="app__container",
        ),
    ]
)


def df_row_from_hover(hoverData):
    """ Returns row for hover point as a Pandas Series. """

    try:
        point_number = hoverData["points"][0]["pointNumber"]
        molecule_name = str(FIGURE["data"][0]["text"][point_number]).strip()
        return df.loc[df["NAME"] == molecule_name]
    except KeyError as error:
        print(error)
        return pd.Series()


@app.callback(
    Output("clickable-graph", "figure"),
    [Input("chem_dropdown", "value"), Input("charts_radio", "value")],
)
def highlight_molecule(chem_dropdown_values, plot_type):
    """
    Selected chemical dropdown values handler.

    :params chem_dropdown_values: selected dropdown values
    :params plot_type: selected plot graph
    """

    return create_plot(
        x=df["PKA"],
        y=df["LOGP"],
        z=df["SOL"],
        size=df["MW"],
        color=df["MW"],
        name=df["NAME"],
        markers=chem_dropdown_values,
        plot_type=plot_type,
    )


@app.callback(Output("table-element", "children"), [Input("chem_dropdown", "value")])
def update_table(chem_dropdown_value):
    """
    Update the table rows.

    :params chem_dropdown_values: selected dropdown values
    """

    return make_dash_table(chem_dropdown_value, df)


@app.callback(
    [
        Output("chem_name", "children"),
        Output("chem_name", "href"),
        Output("chem_img", "src"),
        Output("chem_desc", "children"),
    ],
    [Input("clickable-graph", "hoverData")],
)
def chem_info_on_hover(hoverData):
    """
    Display chemical information on graph hover.
    Update the image, link, description.

    :params hoverData: data on graph hover
    """

    if hoverData is None:
        raise PreventUpdate

    try:
        row = df_row_from_hover(hoverData)
        if row.empty:
            raise Exception
        return (
            row["NAME"].iloc[0],
            row["PAGE"].iloc[0],
            row["IMG_URL"].iloc[0],
            row["DESC"].iloc[0],
        )

    except Exception as error:
        print(error)
        raise PreventUpdate


if __name__ == "__main__":
    app.run_server(debug=True)
