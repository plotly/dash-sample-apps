# -*- coding: utf-8 -*-
import base64
import io
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import pandas as pd
import pathlib

from dash.dependencies import Input, Output, State
from scipy import stats

group_colors = {"control": "light blue", "reference": "red"}

app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)

PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()
default_study_data = pd.read_csv(DATA_PATH.joinpath("study.csv"))

error_study_data = pd.DataFrame(
    dict(study_id=["error"], group_id=[""], group_type=[""], reading_value=[0])
)

app.layout = html.Div(
    className="",
    children=[
        html.Div(id="error-message"),

        # Top Banner
        html.Div(
            className="study-browser-banner row",
            children=[
                html.H2(
                    className="h2-title",
                    children="ANIMAL STUDY BROWSER"
                ),
                html.Img(
                    className="logo",
                    src=app.get_asset_url("dash-logo.png"),
                ),
            ]
        ),

        # Body of the App
        html.Div(
            className="row",
            children=[
                
                html.Div(
                    className="four columns card",
                    children=[
                        html.Div(
                            className="bg-white",
                            children=[
                                html.H6("Test Articles"),
                                dcc.Dropdown(id="study-dropdown"),
                                html.H6("Choose the type of plot"),
                                dcc.RadioItems(
                                    id="chart-type",
                                    options=[
                                        {"label": "Box Plot", "value": "box"},
                                        {"label": "Violin Plot", "value": "violin"},
                                    ],
                                    value="violin",
                                    labelStyle={
                                        'display': 'inline-block',
                                        'padding':'12px 12px 12px 0px'    
                                    }
                                ),
                                html.H6("CSV File"),
                                dcc.Upload(
                                    id="upload-data",
                                    className="upload",
                                    children=html.Div(["Drag and Drop or ", html.A("Select Files")])
                                )
                            ]
                        )
                    ]
                ),

                # Graph
                html.Div(
                    className="eight columns card",
                    children=[
                        html.Div(
                            className="bg-white",
                            children=[
                                html.H5("Animal data plot"),
                                dcc.Graph(id="plot")
                            ]
                        )
                    ]
                )

            ]
        ),

    ],
)


@app.callback(Output("error-message", "children"), [Input("upload-data", "contents")])
def update_error(contents):
    if contents:
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)
        try:
            study_data = pd.read_csv(io.StringIO(decoded.decode("utf-8")))
        except pd.errors.ParserError:
            return html.Div(
                className="alert",
                children=["That doesn't seem to be a valid csv file!"],
            )
    else:
        study_data = default_study_data

    missing_columns = {
        "group_id",
        "group_type",
        "reading_value",
        "study_id",
    }.difference(study_data.columns)

    if missing_columns:
        return html.Div(
            className="alert", children=["Missing columns: " + str(missing_columns)]
        )

    return None


@app.callback(
    [
        Output("study-dropdown", "options"),
        Output("study-dropdown", "value")
    ],
    [Input("error-message", "children")],
    [State("upload-data", "contents")],
)
def update_dropdown(error_message, contents):
    if error_message:
        study_data = error_study_data
    elif contents:
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)
        study_data = pd.read_csv(io.StringIO(decoded.decode("utf-8")))
    else:
        study_data = default_study_data

    options = []
    if "test_article" in study_data.columns:
        test_articles = study_data.test_article.unique()
        for test_article in test_articles:
            studies = study_data.study_id[
                study_data.test_article == test_article
            ].unique()
            for study in studies:
                options.append(
                    {
                        "label": "{} (study: {})".format(test_article, study),
                        "value": study,
                    }
                )
    else:
        studies = study_data.study_id.unique()
        for study in studies:
            options.append({"label": study, "value": study})

    options.sort(key=lambda item: item["label"])
    value = options[0]["value"]

    return options, value


@app.callback(
    Output("plot", "figure"),
    [Input("chart-type", "value"), Input("study-dropdown", "value")],
    [State("upload-data", "contents"), State("error-message", "children")],
)
def update_output(chart_type, study, contents, error_message):
    if error_message:
        study_data = error_study_data
    elif contents:
        content_type, content_string = contents.split(",")
        decoded = base64.b64decode(content_string)
        study_data = pd.read_csv(io.StringIO(decoded.decode("utf-8")))
        study_data["reading_value"] = pd.to_numeric(
            study_data["reading_value"], errors="coerce"
        )
    else:
        study_data = default_study_data

    if study is None:
        study = study_data.study_id[0]

    study_data = study_data[study_data.study_id == study]
    vehicle_readings = study_data["reading_value"][
        study_data["group_type"] == "control"
    ]
    data_range = study_data["reading_value"].max() - study_data["reading_value"].min()

    test_stats = {}
    box_data = []
    violin_data = []
    for i, group_id in enumerate(study_data.group_id.unique()):
        try:
            group_name = study_data["group_name"][
                study_data.group_id == group_id
            ].values[0]
        except KeyError:
            group_name = group_id

        group_type = study_data["group_type"][study_data.group_id == group_id].values[0]
        y_data = study_data["reading_value"][study_data.group_id == group_id]

        try:
            subject_ids = study_data["subject_id"][study_data.group_id == group_id]
        except KeyError:
            subject_ids = None

        t, p = stats.ttest_ind(vehicle_readings, y_data)
        test_stats[group_id] = {
            "t": t,
            "p": p,
            "pf": "p={:0.3f}".format(p) if p >= 0.001 else "p<0.001",
            "astrix": "***"
            if p <= 0.001
            else "**"
            if p <= 0.01
            else "*"
            if p <= 0.05
            else "",
            "max_y": y_data.max(),
            "index": i,
        }

        box_data.append(
            go.Box(
                y=y_data,
                name=group_name,
                text=subject_ids,
                hoveron="points",
                boxmean=True,
                showlegend=False,
                boxpoints="all",
                pointpos=0,
                line={"color": group_colors.get(group_type, "green")},
            )
        )
        violin_data.append(
            go.Violin(
                y=y_data,
                name=group_name,
                text=subject_ids,
                hoveron="points",
                meanline={"visible": True},
                showlegend=False,
                points="all",
                pointpos=0,
                line={"color": group_colors.get(group_type, "green")},
            )
        )

    chart_data = {"box": box_data, "violin": violin_data}

    if "reading_name" in study_data.columns:
        reading_name = study_data["reading_name"].unique()[0]
    else:
        reading_name = None

    if not vehicle_readings.empty:
        ref_groups = set(
            study_data.group_id[study_data.group_type == "reference"].unique()
        )
        control_groups = set(
            study_data.group_id[study_data.group_type == "control"].unique()
        )
        all_groups = set(study_data.group_id.unique())
        groups_to_annotate = all_groups - ref_groups - control_groups
        annotations = [
            dict(
                x=test_stats[group_id]["index"],
                y=test_stats[group_id]["max_y"]
                + data_range / (4 if chart_type == "violin" else 10),
                text="{}<br>{}".format(
                    test_stats[group_id]["astrix"], test_stats[group_id]["pf"]
                ),
                showarrow=False,
            )
            for group_id in groups_to_annotate
        ]
    else:
        annotations = None

    figure = go.Figure(
        data=chart_data[chart_type],
        layout=go.Layout(
            margin=go.layout.Margin(t=0, r=50, b=50, l=50),
            yaxis=dict(
                title=dict(text=reading_name)
            ),
            annotations=annotations,
        ),
    )

    return figure


if __name__ == "__main__":
    app.run_server(debug=True)
