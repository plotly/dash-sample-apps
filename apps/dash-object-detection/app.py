from textwrap import dedent
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_player as player
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
import pathlib

FRAMERATE = 24.0

app = dash.Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)
server = app.server
app.config.suppress_callback_exceptions = True

BASE_PATH = pathlib.Path(__file__).parent.resolve()
DATA_PATH = BASE_PATH.joinpath("data").resolve()


def load_data(path):
    """Load data about a specific footage (given by the path). It returns a dictionary of useful variables such as
    the dataframe containing all the detection and bounds localization, the number of classes inside that footage,
    the matrix of all the classes in string, the given class with padding, and the root of the number of classes,
    rounded."""

    # Load the dataframe containing all the processed object detections inside the video
    video_info_df = pd.read_csv(DATA_PATH.joinpath(path))

    # The list of classes, and the number of classes
    classes_list = video_info_df["class_str"].value_counts().index.tolist()
    n_classes = len(classes_list)

    # Gets the smallest value needed to add to the end of the classes list to get a square matrix
    root_round = np.ceil(np.sqrt(len(classes_list)))
    total_size = root_round ** 2
    padding_value = int(total_size - n_classes)
    classes_padded = np.pad(classes_list, (0, padding_value), mode="constant")

    # The padded matrix containing all the classes inside a matrix
    classes_matrix = np.reshape(classes_padded, (int(root_round), int(root_round)))

    # Flip it for better looks
    classes_matrix = np.flip(classes_matrix, axis=0)

    data_dict = {
        "video_info_df": video_info_df,
        "n_classes": n_classes,
        "classes_matrix": classes_matrix,
        "classes_padded": classes_padded,
        "root_round": root_round,
    }

    if True:
        print(f"{path} loaded.")

    return data_dict


def markdown_popup():
    return html.Div(
        id="markdown",
        className="modal",
        style={"display": "none"},
        children=(
            html.Div(
                className="markdown-container",
                children=[
                    html.Div(
                        className="close-container",
                        children=html.Button(
                            "Close",
                            id="markdown_close",
                            n_clicks=0,
                            className="closeButton",
                        ),
                    ),
                    html.Div(
                        className="markdown-text",
                        children=[
                            dcc.Markdown(
                                children=dedent(
                                    """
                                ##### What am I looking at?
                                
                                This app enhances visualization of objects detected using state-of-the-art Mobile Vision Neural Networks.
                                Most user generated videos are dynamic and fast-paced, which might be hard to interpret. A confidence
                                heatmap stays consistent through the video and intuitively displays the model predictions. The pie chart
                                lets you interpret how the object classes are divided, which is useful when analyzing videos with numerous
                                and differing objects.

                                ##### More about this Dash app
                                
                                The purpose of this demo is to explore alternative visualization methods for object detection. Therefore,
                                the visualizations, predictions and videos are not generated in real time, but done beforehand. To read
                                more about it, please visit the [project repo](https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-object-detection).

                                """
                                )
                            )
                        ],
                    ),
                ],
            )
        ),
    )


# Main App
app.layout = html.Div(
    children=[
        dcc.Interval(id="interval-updating-graphs", interval=1000, n_intervals=0),
        html.Div(id="top-bar", className="row"),
        html.Div(
            className="container",
            children=[
                html.Div(
                    id="left-side-column",
                    className="eight columns",
                    children=[
                        html.Img(
                            id="logo-mobile", src=app.get_asset_url("dash-logo.png")
                        ),
                        html.Div(
                            id="header-section",
                            children=[
                                html.H4("Object Detection Explorer"),
                                html.P(
                                    "To get started, select the footage you want to view, and choose the display mode (with or without "
                                    "bounding boxes). Then, you can start playing the video, and the result of objects detected "
                                    "will be displayed in accordance to the current video-playing time."
                                ),
                                html.Button(
                                    "Learn More", id="learn-more-button", n_clicks=0
                                ),
                            ],
                        ),
                        html.Div(
                            className="video-outer-container",
                            children=html.Div(
                                className="video-container",
                                children=player.DashPlayer(
                                    id="video-display",
                                    url="https://www.youtube.com/watch?v=gPtn6hD7o8g",
                                    controls=True,
                                    playing=False,
                                    volume=1,
                                    width="100%",
                                    height="100%",
                                ),
                            ),
                        ),
                        html.Div(
                            className="control-section",
                            children=[
                                html.Div(
                                    className="control-element",
                                    children=[
                                        html.Div(
                                            children=["Minimum Confidence Threshold:"]
                                        ),
                                        html.Div(
                                            dcc.Slider(
                                                id="slider-minimum-confidence-threshold",
                                                min=20,
                                                max=80,
                                                marks={
                                                    i: f"{i}%"
                                                    for i in range(20, 81, 10)
                                                },
                                                value=30,
                                                updatemode="drag",
                                            )
                                        ),
                                    ],
                                ),
                                html.Div(
                                    className="control-element",
                                    children=[
                                        html.Div(children=["Footage Selection:"]),
                                        dcc.Dropdown(
                                            id="dropdown-footage-selection",
                                            options=[
                                                {
                                                    "label": "Drone recording of canal festival",
                                                    "value": "DroneCanalFestival",
                                                },
                                                {
                                                    "label": "Drone recording of car festival",
                                                    "value": "car_show_drone",
                                                },
                                                {
                                                    "label": "Drone recording of car festival #2",
                                                    "value": "DroneCarFestival2",
                                                },
                                                {
                                                    "label": "Drone recording of a farm",
                                                    "value": "FarmDrone",
                                                },
                                                {
                                                    "label": "Lion fighting Zebras",
                                                    "value": "zebra",
                                                },
                                                {
                                                    "label": "Man caught by a CCTV",
                                                    "value": "ManCCTV",
                                                },
                                                {
                                                    "label": "Man driving expensive car",
                                                    "value": "car_footage",
                                                },
                                                {
                                                    "label": "Restaurant Robbery",
                                                    "value": "RestaurantHoldup",
                                                },
                                            ],
                                            value="car_show_drone",
                                            clearable=False,
                                        ),
                                    ],
                                ),
                                html.Div(
                                    className="control-element",
                                    children=[
                                        html.Div(children=["Video Display Mode:"]),
                                        dcc.Dropdown(
                                            id="dropdown-video-display-mode",
                                            options=[
                                                {
                                                    "label": "Regular Display",
                                                    "value": "regular",
                                                },
                                                {
                                                    "label": "Display with Bounding Boxes",
                                                    "value": "bounding_box",
                                                },
                                            ],
                                            value="bounding_box",
                                            searchable=False,
                                            clearable=False,
                                        ),
                                    ],
                                ),
                                html.Div(
                                    className="control-element",
                                    children=[
                                        html.Div(children=["Graph View Mode:"]),
                                        dcc.Dropdown(
                                            id="dropdown-graph-view-mode",
                                            options=[
                                                {
                                                    "label": "Visual Mode",
                                                    "value": "visual",
                                                },
                                                {
                                                    "label": "Detection Mode",
                                                    "value": "detection",
                                                },
                                            ],
                                            value="visual",
                                            searchable=False,
                                            clearable=False,
                                        ),
                                    ],
                                ),
                            ],
                        ),
                    ],
                ),
                html.Div(
                    id="right-side-column",
                    className="four columns",
                    children=[
                        html.Div(
                            className="img-container",
                            children=html.Img(
                                id="logo-web", src=app.get_asset_url("dash-logo.png")
                            ),
                        ),
                        html.Div(id="div-visual-mode"),
                        html.Div(id="div-detection-mode"),
                    ],
                ),
            ],
        ),
        markdown_popup(),
    ]
)


# Data Loading
@app.server.before_first_request
def load_all_footage():
    global data_dict, url_dict

    # Load the dictionary containing all the variables needed for analysis
    data_dict = {
        "james_bond": load_data("james_bond_object_data.csv"),
        "zebra": load_data("Zebra_object_data.csv"),
        "car_show_drone": load_data("CarShowDrone_object_data.csv"),
        "car_footage": load_data("CarFootage_object_data.csv"),
        "DroneCanalFestival": load_data("DroneCanalFestivalDetectionData.csv"),
        "DroneCarFestival2": load_data("DroneCarFestival2DetectionData.csv"),
        "FarmDrone": load_data("FarmDroneDetectionData.csv"),
        "ManCCTV": load_data("ManCCTVDetectionData.csv"),
        "RestaurantHoldup": load_data("RestaurantHoldupDetectionData.csv"),
    }

    url_dict = {
        "regular": {
            "james_bond": "https://www.youtube.com/watch?v=g9S5GndUhko",
            "zebra": "https://www.youtube.com/watch?v=TVvtD3AVt10",
            "car_show_drone": "https://www.youtube.com/watch?v=gPtn6hD7o8g",
            "car_footage": "https://www.youtube.com/watch?v=qX3bDxHuq6I",
            "DroneCanalFestival": "https://youtu.be/0oucTt2OW7M",
            "DroneCarFestival2": "https://youtu.be/vhJ7MHsJvwY",
            "FarmDrone": "https://youtu.be/aXfKuaP8v_A",
            "ManCCTV": "https://youtu.be/BYZORBIxgbc",
            "RestaurantHoldup": "https://youtu.be/WDin4qqgpac",
        },
        "bounding_box": {
            "james_bond": "https://www.youtube.com/watch?v=g9S5GndUhko",
            "zebra": "https://www.youtube.com/watch?v=G2pbZgyWQ5E",
            "car_show_drone": "https://www.youtube.com/watch?v=9F5FdcVmLOY",
            "car_footage": "https://www.youtube.com/watch?v=EhnNosq1Lrc",
            "DroneCanalFestival": "https://youtu.be/6ZZmsnwk2HQ",
            "DroneCarFestival2": "https://youtu.be/2Gr4RQ-JHIs",
            "FarmDrone": "https://youtu.be/pvvW5yZlpyc",
            "ManCCTV": "https://youtu.be/1oMrHLrtOZw",
            "RestaurantHoldup": "https://youtu.be/HOIKOwixYEY",
        },
    }


# Footage Selection
@app.callback(
    Output("video-display", "url"),
    [
        Input("dropdown-footage-selection", "value"),
        Input("dropdown-video-display-mode", "value"),
    ],
)
def select_footage(footage, display_mode):
    # Find desired footage and update player video
    url = url_dict[display_mode][footage]
    return url


# Learn more popup
@app.callback(
    Output("markdown", "style"),
    [Input("learn-more-button", "n_clicks"), Input("markdown_close", "n_clicks")],
)
def update_click_output(button_click, close_click):
    ctx = dash.callback_context
    prop_id = ""
    if ctx.triggered:
        prop_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if prop_id == "learn-more-button":
        return {"display": "block"}
    else:
        return {"display": "none"}


@app.callback(
    Output("div-visual-mode", "children"), [Input("dropdown-graph-view-mode", "value")]
)
def update_output(dropdown_value):
    if dropdown_value == "visual":
        return [
            html.Div(
                children=[
                    html.P(
                        children="Confidence Level of Object Presence",
                        className="plot-title",
                    ),
                    dcc.Graph(
                        id="heatmap-confidence",
                        figure=dict(
                            layout={
                                "paper_bgcolor": "rgb(249,249,249)",
                                "plot_bgcolor": "rgb(249,249,249)",
                            }
                        ),
                    ),
                    html.P(children="Object Count", className="plot-title"),
                    dcc.Graph(
                        id="pie-object-count",
                        figure=dict(
                            layout={
                                "paper_bgcolor": "rgb(249,249,249)",
                                "plot_bgcolor": "rgb(249,249,249)",
                            }
                        ),
                    ),
                ]
            )
        ]
    return []


@app.callback(
    Output("div-detection-mode", "children"),
    [Input("dropdown-graph-view-mode", "value")],
)
def update_detection_mode(value):
    if value == "detection":
        return [
            html.Div(
                children=[
                    html.P(
                        children="Detection Score of Most Probable Objects",
                        className="plot-title",
                    ),
                    dcc.Graph(id="bar-score-graph"),
                ]
            )
        ]
    return []


# Updating Figures
@app.callback(
    Output("bar-score-graph", "figure"),
    [Input("interval-updating-graphs", "n_intervals")],
    [
        State("video-display", "currentTime"),
        State("dropdown-footage-selection", "value"),
        State("slider-minimum-confidence-threshold", "value"),
    ],
)
def update_score_bar(n, current_time, footage, threshold):
    layout = go.Layout(
        showlegend=False,
        paper_bgcolor="rgb(249,249,249)",
        plot_bgcolor="rgb(249,249,249)",
        xaxis={"automargin": True},
        yaxis={"title": "Score", "automargin": True, "range": [0, 1]},
    )

    if current_time is not None:
        current_frame = round(current_time * FRAMERATE)

        if n > 0 and current_frame > 0:
            video_info_df = data_dict[footage]["video_info_df"]

            # Select the subset of the dataset that correspond to the current frame
            frame_df = video_info_df[video_info_df["frame"] == current_frame]

            # Select only the frames above the threshold
            threshold_dec = threshold / 100  # Threshold in decimal
            frame_df = frame_df[frame_df["score"] > threshold_dec]

            # Select up to 8 frames with the highest scores
            frame_df = frame_df[: min(8, frame_df.shape[0])]

            # Add count to object names (e.g. person --> person 1, person --> person 2)
            objects = frame_df["class_str"].tolist()
            object_count_dict = {
                x: 0 for x in set(objects)
            }  # Keeps count of the objects
            objects_wc = []  # Object renamed with counts
            for object in objects:
                object_count_dict[object] += 1  # Increment count
                objects_wc.append(f"{object} {object_count_dict[object]}")

            colors = list("rgb(250,79,86)" for i in range(len(objects_wc)))

            # Add text information
            y_text = [
                f"{round(value * 100)}% confidence"
                for value in frame_df["score"].tolist()
            ]

            figure = go.Figure(
                {
                    "data": [
                        {
                            "hoverinfo": "x+text",
                            "name": "Detection Scores",
                            "text": y_text,
                            "type": "bar",
                            "x": objects_wc,
                            "marker": {"color": colors},
                            "y": frame_df["score"].tolist(),
                        }
                    ],
                    "layout": {
                        "showlegend": False,
                        "autosize": False,
                        "paper_bgcolor": "rgb(249,249,249)",
                        "plot_bgcolor": "rgb(249,249,249)",
                        "xaxis": {"automargin": True, "tickangle": -45},
                        "yaxis": {
                            "automargin": True,
                            "range": [0, 1],
                            "title": {"text": "Score"},
                        },
                    },
                }
            )
            return figure

    return go.Figure(data=[go.Bar()], layout=layout)  # Returns empty bar


@app.callback(
    Output("pie-object-count", "figure"),
    [Input("interval-updating-graphs", "n_intervals")],
    [
        State("video-display", "currentTime"),
        State("dropdown-footage-selection", "value"),
        State("slider-minimum-confidence-threshold", "value"),
    ],
)
def update_object_count_pie(n, current_time, footage, threshold):
    layout = go.Layout(
        showlegend=True,
        paper_bgcolor="rgb(249,249,249)",
        plot_bgcolor="rgb(249,249,249)",
        autosize=False,
        margin=dict(l=10, r=10, t=15, b=15),
    )

    if current_time is not None:
        current_frame = round(current_time * FRAMERATE)

        if n > 0 and current_frame > 0:
            video_info_df = data_dict[footage]["video_info_df"]

            # Select the subset of the dataset that correspond to the current frame
            frame_df = video_info_df[video_info_df["frame"] == current_frame]

            # Select only the frames above the threshold
            threshold_dec = threshold / 100  # Threshold in decimal
            frame_df = frame_df[frame_df["score"] > threshold_dec]

            # Get the count of each object class
            class_counts = frame_df["class_str"].value_counts()

            classes = class_counts.index.tolist()  # List of each class
            counts = class_counts.tolist()  # List of each count

            text = [f"{count} detected" for count in counts]

            # Set colorscale to piechart
            colorscale = [
                "#fa4f56",
                "#fe6767",
                "#ff7c79",
                "#ff908b",
                "#ffa39d",
                "#ffb6b0",
                "#ffc8c3",
                "#ffdbd7",
                "#ffedeb",
                "#ffffff",
            ]

            pie = go.Pie(
                labels=classes,
                values=counts,
                text=text,
                hoverinfo="text+percent",
                textinfo="label+percent",
                marker={"colors": colorscale[: len(classes)]},
            )
            return go.Figure(data=[pie], layout=layout)

    return go.Figure(data=[go.Pie()], layout=layout)  # Returns empty pie chart


@app.callback(
    Output("heatmap-confidence", "figure"),
    [Input("interval-updating-graphs", "n_intervals")],
    [
        State("video-display", "currentTime"),
        State("dropdown-footage-selection", "value"),
        State("slider-minimum-confidence-threshold", "value"),
    ],
)
def update_heatmap_confidence(n, current_time, footage, threshold):
    layout = go.Layout(
        showlegend=False,
        paper_bgcolor="rgb(249,249,249)",
        plot_bgcolor="rgb(249,249,249)",
        autosize=False,
        margin=dict(l=10, r=10, b=20, t=20, pad=4),
    )

    if current_time is not None:
        current_frame = round(current_time * FRAMERATE)

        if n > 0 and current_frame > 0:
            # Load variables from the data dictionary
            video_info_df = data_dict[footage]["video_info_df"]
            classes_padded = data_dict[footage]["classes_padded"]
            root_round = data_dict[footage]["root_round"]
            classes_matrix = data_dict[footage]["classes_matrix"]

            # Select the subset of the dataset that correspond to the current frame
            frame_df = video_info_df[video_info_df["frame"] == current_frame]

            # Select only the frames above the threshold
            threshold_dec = threshold / 100
            frame_df = frame_df[frame_df["score"] > threshold_dec]

            # Remove duplicate, keep the top result
            frame_no_dup = frame_df[["class_str", "score"]].drop_duplicates("class_str")
            frame_no_dup.set_index("class_str", inplace=True)

            # The list of scores
            score_list = []
            for el in classes_padded:
                if el in frame_no_dup.index.values:
                    score_list.append(frame_no_dup.loc[el][0])
                else:
                    score_list.append(0)

            # Generate the score matrix, and flip it for visual
            score_matrix = np.reshape(score_list, (-1, int(root_round)))
            score_matrix = np.flip(score_matrix, axis=0)

            # We set the color scale to white if there's nothing in the frame_no_dup
            if frame_no_dup.shape != (0, 1):
                colorscale = [[0, "#f9f9f9"], [1, "#fa4f56"]]
            else:
                colorscale = [[0, "#f9f9f9"], [1, "#f9f9f9"]]

            hover_text = [f"{score * 100:.2f}% confidence" for score in score_list]
            hover_text = np.reshape(hover_text, (-1, int(root_round)))
            hover_text = np.flip(hover_text, axis=0)

            # Add linebreak for multi-word annotation
            classes_matrix = classes_matrix.astype(dtype="|U40")

            for index, row in enumerate(classes_matrix):
                row = list(map(lambda x: "<br>".join(x.split()), row))
                classes_matrix[index] = row

            # Set up annotation text
            annotation = []
            for y_cord in range(int(root_round)):
                for x_cord in range(int(root_round)):
                    annotation_dict = dict(
                        showarrow=False,
                        text=classes_matrix[y_cord][x_cord],
                        xref="x",
                        yref="y",
                        x=x_cord,
                        y=y_cord,
                    )
                    if score_matrix[y_cord][x_cord] > 0:
                        annotation_dict["font"] = {"color": "#F9F9F9", "size": "11"}
                    else:
                        annotation_dict["font"] = {"color": "#606060", "size": "11"}
                    annotation.append(annotation_dict)

            # Generate heatmap figure

            figure = {
                "data": [
                    {
                        "colorscale": colorscale,
                        "showscale": False,
                        "hoverinfo": "text",
                        "text": hover_text,
                        "type": "heatmap",
                        "zmin": 0,
                        "zmax": 1,
                        "xgap": 1,
                        "ygap": 1,
                        "z": score_matrix,
                    }
                ],
                "layout": {
                    "showlegend": False,
                    "autosize": False,
                    "paper_bgcolor": "rgb(249,249,249)",
                    "plot_bgcolor": "rgb(249,249,249)",
                    "margin": {"l": 10, "r": 10, "b": 20, "t": 20, "pad": 2},
                    "annotations": annotation,
                    "xaxis": {
                        "showticklabels": False,
                        "showgrid": False,
                        "showline": False,
                        "zeroline": False,
                        "side": "top",
                        "ticks": "",
                    },
                    "yaxis": {
                        "showticklabels": False,
                        "showgrid": False,
                        "showline": False,
                        "zeroline": False,
                        "side": "left",
                        "ticks": "",
                    },
                },
            }

            return figure

    # Returns empty figure
    return go.Figure(data=[go.Pie()], layout=layout)


# Running the server
if __name__ == "__main__":
    app.run_server(debug=True, port=8053)
