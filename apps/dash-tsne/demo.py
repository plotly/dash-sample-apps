import base64
import io
import pathlib

import numpy as np
import dash
import dash_core_components as dcc
import dash_html_components as html
from PIL import Image
from io import BytesIO
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate

import pandas as pd
import plotly.graph_objs as go
import scipy.spatial.distance as spatial_distance

# get relative data folder
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("data").resolve()

data_dict = {
    "mnist_3000": pd.read_csv(DATA_PATH.joinpath("mnist_3000_input.csv")),
    "wikipedia_3000": pd.read_csv(DATA_PATH.joinpath("wikipedia_3000.csv")),
    "twitter_3000": pd.read_csv(
        DATA_PATH.joinpath("twitter_3000.csv"), encoding="ISO-8859-1"
    ),
}

# Import datasets here for running the Local version
IMAGE_DATASETS = "mnist_3000"
WORD_EMBEDDINGS = ("wikipedia_3000", "twitter_3000")


with open(PATH.joinpath("demo_intro.md"), "r") as file:
    demo_intro_md = file.read()

with open(PATH.joinpath("demo_description.md"), "r") as file:
    demo_description_md = file.read()


def numpy_to_b64(array, scalar=True):
    # Convert from 0-1 to 0-255
    if scalar:
        array = np.uint8(255 * array)

    im_pil = Image.fromarray(array)
    buff = BytesIO()
    im_pil.save(buff, format="png")
    im_b64 = base64.b64encode(buff.getvalue()).decode("utf-8")

    return im_b64


# Methods for creating components in the layout code
def Card(children, **kwargs):
    return html.Section(children, className="card-style")


def NamedSlider(name, short, min, max, step, val, marks=None):
    if marks:
        step = None
    else:
        marks = {i: i for i in range(min, max + 1, step)}

    return html.Div(
        style={"margin": "25px 5px 30px 0px"},
        children=[
            f"{name}:",
            html.Div(
                style={"margin-left": "5px"},
                children=[
                    dcc.Slider(
                        id=f"slider-{short}",
                        min=min,
                        max=max,
                        marks=marks,
                        step=step,
                        value=val,
                    )
                ],
            ),
        ],
    )


def NamedInlineRadioItems(name, short, options, val, **kwargs):
    return html.Div(
        id=f"div-{short}",
        style={"display": "inline-block"},
        children=[
            f"{name}:",
            dcc.RadioItems(
                id=f"radio-{short}",
                options=options,
                value=val,
                labelStyle={"display": "inline-block", "margin-right": "7px"},
                style={"display": "inline-block", "margin-left": "7px"},
            ),
        ],
    )


def create_layout(app):
    # Actual layout of the app
    return html.Div(
        className="row",
        style={"max-width": "100%", "font-size": "1.5rem", "padding": "0px 0px"},
        children=[
            # Header
            html.Div(
                className="row header",
                id="app-header",
                style={"background-color": "#f9f9f9"},
                children=[
                    html.Div(
                        [
                            html.Img(
                                src=app.get_asset_url("dash-logo.png"),
                                className="logo",
                                id="plotly-image",
                            )
                        ],
                        className="three columns header_img",
                    ),
                    html.Div(
                        [
                            html.H3(
                                "t-SNE Explorer",
                                className="header_title",
                                id="app-title",
                            )
                        ],
                        className="nine columns header_title_container",
                    ),
                ],
            ),
            # Demo Description
            html.Div(
                className="row background",
                id="demo-explanation",
                style={"padding": "50px 45px"},
                children=[
                    html.Div(
                        id="description-text", children=dcc.Markdown(demo_intro_md)
                    ),
                    html.Div(
                        html.Button(id="learn-more-button", children=["Learn More"])
                    ),
                ],
            ),
            # Body
            html.Div(
                className="row background",
                style={"padding": "10px"},
                children=[
                    html.Div(
                        className="three columns",
                        children=[
                            Card(
                                [
                                    dcc.Dropdown(
                                        id="dropdown-dataset",
                                        searchable=False,
                                        clearable=False,
                                        options=[
                                            {
                                                "label": "MNIST Digits",
                                                "value": "mnist_3000",
                                            },
                                            {
                                                "label": "Twitter (GloVe)",
                                                "value": "twitter_3000",
                                            },
                                            {
                                                "label": "Wikipedia (GloVe)",
                                                "value": "wikipedia_3000",
                                            },
                                        ],
                                        placeholder="Select a dataset",
                                        value="mnist_3000",
                                    ),
                                    NamedSlider(
                                        name="Number Of Iterations",
                                        short="iterations",
                                        min=250,
                                        max=1000,
                                        step=None,
                                        val=500,
                                        marks={
                                            i: str(i) for i in [250, 500, 750, 1000]
                                        },
                                    ),
                                    NamedSlider(
                                        name="Perplexity",
                                        short="perplexity",
                                        min=3,
                                        max=100,
                                        step=None,
                                        val=30,
                                        marks={i: str(i) for i in [3, 10, 30, 50, 100]},
                                    ),
                                    NamedSlider(
                                        name="Initial PCA Dimensions",
                                        short="pca-dimension",
                                        min=25,
                                        max=100,
                                        step=None,
                                        val=50,
                                        marks={i: str(i) for i in [25, 50, 100]},
                                    ),
                                    NamedSlider(
                                        name="Learning Rate",
                                        short="learning-rate",
                                        min=10,
                                        max=200,
                                        step=None,
                                        val=100,
                                        marks={i: str(i) for i in [10, 50, 100, 200]},
                                    ),
                                    html.Div(
                                        id="div-wordemb-controls",
                                        style={"display": "none"},
                                        children=[
                                            NamedInlineRadioItems(
                                                name="Display Mode",
                                                short="wordemb-display-mode",
                                                options=[
                                                    {
                                                        "label": " Regular",
                                                        "value": "regular",
                                                    },
                                                    {
                                                        "label": " Top-100 Neighbors",
                                                        "value": "neighbors",
                                                    },
                                                ],
                                                val="regular",
                                            ),
                                            dcc.Dropdown(
                                                id="dropdown-word-selected",
                                                placeholder="Select word to display its neighbors",
                                                style={"background-color": "#f2f3f4"},
                                            ),
                                        ],
                                    ),
                                ]
                            )
                        ],
                    ),
                    html.Div(
                        className="six columns",
                        children=[
                            dcc.Graph(id="graph-3d-plot-tsne", style={"height": "98vh"})
                        ],
                    ),
                    html.Div(
                        className="three columns",
                        id="euclidean-distance",
                        children=[
                            Card(
                                style={"padding": "5px"},
                                children=[
                                    html.Div(
                                        id="div-plot-click-message",
                                        style={
                                            "text-align": "center",
                                            "margin-bottom": "7px",
                                            "font-weight": "bold",
                                        },
                                    ),
                                    html.Div(id="div-plot-click-image"),
                                    html.Div(id="div-plot-click-wordemb"),
                                ],
                            )
                        ],
                    ),
                ],
            ),
        ],
    )


def demo_callbacks(app):
    def generate_figure_image(groups, layout):
        data = []

        for idx, val in groups:
            scatter = go.Scatter3d(
                name=idx,
                x=val["x"],
                y=val["y"],
                z=val["z"],
                text=[idx for _ in range(val["x"].shape[0])],
                textposition="top center",
                mode="markers",
                marker=dict(size=3, symbol="circle"),
            )
            data.append(scatter)

        figure = go.Figure(data=data, layout=layout)

        return figure

    # Scatter Plot of the t-SNE datasets
    def generate_figure_word_vec(
        embedding_df, layout, wordemb_display_mode, selected_word, dataset
    ):

        try:
            # Regular displays the full scatter plot with only circles
            if wordemb_display_mode == "regular":
                plot_mode = "markers"

            # Nearest Neighbors displays only the 200 nearest neighbors of the selected_word, in text rather than circles
            elif wordemb_display_mode == "neighbors":
                if not selected_word:
                    return go.Figure()

                plot_mode = "text"

                # Get the nearest neighbors indices using Euclidean distance
                vector = data_dict[dataset].set_index("0")
                selected_vec = vector.loc[selected_word]

                def compare_pd(vector):
                    return spatial_distance.euclidean(vector, selected_vec)

                # vector.apply takes compare_pd function as the first argument
                distance_map = vector.apply(compare_pd, axis=1)
                neighbors_idx = distance_map.sort_values()[:100].index

                # Select those neighbors from the embedding_df
                embedding_df = embedding_df.loc[neighbors_idx]

            scatter = go.Scatter3d(
                name=str(embedding_df.index),
                x=embedding_df["x"],
                y=embedding_df["y"],
                z=embedding_df["z"],
                text=embedding_df.index,
                textposition="middle center",
                showlegend=False,
                mode=plot_mode,
                marker=dict(size=3, color="#3266c1", symbol="circle"),
            )

            figure = go.Figure(data=[scatter], layout=layout)

            return figure
        except KeyError as error:
            print(error)
            raise PreventUpdate

    # Callback function for the learn-more button
    @app.callback(
        [
            Output("description-text", "children"),
            Output("learn-more-button", "children"),
        ],
        [Input("learn-more-button", "n_clicks")],
    )
    def learn_more(n_clicks):
        # If clicked odd times, the instructions will show; else (even times), only the header will show
        if n_clicks is None:
            n_clicks = 0
        if (n_clicks % 2) == 1:
            n_clicks += 1
            return (
                html.Div(
                    style={"padding-right": "15%"},
                    children=[dcc.Markdown(demo_description_md)],
                ),
                "Close",
            )
        else:
            n_clicks += 1
            return (
                html.Div(
                    style={"padding-right": "15%"},
                    children=[dcc.Markdown(demo_intro_md)],
                ),
                "Learn More",
            )

    @app.callback(
        Output("div-wordemb-controls", "style"), [Input("dropdown-dataset", "value")]
    )
    def show_wordemb_controls(dataset):
        if dataset in WORD_EMBEDDINGS:
            return None
        else:
            return {"display": "none"}

    @app.callback(
        Output("dropdown-word-selected", "disabled"),
        [Input("radio-wordemb-display-mode", "value")],
    )
    def disable_word_selection(mode):
        return not mode == "neighbors"

    @app.callback(
        Output("dropdown-word-selected", "options"),
        [Input("dropdown-dataset", "value")],
    )
    def fill_dropdown_word_selection_options(dataset):
        if dataset in WORD_EMBEDDINGS:
            return [
                {"label": i, "value": i} for i in data_dict[dataset].iloc[:, 0].tolist()
            ]
        else:
            return []

    @app.callback(
        Output("graph-3d-plot-tsne", "figure"),
        [
            Input("dropdown-dataset", "value"),
            Input("slider-iterations", "value"),
            Input("slider-perplexity", "value"),
            Input("slider-pca-dimension", "value"),
            Input("slider-learning-rate", "value"),
            Input("dropdown-word-selected", "value"),
            Input("radio-wordemb-display-mode", "value"),
        ],
    )
    def display_3d_scatter_plot(
        dataset,
        iterations,
        perplexity,
        pca_dim,
        learning_rate,
        selected_word,
        wordemb_display_mode,
    ):
        if dataset:
            path = f"demo_embeddings/{dataset}/iterations_{iterations}/perplexity_{perplexity}/pca_{pca_dim}/learning_rate_{learning_rate}"

            try:

                data_url = [
                    "demo_embeddings",
                    str(dataset),
                    "iterations_" + str(iterations),
                    "perplexity_" + str(perplexity),
                    "pca_" + str(pca_dim),
                    "learning_rate_" + str(learning_rate),
                    "data.csv",
                ]
                full_path = PATH.joinpath(*data_url)
                embedding_df = pd.read_csv(
                    full_path, index_col=0, encoding="ISO-8859-1"
                )

            except FileNotFoundError as error:
                print(
                    error,
                    "\nThe dataset was not found. Please generate it using generate_demo_embeddings.py",
                )
                return go.Figure()

            # Plot layout
            axes = dict(title="", showgrid=True, zeroline=False, showticklabels=False)

            layout = go.Layout(
                margin=dict(l=0, r=0, b=0, t=0),
                scene=dict(xaxis=axes, yaxis=axes, zaxis=axes),
            )

            # For Image datasets
            if dataset in IMAGE_DATASETS:
                embedding_df["label"] = embedding_df.index

                groups = embedding_df.groupby("label")
                figure = generate_figure_image(groups, layout)

            # Everything else is word embeddings
            elif dataset in WORD_EMBEDDINGS:
                figure = generate_figure_word_vec(
                    embedding_df=embedding_df,
                    layout=layout,
                    wordemb_display_mode=wordemb_display_mode,
                    selected_word=selected_word,
                    dataset=dataset,
                )

            else:
                figure = go.Figure()

            return figure

    @app.callback(
        Output("div-plot-click-image", "children"),
        [
            Input("graph-3d-plot-tsne", "clickData"),
            Input("dropdown-dataset", "value"),
            Input("slider-iterations", "value"),
            Input("slider-perplexity", "value"),
            Input("slider-pca-dimension", "value"),
            Input("slider-learning-rate", "value"),
        ],
    )
    def display_click_image(
        clickData, dataset, iterations, perplexity, pca_dim, learning_rate
    ):
        if dataset in IMAGE_DATASETS and clickData:
            # Load the same dataset as the one displayed

            try:
                data_url = [
                    "demo_embeddings",
                    str(dataset),
                    "iterations_" + str(iterations),
                    "perplexity_" + str(perplexity),
                    "pca_" + str(pca_dim),
                    "learning_rate_" + str(learning_rate),
                    "data.csv",
                ]

                full_path = PATH.joinpath(*data_url)
                embedding_df = pd.read_csv(full_path, encoding="ISO-8859-1")

            except FileNotFoundError as error:
                print(
                    error,
                    "\nThe dataset was not found. Please generate it using generate_demo_embeddings.py",
                )
                return

            # Convert the point clicked into float64 numpy array
            click_point_np = np.array(
                [clickData["points"][0][i] for i in ["x", "y", "z"]]
            ).astype(np.float64)
            # Create a boolean mask of the point clicked, truth value exists at only one row
            bool_mask_click = (
                embedding_df.loc[:, "x":"z"].eq(click_point_np).all(axis=1)
            )
            # Retrieve the index of the point clicked, given it is present in the set
            if bool_mask_click.any():
                clicked_idx = embedding_df[bool_mask_click].index[0]

                # Retrieve the image corresponding to the index
                image_vector = data_dict[dataset].iloc[clicked_idx]
                if dataset == "cifar_gray_3000":
                    image_np = image_vector.values.reshape(32, 32).astype(np.float64)
                else:
                    image_np = image_vector.values.reshape(28, 28).astype(np.float64)

                # Encode image into base 64
                image_b64 = numpy_to_b64(image_np)

                return html.Img(
                    src="data:image/png;base64, " + image_b64,
                    style={"height": "25vh", "display": "block", "margin": "auto"},
                )
        return None

    @app.callback(
        Output("div-plot-click-wordemb", "children"),
        [Input("graph-3d-plot-tsne", "clickData"), Input("dropdown-dataset", "value")],
    )
    def display_click_word_neighbors(clickData, dataset):
        if dataset in WORD_EMBEDDINGS and clickData:
            selected_word = clickData["points"][0]["text"]

            try:
                # Get the nearest neighbors indices using Euclidean distance
                vector = data_dict[dataset].set_index("0")
                selected_vec = vector.loc[selected_word]

                def compare_pd(vector):
                    return spatial_distance.euclidean(vector, selected_vec)

                # vector.apply takes compare_pd function as the first argument
                distance_map = vector.apply(compare_pd, axis=1)
                nearest_neighbors = distance_map.sort_values()[1:6]

                trace = go.Bar(
                    x=nearest_neighbors.values,
                    y=nearest_neighbors.index,
                    width=0.5,
                    orientation="h",
                    marker=dict(color="rgb(50, 102, 193)"),
                )

                layout = go.Layout(
                    title=f'5 nearest neighbors of "{selected_word}"',
                    xaxis=dict(title="Euclidean Distance"),
                    margin=go.layout.Margin(l=60, r=60, t=35, b=35),
                )

                fig = go.Figure(data=[trace], layout=layout)

                return dcc.Graph(
                    id="graph-bar-nearest-neighbors-word",
                    figure=fig,
                    style={"height": "25vh"},
                    config={"displayModeBar": False},
                )
            except KeyError as error:
                raise PreventUpdate
        return None

    @app.callback(
        Output("div-plot-click-message", "children"),
        [Input("graph-3d-plot-tsne", "clickData"), Input("dropdown-dataset", "value")],
    )
    def display_click_message(clickData, dataset):
        # Displays message shown when a point in the graph is clicked, depending whether it's an image or word
        if dataset in IMAGE_DATASETS:
            if clickData:
                return "Image Selected"
            else:
                return "Click a data point on the scatter plot to display its corresponding image."

        elif dataset in WORD_EMBEDDINGS:
            if clickData:
                return None
            else:
                return "Click a word on the plot to see its top 5 neighbors."
