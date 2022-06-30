import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash import html, dcc
from utils.figures import make_default_figure
from constants import class_labels, DEFAULT_STROKE_WIDTH, SEG_FEATURE_TYPES
from utils.helper_functions import class_to_color

# Modal
with open("explanations.md", "r") as f:
    howto_md = f.read()

modal_overlay = dbc.Modal(
    [
        dbc.ModalBody(html.Div([dcc.Markdown(howto_md)], id="howto-md")),
        dbc.ModalFooter(dbc.Button("Close", id="howto-close", className="howto-bn")),
    ],
    id="modal",
    size="lg",
)

button_howto = dbc.Button(
    "Learn more",
    id="howto-open",
    outline=True,
    color="info",
    # Turn off lowercase transformation for class .button in stylesheet
    style={"textTransform": "none"},
)

button_github = dbc.Button(
    "View Code on github",
    outline=True,
    color="primary",
    href="https://github.com/plotly/dash-sample-apps/tree/master/apps/dash-image-segmentation",
    id="gh-link",
    style={"text-transform": "none"},
)


# def header(
#     app, header_color, header, subheader=None, header_background_color="transparent"
# ):
#     left_headers = html.Div(
#         [
#             html.Div(header, className="header-title"),
#             html.Div(subheader, className="subheader-title"),
#         ],
#         style={"color": header_color},
#     )

#     logo = html.Img(src=app.get_asset_url("images/plotly-logo-light-theme.png"))
#     logo_link = html.A(logo, href="https://plotly.com/get-demo/", target="_blank")
#     demo_link = html.A(
#         "LEARN MORE",
#         href="https://plotly.com/dash/",
#         target="_blank",
#         className="demo-button",
#     )
#     right_logos = html.Div([demo_link, logo_link], className="header-logos")

#     return html.Div(
#         [left_headers, right_logos],
#         className="header",
#         style={"background-color": header_background_color},
#     )


# # DBC Header
# dbc_header = dbc.Navbar(
#     dbc.Container(
#         [
#             dbc.Row(
#                 [
#                     dbc.Col(
#                         html.Img(
#                             id="logo",
#                             src=("assets/images/plotly-logo-dark-theme.png"),
#                             height="30px",
#                         ),
#                         md="auto",
#                     ),
#                     dbc.Col(
#                         [
#                             html.Div(
#                                 [
#                                     html.H3("Interactive Machine Learning"),
#                                     html.P("Image segmentation"),
#                                 ],
#                                 id="app-title",
#                             )
#                         ],
#                         md=True,
#                         align="center",
#                     ),
#                 ],
#                 align="center",
#             ),
#             dbc.Row(
#                 [
#                     dbc.Col(
#                         [
#                             dbc.NavbarToggler(id="navbar-toggler"),
#                             dbc.Collapse(
#                                 dbc.Nav(
#                                     [
#                                         dbc.NavItem(button_howto),
#                                         dbc.NavItem(button_github),
#                                     ],
#                                     navbar=True,
#                                 ),
#                                 id="navbar-collapse",
#                                 navbar=True,
#                             ),
#                             modal_overlay,
#                         ],
#                         md=2,
#                     ),
#                 ],
#                 align="center",
#             ),
#         ],
#         fluid=True,
#     ),
#     dark=True,
#     color="dark",
#     sticky="top",
# )

header_items = dmc.Group(
    position="apart",
    children=[
        dmc.Image(
            src="assets/images/plotly-logo-light-theme.png", width=200, height=40
        ),
        dmc.Text(
            "Dash Image Segmentation",
            color="gray",
            size="xl",
            weight=600,
            transform="capitalize",
        ),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem(
                    "Behind the App", href="https://plotly.com/dash/design-kit/"
                ),
                dbc.DropdownMenuItem(
                    "View Code on Github",
                    href="https://github.com/plotly/dash-sample-apps/blob/main/apps/dash-image-segmentation/app.py",
                ),
            ],
            nav=True,
            in_navbar=True,
            label="Learn More",
        ),
    ],
)


# Description
description = dbc.Col(
    [
        dbc.Card(
            id="description-card",
            children=[
                dbc.CardHeader("Explanation"),
                dbc.CardBody(
                    [
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        html.Img(
                                            src="assets/images/segmentation_img_example_marks.jpg",
                                            width="200px",
                                        )
                                    ],
                                    md="auto",
                                ),
                                dbc.Col(
                                    html.P(
                                        "This is an example of interactive machine learning for image classification. "
                                        "To train the classifier, draw some marks on the picture using different colors for "
                                        'different parts, like in the example image. Then enable "Show segmentation" to see the '
                                        "classes a Random Forest Classifier gave to regions of the image, based on the marks you "
                                        "used as a guide. You may add more marks to clarify parts of the image where the "
                                        "classifier was not successful and the classification will update."
                                    ),
                                    md=True,
                                ),
                            ]
                        ),
                    ]
                ),
            ],
        )
    ],
    md=12,
)

# Image Segmentation
segmentation = [
    dbc.Card(
        id="segmentation-card",
        children=[
            dbc.CardHeader("Viewer"),
            dbc.CardBody(
                [
                    # Wrap dcc.Loading in a div to force transparency when loading
                    html.Div(
                        id="transparent-loader-wrapper",
                        children=[
                            dcc.Loading(
                                id="segmentations-loading",
                                type="cube",
                                children=[
                                    # Graph
                                    dcc.Graph(
                                        id="graph",
                                        figure=make_default_figure(),
                                        config={
                                            "modeBarButtonsToAdd": [
                                                "drawrect",
                                                "drawopenpath",
                                                "eraseshape",
                                            ]
                                        },
                                    ),
                                ],
                            )
                        ],
                    ),
                ]
            ),
            dbc.CardFooter(
                [
                    # Download links
                    html.A(
                        id="download",
                        download="classifier.json",
                    ),
                    html.Div(
                        children=[
                            dbc.ButtonGroup(
                                [
                                    dbc.Button(
                                        "Download classified image",
                                        id="download-image-button",
                                        outline=True,
                                    ),
                                    dbc.Button(
                                        "Download classifier",
                                        id="download-button",
                                        outline=True,
                                    ),
                                ],
                                size="lg",
                                style={"width": "100%"},
                            ),
                        ],
                    ),
                    html.A(
                        id="download-image",
                        download="classified-image.png",
                    ),
                ]
            ),
        ],
    )
]

# sidebar
sidebar = [
    dbc.Card(
        id="sidebar-card",
        children=[
            dbc.CardHeader("Tools"),
            dbc.CardBody(
                [
                    html.H6("Label class", className="card-title"),
                    # Label class chosen with buttons
                    html.Div(
                        id="label-class-buttons",
                        children=[
                            dbc.Button(
                                "%2d" % (n,),
                                id={"type": "label-class-button", "index": n},
                                style={"background-color": class_to_color(c)},
                            )
                            for n, c in enumerate(class_labels)
                        ],
                    ),
                    html.Hr(),
                    dbc.Form(
                        [
                            dbc.Row(
                                [
                                    dbc.Label(
                                        "Width of annotation paintbrush",
                                        html_for="stroke-width",
                                    ),
                                    # Slider for specifying stroke width
                                    dcc.Slider(
                                        id="stroke-width",
                                        min=0,
                                        max=6,
                                        step=1,
                                        value=DEFAULT_STROKE_WIDTH,
                                    ),
                                ]
                            ),
                            dbc.Row(
                                [
                                    html.H6(
                                        id="stroke-width-display",
                                        className="card-title",
                                    ),
                                    dbc.Label(
                                        "Blurring parameter",
                                        html_for="sigma-range-slider",
                                    ),
                                    dcc.RangeSlider(
                                        id="sigma-range-slider",
                                        min=0.01,
                                        max=20,
                                        step=0.01,
                                        value=[0.5, 16],
                                    ),
                                ]
                            ),
                            dbc.Row(
                                [
                                    dbc.Label(
                                        "Select features",
                                        html_for="segmentation-features",
                                    ),
                                    dcc.Checklist(
                                        id="segmentation-features",
                                        options=[
                                            {"label": l.capitalize(), "value": l}
                                            for l in SEG_FEATURE_TYPES
                                        ],
                                        value=["intensity", "edges"],
                                    ),
                                ]
                            ),
                            # Indicate showing most recently computed segmentation
                            dcc.Checklist(
                                id="show-segmentation",
                                options=[
                                    {
                                        "label": "Show segmentation",
                                        "value": "Show segmentation",
                                    }
                                ],
                                value=[],
                            ),
                        ]
                    ),
                ]
            ),
        ],
    ),
]

meta = [
    html.Div(
        id="no-display",
        children=[
            # Store for user created masks
            # data is a list of dicts describing shapes
            dcc.Store(id="masks", data={"shapes": []}),
            dcc.Store(id="classifier-store", data={}),
            dcc.Store(id="classified-image-store", data=""),
            dcc.Store(id="features_hash", data=""),
        ],
    ),
    html.Div(id="download-dummy"),
    html.Div(id="download-image-dummy"),
]
