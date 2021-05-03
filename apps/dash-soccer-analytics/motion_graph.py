import pandas as pd
import plotly.express as px
import os
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import dash
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from time import perf_counter
from datetime import datetime

# Disable non-applicable warnings
pd.options.mode.chained_assignment = None  # default='warn'

# Plots moving football data in plotly using built-in animation functionality
# This graph is based on Metrica Elite data and would need modified (minimally) for any other source data
# Execute this script in order to create the .json file that will be used by the main app. This file
# does NOT get called by the main app. Its only purpose is to pre-process tracking data so that it can be
# displayed in a timely fashion. Otherwise it just takes way too long to render more than
# 5 minuts of actino on demand
def game_simulator(file, half, start, stop):
    # start perf monitoring so we see how long it takes
    t1_start = perf_counter()
    print("Process started at: ", datetime.now())

    file_plus_path = "data/" + file
    df = pd.read_csv(file_plus_path, error_bad_lines=False)

    df["time"] = df["time"].astype(float)
    df["half"] = df["half"].astype(int)

    # Filter dataframe to include only activity within the time window provided
    df = df[(df["half"] == half) & (df["time"] >= start) & (df["time"] <= stop)]

    # Set the marker size for scatterplot points
    df["size"] = 10

    # Make the ball marker size smaller than the other markers
    df.loc[df["jersey_number"] == 0, "size"] = 3

    # Hide substitutes and ball out of bounds by making their size 0
    df.loc[df["x"] == None, "size"] = 0
    df.loc[df["y"] == None, "size"] = 0
    df.loc[df["x"] == None, "jersey_number"] = ""
    df.loc[df["y"] == None, "jersey_number"] = ""

    # Limit the dataframe to include only columns that we're going to use
    df = df[["half", "time", "x", "y", "team", "size", "jersey_number"]]

    # Replace nan with None id df
    df = df.where(pd.notnull(df), None)

    # Rename teamId column so it looks nicer when displayed on Legend
    df = df.rename(columns={"team": "Team"})

    # Limit time column to two decimals
    df["time"] = df["time"].astype(float)
    df["time"] = df.time.round(decimals=2)

    # else:
    # df = df.sort_values(by=['half', 'time'])

    colour0 = "#009BFF"
    colour1 = "grey"
    colour_ball = "red"

    # Make sure jersey number is an int not a decimal
    df["jersey_number"] = (
        df["jersey_number"].astype("str").replace("\.0", "", regex=True)
    )

    # null out jersey number field for ball or it shows a zero on the marker, which looks stupid
    df["jersey_number"].replace("0", "", inplace=True)

    color_discrete_map = {"Home": colour0, "Away": colour1, "Ball": colour_ball}

    # Plotly Express version
    fig = px.scatter(
        df,
        x="x",
        y="y",
        color="Team",
        hover_name="jersey_number",
        animation_frame="time",
        animation_group="jersey_number",
        range_x=[-0.05, 1.05],
        range_y=[-0.05, 1.05],
        size="size",
        size_max=10,
        opacity=0.8,
        color_discrete_map=color_discrete_map,
        text="jersey_number",
        hover_data={
            "x": False,
            "y": False,
            "time": False,
            "size": False,
            "Team": False,
            "jersey_number": False,
        },
    )

    # Add corner flags to prevent zoom and pitch distortion
    fig.add_scatter(
        x=[0, 0, 1, 1],
        y=[0, 1, 0, 1],
        mode="markers",
        marker=dict(size=1, color="grey"),
        name="Flags",
    )

    # Make jersey number really small inside markers
    fig.update_traces(
        textfont_size=7, textfont_color="white", hovertemplate=None, hoverinfo="none"
    )
    fig.update_yaxes(autorange="reversed")

    fig.update_layout(
        xaxis=dict(range=[-0.05, 1.05]),
        yaxis=dict(range=[-0.05, 1.05]),
        coloraxis_showscale=False,
    )

    # Remove side color scale and hide zero and gridlines
    fig.update_layout(
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False),
    )

    # Disable axis ticks and labels
    fig.update_xaxes(showticklabels=False, title_text="")
    fig.update_yaxes(showticklabels=False, title_text="")

    image_file = "assets/Pitch.png"
    image_path = os.path.join(os.getcwd(), image_file)

    from PIL import Image

    img = Image.open(image_path)

    fig.add_layout_image(
        dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=0,
            sizex=1,
            sizey=1,
            sizing="stretch",
            opacity=0.7,
            layer="below",
        )
    )

    ############  Playback Speed Setting  ##################
    fig.layout.updatemenus[0].buttons[0].args[1]["frame"]["duration"] = 200
    ###############################################

    fig["layout"]["sliders"][0]["pad"]["t"] = 0
    fig["layout"]["updatemenus"][0]["pad"]["t"] = 0

    pio.templates["custom_dark"] = go.layout.Template()
    pio.templates["custom_dark"]["layout"]["paper_bgcolor"] = "#282828"
    pio.templates["custom_dark"]["layout"]["plot_bgcolor"] = "#282828"

    fig.update_layout(
        template="custom_dark",
        xaxis=dict(showgrid=False, showticklabels=False),
        # plot_bgcolor='#282828',
        # paper_bgcolor='#282828'
    )

    # slider format and adjustments for aesthetic purposes
    fig["layout"]["sliders"][0]["pad"] = dict(r=0, t=0.0,)
    fig["layout"]["sliders"][0]["minorticklen"] = 2
    fig["layout"]["sliders"][0]["ticklen"] = 5
    fig["layout"]["sliders"][0]["tickcolor"] = "grey"
    fig["layout"]["sliders"][0]["font"]["color"] = "grey"
    fig["layout"]["sliders"][0]["bgcolor"] = "grey"
    fig["layout"]["sliders"][0]["bordercolor"] = "grey"
    fig["layout"]["template"]["data"]["scatter"][0]["marker"]["line"][
        "color"
    ] = "lightgrey"
    fig["layout"]["template"]["data"]["scatter"][0]["marker"]["opacity"] = 0.9

    fig.update_layout(margin=dict(l=20, r=20, b=20, t=20))

    fig.update_layout(
        legend_orientation="v", transition={"duration": 0, "ordering": "traces first"}
    )

    # Make sure pitch background image shape doesn't get distorted
    fig.update_yaxes(scaleanchor="x", scaleratio=0.65)

    t1_stop = perf_counter()
    print("Process took " + str((t1_stop - t1_start) / 60) + " minutes of time")

    fig.update_layout(legend=dict(yanchor="top", y=0.95, xanchor="left", x=-0.08))
    fig.update_layout(
        updatemenus=[
            dict(
                type="buttons",
                showactive=False,
                y=-0.14,
                x=-0.08,
                xanchor="left",
                yanchor="bottom",
            )
        ]
    )
    fig.update_layout(autosize=True, hovermode="closest")
    # fig.update_layout(showlegend=False)
    fig.update_layout(legend=dict(font=dict(family="Arial", size=10, color="grey")))

    # Hide corner flag trace in the legend
    for trace in fig["data"]:
        if trace["name"] == "Flags":
            trace["showlegend"] = False

    export = input("Do you wish to export the graph to json (y/n)?:")
    if export == "y":
        export_file_name = input(
            "Please enter a name for the json file to be exported (ending with .json): "
        )
        export_file_name = "data/" + export_file_name
        with open(export_file_name, "w") as f:
            # json_data = fig.to_json()
            pio.write_json(fig, f)
            f.close()

    return fig


app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.CYBORG],
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)
filename = input("Filename: ")
half = input("Half: ")
half = int(half)
start_time = input("Start Time: ")
start_time = int(start_time)
end_time = input("End Time: ")
end_time = int(end_time)

app.layout = html.Div(
    [
        dcc.Graph(
            figure=game_simulator(filename, int(half), int(start_time), int(end_time))
        )
    ]
)

start_server = input("Do you want to display the graph (y/n)?: ")
if start_server == "y":
    app.run_server(debug=True, use_reloader=False)
