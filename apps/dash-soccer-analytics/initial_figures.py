import plotly.express as px
import os
import plotly.graph_objects as go

# Create initial placeholder figure for game simulator
def initial_figure_simulator():
    # fig = px.scatter(x=[0, 0, 105, 105], y=[69, -2, 69, -2])
    fig = px.scatter(x=[0, 0, 1, 1], y=[0, 1, 0, 1])
    fig.update_layout(xaxis=dict(range=[0, 1]))
    fig.update_layout(yaxis=dict(range=[0, 1]))
    fig.update_traces(marker=dict(color="white", size=6))

    # Remove side color scale and hide zero and gridlines
    fig.update_layout(
        coloraxis_showscale=False,
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False),
        autosize=True
        # width=900,
        # height=600
    )

    # Disable axis ticks and labels
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
    fig.update_xaxes(title_text="")
    fig.update_yaxes(title_text="")
    fig.update_layout(margin=dict(l=80, r=80, b=10, t=20))
    fig.update_layout(modebar=dict(bgcolor="rgba(0, 0, 0, 0)"))
    image_file = "assets/Pitch.png"
    image_path = os.path.join(os.getcwd(), image_file)

    # Import and use pre-fabricated football pitch image
    from PIL import Image

    img = Image.open(image_path)

    fig.add_layout_image(
        dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=1,
            sizex=1,
            sizey=1,
            sizing="stretch",
            opacity=0.7,
            layer="below",
        )
    )

    fig.update_yaxes(scaleanchor="x", scaleratio=0.70)

    fig.update_layout(autosize=True)

    fig.update_layout(
        template="plotly_dark",
        xaxis=dict(showgrid=False, showticklabels=False),
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
    )
    return fig


# Create initial placeholder figure for event plot
def initial_figure_events():
    # fig = px.scatter(x=[0], y=[0])
    fig = px.scatter(x=[0, 0, 1, 1], y=[0, 1, 0, 1])
    fig.update_traces(marker=dict(color="white", size=6))
    fig.update_layout(yaxis=dict(range=[0, 1]))
    fig.update_layout(xaxis=dict(range=[0, 1]))
    fig.update_layout(margin=dict(l=10, r=100, b=10, t=45))

    fig.update_layout(modebar=dict(bgcolor="rgba(0, 0, 0, 0)"))

    # Remove side color scale and hide zero and gridlines
    fig.update_layout(
        coloraxis_showscale=False,
        xaxis=dict(showgrid=False, zeroline=False),
        yaxis=dict(showgrid=False, zeroline=False),
    )

    fig.update_layout(
        template="plotly_dark",
        xaxis=dict(showgrid=False, showticklabels=False),
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
    )

    # Disable axis ticks and labels
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)
    fig.update_xaxes(title_text="")
    fig.update_yaxes(title_text="")
    fig.update_xaxes(fixedrange=True)
    fig.update_yaxes(scaleanchor="x", scaleratio=0.70)
    fig.update_layout(margin=dict(l=10, r=30, b=30, t=30), autosize=True)
    image_file = "assets/Pitch.png"

    from PIL import Image

    img = Image.open(image_file)
    fig.add_layout_image(
        dict(
            source=img,
            xref="x",
            yref="y",
            x=0,
            y=1,
            sizex=1,
            sizey=1,
            sizing="stretch",
            opacity=0.8,
            layer="below",
        )
    )

    return fig


# Create initial placeholder figure for radar plot
def initial_figure_radar():
    fig = go.Figure()
    r_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    fig.add_trace(
        go.Scatterpolar(
            r=r_values,
            theta=[
                "Set<br>Plays",
                "Passes",
                "Balls<br>Lost",
                "Recoveries",
                "Challenges",
                "Shots",
                "Interceptions",
                "Crosses",
                "Long<br>Balls",
                "Free<br>Kicks",
            ],
            fill="toself",
            opacity=0.50,
        )
    )

    fig.update_layout(
        template="plotly_dark",
        plot_bgcolor="rgba(0, 0, 0, 0)",
        paper_bgcolor="rgba(0, 0, 0, 0)",
    )

    fig.update_layout(modebar=dict(bgcolor="rgba(0, 0, 0, 0)"))

    fig.update_layout(
        polar=dict(
            bgcolor="#282828",
            radialaxis=dict(visible=True, range=[0, 1], showticklabels=False,),
        ),
        showlegend=False,
    )
    fig.update_layout(autosize=True)
    fig.update_layout(margin=dict(l=60, r=60, b=30, t=45))
    return fig
