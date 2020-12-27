from aktools.conf import *
from aktools.tools import *
from aktools.graphics import *


@app.callback(
    Output("live-graphs_host", "figure"),
    Input(component_id="hosts_dropdown", component_property="value"),
    Input("graph-update", "n_intervals"),
)
def update_output(value, interval):

    global connected_hosts, received_hosts, active_hosts

    init_datetime = int((time.time() - 5) * 1000)
    end_datetime = int(time.time() * 1000)
    past_files = sorted(
        [
            filename
            for filename in glob.glob("output/*.txt")
            if os.path.getmtime(filename) >= init_datetime / 1000 - log_ofo_time
        ],
        key=os.path.getmtime,
    )[::-1]

    for filename in past_files:
        connected_hosts += connected_to(filename, init_datetime, end_datetime, value)
        received_hosts += received_from(filename, init_datetime, end_datetime, value)
        active_hosts += generated_conn(filename, init_datetime, end_datetime)

    fig = make_subplots(
        rows=2,
        cols=2,
        specs=[[{"type": "domain"}, {"type": "domain"}], [{"colspan": 2}, None]],
        subplot_titles=(
            "Generated connections",
            "received connections",
            "total number of connections of all hosts",
        ),
    )

    fig.add_trace(
        go.Pie(
            labels=list(connected_hosts.keys()),
            values=list(connected_hosts.values()),
            textinfo="label+value",
            name="connected to",
            hole=0.65,
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Pie(
            labels=list(received_hosts.keys()),
            values=list(received_hosts.values()),
            textinfo="label+value",
            name="received from",
            hole=0.65,
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Bar(
            x=list(active_hosts.keys()),
            y=list(active_hosts.values()),
            name="All connections",
            marker=dict(color="orange", coloraxis="coloraxis"),
        ),
        row=2,
        col=1,
    )

    # fig.update_layout(
    # title_text="Host that connected and received connection from the selected host")

    return fig


if __name__ == "__main__":

    app.run_server(host="127.0.0.1", port=8089)
