from dash import Dash, html, dcc, Dash, Input, Output, State, callback, callback_context, dash_table
import pandas as pd

from constants import df, params, state_dict, max_length, suffix_button_id, suffix_sparkline_graph, suffix_count, suffix_ooc_n, suffix_ooc_g, suffix_indicator
from utils.helper_functions import update_sparkline, init_value_setter_store, populate_ooc, update_count
from utils.graphs import generate_graph, build_chart_panel
import utils.components as comp

app = Dash(
    __name__,
    title="Manufacturing SPC Dashboard",
    update_title=None
)
server = app.server

app.layout = html.Div(
    id="big-app-container",
    children=[
        comp.build_banner(),
        dcc.Interval(
            id="interval-component",
            interval=2 * 1000,  # in milliseconds
            n_intervals=50,  # start at batch 50
            disabled=True,
        ),
        html.Div(
            id="app-container",
            children=[
                comp.build_tabs(),
                html.Div(id="app-content"), # Main app
            ],
        ),
        dcc.Store(id="value-setter-store", data=init_value_setter_store(df)),
        dcc.Store(id="n-interval-stage", data=50),
        comp.generate_modal(),
    ],
)


@callback(
    Output("app-content", "children"), 
    Output("interval-component", "n_intervals"),
    Input("app-tabs", "value"),
    State("n-interval-stage", "data"),
)
def render_tab_content(tab_switch, stopped_interval):
    if tab_switch == "tab1":
        return comp.build_tab_1(), stopped_interval
    return html.Div(
            id="status-container",
            children=[
                comp.build_quick_stats_panel(),
                html.Div(
                    id="graphs-container",
                    children=[comp.build_top_panel(stopped_interval), build_chart_panel()],
                ),
            ],
        ), stopped_interval


@callback(
    Output("n-interval-stage", "data"),
    Input("app-tabs", "value"),
    State("interval-component", "n_intervals"),
    State("interval-component", "disabled"),
    State("n-interval-stage", "data"),
)
def update_interval_state(tab_switch, cur_interval, disabled, cur_stage):
    if disabled:
        return cur_interval

    if tab_switch == "tab1":
        return cur_interval
    return cur_stage


@callback(
    Output("interval-component", "disabled"), 
    Output("stop-button", "buttonText"),
    Input("stop-button", "n_clicks"),
    State("interval-component", "disabled"),
)
def stop_production(n_clicks, current):
    """
        Callbacks for stopping interval update
    """
    if n_clicks == 0:
        return True, "start"
    return not current, "stop" if current else "start"


@callback(
    Output("markdown", "style"),
    Input("learn-more-button", "n_clicks"), 
    Input("markdown_close", "n_clicks"),
)
def update_click_output(button_click, close_click):
    """
        Callbacks for modal popup
    """
    ctx = callback_context

    if ctx.triggered:
        prop_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if prop_id == "learn-more-button":
            return {"display": "block"}

    return {"display": "none"}


@callback(
    Output("progress-gauge", "value"),
    Input("interval-component", "n_intervals"),
)
def update_gauge(interval):
    """
        Update progress gauge
    """
    if interval < max_length:
        total_count = interval
    else:
        total_count = max_length

    return int(total_count)


@callback(
    Output("ud_usl_input", "value"),
    Output("ud_lsl_input", "value"),
    Output("ud_ucl_input", "value"),
    Output("ud_lcl_input", "value"),
    Input("metric-select-dropdown", "value"),
    State("value-setter-store", "data"),
)
def build_value_setter_panel(dd_select, state_value):
    """
        Update values based on store data and dropdown selection
    """
    return state_value[dd_select]["usl"], state_value[dd_select]["lsl"], state_value[dd_select]["ucl"], state_value[dd_select]["lcl"],


@callback(
    Output("value-setter-store", "data"),
    Input("value-setter-set-btn", "n_clicks"),
    State("metric-select-dropdown", "value"),
    State("value-setter-store", "data"),
    State("ud_usl_input", "value"),
    State("ud_lsl_input", "value"),
    State("ud_ucl_input", "value"),
    State("ud_lcl_input", "value"),
)
def set_value_setter_store(set_btn, param, data, usl, lsl, ucl, lcl):
    """
        Update stored data via click and recalculate ooc in case of param updates
    """
    if set_btn is not None:
        data[param]["usl"] = usl
        data[param]["lsl"] = lsl
        data[param]["ucl"] = ucl
        data[param]["lcl"] = lcl

        data[param]["ooc"] = populate_ooc(df[param], ucl, lcl)
    return data


@callback(
    Output("value-setter-view-output", "children"),
    Input("value-setter-view-btn", "n_clicks"),
    Input("metric-select-dropdown", "value"),
    Input("value-setter-store", "data"),
)
def show_current_specs(n_clicks, dd_select, store_data):
    if n_clicks > 0:
        curr_col_data = store_data[dd_select]
        new_df_dict = {
            "Specs": [
                "Upper Specification Limit",
                "Lower Specification Limit",
                "Upper Control Limit",
                "Lower Control Limit",
            ],
            "Current Setup": [
                curr_col_data["usl"],
                curr_col_data["lsl"],
                curr_col_data["ucl"],
                curr_col_data["lcl"],
            ],
        }
        new_df = pd.DataFrame.from_dict(new_df_dict)
        return dash_table.DataTable(
            style_header={"fontWeight": "bold", "color": "inherit"},
            style_as_list_view=True,
            fill_width=True,
            style_cell_conditional=[
                {"if": {"column_id": "Specs"}, "textAlign": "left"}
            ],
            style_cell={
                "backgroundColor": "#1e2130",
                "fontFamily": "Open Sans",
                "padding": "0 2rem",
                "color": "darkgray",
                "border": "none",
            },
            css=[
                {"selector": "tr:hover td", "rule": "color: #91dfd2 !important;"},
                {"selector": "td", "rule": "border: none !important;"},
                {
                    "selector": ".dash-cell.focused",
                    "rule": "background-color: #1e2130 !important;",
                },
                {"selector": "table", "rule": "--accent: #1e2130;"},
                {"selector": "tr", "rule": "background-color: transparent"},
            ],
            data=new_df.to_dict("rows"),
            columns=[{"id": c, "name": c} for c in ["Specs", "Current Setup"]],
        )


# decorator for list of output
def create_callback(param):
    def callback(interval, stored_data):
        count, ooc_n, ooc_g_value, indicator = update_count(max_length, interval, param, stored_data)
        spark_line_data = update_sparkline(state_dict, max_length, interval, param)
        return count, spark_line_data, ooc_n, ooc_g_value, indicator

    return callback


for param in params[1:]:
    update_param_row_function = create_callback(param)
    callback(
        output=[
            Output(param + suffix_count, "children"),
            Output(param + suffix_sparkline_graph, "extendData"),
            Output(param + suffix_ooc_n, "children"),
            Output(param + suffix_ooc_g, "value"),
            Output(param + suffix_indicator, "color"),
        ],
        inputs=[Input("interval-component", "n_intervals")],
        state=[State("value-setter-store", "data")],
    )(update_param_row_function)


@callback(
    Output("control-chart-live", "figure"),
    Input("interval-component", "n_intervals"),
    Input(params[1] + suffix_button_id, "n_clicks"),
    Input(params[2] + suffix_button_id, "n_clicks"),
    Input(params[3] + suffix_button_id, "n_clicks"),
    Input(params[4] + suffix_button_id, "n_clicks"),
    Input(params[5] + suffix_button_id, "n_clicks"),
    Input(params[6] + suffix_button_id, "n_clicks"),
    Input(params[7] + suffix_button_id, "n_clicks"),
    State("value-setter-store", "data"), 
    State("control-chart-live", "figure"),
)
def update_control_chart(interval, n1, n2, n3, n4, n5, n6, n7, data, cur_fig):
    """
        button to choose/update figure based on click 
    """
    ctx = callback_context # Find which one has been triggered

    if not ctx.triggered:
        return generate_graph(interval, data, params[1])

    if ctx.triggered:
        # Get most recently triggered id and prop_type
        splitted = ctx.triggered[0]["prop_id"].split(".")
        prop_id = splitted[0]
        prop_type = splitted[1]

        if prop_type == "n_clicks":
            curr_id = cur_fig["data"][0]["name"]
            prop_id = prop_id[:-7]
            if curr_id == prop_id:
                return generate_graph(interval, data, curr_id)
            else:
                return generate_graph(interval, data, prop_id)

        if prop_type == "n_intervals" and cur_fig is not None:
            curr_id = cur_fig["data"][0]["name"]
            return generate_graph(interval, data, curr_id)


@callback(
    Output("piechart", "figure"),
    Input("interval-component", "n_intervals"),
    State("value-setter-store", "data"),
)
def update_piechart(interval, stored_data):
    if interval == 0:
        return {
            "data": [],
            "layout": {
                "font": {"color": "white"},
                "paper_bgcolor": "rgba(0,0,0,0)",
                "plot_bgcolor": "rgba(0,0,0,0)",
            },
        }

    if interval >= max_length:
        total_count = max_length - 1
    else:
        total_count = interval - 1

    values = []
    colors = []
    for param in params[1:]:
        ooc_param = (stored_data[param]["ooc"][total_count] * 100) + 1
        values.append(ooc_param)
        if ooc_param > 6:
            colors.append("#f45060")
        else:
            colors.append("#91dfd2")

    new_figure = {
        "data": [
            {
                "labels": params[1:],
                "values": values,
                "type": "pie",
                "marker": {"colors": colors, "line": dict(color="white", width=2)},
                "hoverinfo": "label",
                "textinfo": "label",
            }
        ],
        "layout": {
            "margin": dict(t=20, b=50),
            "uirevision": True,
            "font": {"color": "white"},
            "showlegend": False,
            "paper_bgcolor": "rgba(0,0,0,0)",
            "plot_bgcolor": "rgba(0,0,0,0)",
            "autosize": True,
        },
    }
    return new_figure


if __name__ == "__main__":
    app.run_server(debug=True)