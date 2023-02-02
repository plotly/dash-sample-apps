import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import plotly

import pandas as pd
import plotly.graph_objs as go
import numpy as np

import geojson
import os
from plotly.subplots import make_subplots

mapbox_access_token = "pk.eyJ1Ijoic3RlZmZlbmhpbGwiLCJhIjoiY2ttc3p6ODlrMG1ybzJwcG10d3hoaDZndCJ9.YE2gGNJiw6deBuFgHRHPjg"

dirname = os.path.dirname(__file__)
path = os.path.join(dirname, "data/")

with open(path + "diaphantinhenglish.geojson") as f:
    vn_map = geojson.load(f)
vn_tctk = pd.read_csv(path + 'tctk.csv', sep=',')
# vn_tctk['xuat_cu_2021'] = -vn_tctk['xuat_cu_2021']

food_options_ = {
    "di_cu_2021": "Di cư thuần",
    "nhap_cu_2021": "Nhập cư",
    "xuat_cu_2021": "Xuất cư",

    # "thu_nhap_2021": "Thu nhập 2021",
    # "du_an_DTNN_2021": "DA ĐTNN 2021",
    # "von_DTNN_2021": "Vốn ĐTNN 2021",
    # "SV_DH_2020": "SV ĐH 2020",
    # "SV_nghe_2020": "SV nghề 2020",
}
labels_ = {
    "di_cu_2021": "Di cư thuần (\u2030)",
    "nhap_cu_2021": "Nhập cư (\u2030)",
    "xuat_cu_2021": "Xuất cư (\u2030)",

    # "thu_nhap_2021": "Thu nhập 2021",
    # "du_an_DTNN_2021": "DA ĐTNN 2021",
    # "von_DTNN_2021": "Vốn ĐTNN 2021",
    # "SV_DH_2020": "SV ĐH 2020",
    # "SV_nghe_2020": "SV nghề 2020",
}
food_options = [dict(label=food_options_[country], value=country) for country in food_options_.keys()]

radio_food_behaviour = dcc.RadioItems(
    id="nutrition_types",
    options=food_options,
    value="di_cu_2021",
    labelStyle={"display": "block", "text-align": "justify"},
)

ages = pd.read_csv(path + 'age.csv', sep=',')
age_fig = go.Figure(
    data=[go.Bar(y=ages["percent"], x=ages['ages'])],
    layout=go.Layout(
        title=go.layout.Title(text="1. Tuổi người di cư (năm 2019)"),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        yaxis=dict(title='%')
    )
)

age_fig.update_yaxes(tickfont_size=7, showgrid=True, gridwidth=1, gridcolor='LightPink', showline=True, linewidth=1,
                     linecolor='black', mirror=True)
age_fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True, tickangle=-40)

training = pd.read_csv(path + 'training.csv', sep=',')

training_fig = go.Figure(
    data=[go.Bar(x=training["percent"], y=training['group'], orientation='h')],
    layout=go.Layout(
        title=go.layout.Title(text="2. Trình độ chuyên môn"),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(title='%')
    )
)

training_fig.update_yaxes(tickfont_size=10, showgrid=True, gridwidth=1, gridcolor='LightPink', showline=True,
                          linewidth=1,
                          linecolor='black', mirror=True)
training_fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)

app = dash.Dash(__name__)


def get_fig_bar(col):
    dat = vn_tctk.sort_values(by='di_cu_2021')
    fig = go.Figure(
        data=[go.Bar(y=dat["Name"], x=dat[col], orientation='h')],
        layout=go.Layout(
            title=go.layout.Title(text=food_options_[col]),
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)'
        )
    )
    fig.update_layout(
        height=700,
        # title_text='GDP and Life Expectancy (Americas, 2007)'
    )
    fig.update_yaxes(tickfont_size=7, showgrid=True, gridwidth=1, gridcolor='LightPink', showline=True, linewidth=1,
                     linecolor='black', mirror=True)
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    return fig


fig_bar = get_fig_bar("di_cu_2021")
## FF ##

# Create app layout
app.layout = html.Div(
    [
        dcc.Store(id="aggregate_data"),
        # empty Div to trigger javascript file for graph resizing
        html.Div(id="output-clientside"),
        html.Div(
            [
                html.Div(
                    [
                        # html.Img(
                        #     src=app.get_asset_url("Nova_IMS.png"),
                        #     id="plotly-image",
                        #     style={
                        #         "height": "60px",
                        #         "width": "auto",
                        #         "margin-bottom": "25px",
                        #     },
                        # )
                    ],
                    className="one-third column",
                ),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H5(
                                    "Hiện trạng di cư giữa các tỉnh thành tại Việt Nam năm 2021",
                                    style={"font-weight": "bold"},
                                ),
                                html.H6(
                                    "Phân tích mối quan hệ giữa di cư và các yếu tố kinh tế - giáo dục",
                                    style={"margin-top": "0px"},
                                ),
                            ]
                        )
                    ],
                    className="three column",
                    id="title",
                ),
                html.Div(
                    # create empty div for align cente
                    className="one-third column",
                ),
            ],
            id="header",
            className="row flex-display",
            style={"margin-bottom": "25px"},
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.H6(
                            "Hiện trạng di cư",
                            style={
                                "margin-top": "0",
                                "font-weight": "bold",
                                "text-align": "center",
                            },
                        ),
                        html.P(
                            "Để tìm tìm nguyên nhân di cư, trước tiên ta cần xem xét phân tích hiện trạng di cư. Theo số liệu di cư năm 2021 của TCTK, tình trạng di cư tại các tỉnh thành phố tại Việt Nam như sau:",
                            className="control_label",
                            style={"text-align": "justify"},
                        ),
                        html.P(),
                        html.P(
                            "Chọn số liệu",
                            className="control_label",
                            style={"text-align": "center", "font-weight": "bold"},
                        ),
                        radio_food_behaviour,
                    ],
                    className="pretty_container four columns",
                    id="cross-filter-options",
                    style={"text-align": "justify"},
                ),
                html.Div(
                    [
                        html.Div(
                            [dcc.Graph(id="choropleth")],
                            className="pretty_container",
                        ),
                    ],
                    id="right-column",
                    className="eight columns",
                ),
            ],
            className="row flex-display",
        ),
        html.Div(
            [
                html.H5(
                    "Đặc điểm người di cư",
                    style={
                        "margin-top": "0",
                        "font-weight": "bold",
                        "text-align": "center",
                    },
                ),
                html.Div(
                    [html.Div(
                        [
                            html.P("1. Tuổi",
                                   className="control_label",
                                   style={"text-align": "center", "font-weight": "bold"},
                                   ),
                            html.Img(
                                src=app.get_asset_url("age.png"),
                                className="bare_container"
                            ),
                        ],
                        className="bare_container six columns",
                    ),
                        html.Div(
                            [
                                html.P("2. Trình độ",
                                       className="control_label",
                                       style={"text-align": "center", "font-weight": "bold"},
                                       ),
                                html.Img(
                                    src=app.get_asset_url("trinhdo.png"),
                                    className="bare_container"
                                ),
                            ],
                            className="bare_container six columns",
                        ),
                    ],
                    className="row no_border_container",
                ),
                html.Div(
                    [html.Div(
                        [
                            html.P("3. Giới tính",
                                   className="control_label",
                                   style={"text-align": "center", "font-weight": "bold"},
                                   ),
                            html.Img(
                                src=app.get_asset_url("gender.png"),
                                className="bare_container"
                            ),
                        ],
                        className="bare_container six columns",
                    ),
                        html.Div(
                            [
                                html.P("4. Hôn nhân",
                                       className="control_label",
                                       style={"text-align": "center", "font-weight": "bold"},
                                       ),
                                html.Img(
                                    src=app.get_asset_url("honnhan.png"),
                                    className="contain",
                                ),
                            ],
                            className="bare_container six columns",
                        ),
                    ],
                    className="row no_border_container",
                ),
            ],
            className="pretty_container",
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.H6(
                            "Tương quan giữa di cư và kinh tế - giáo dục",
                            style={
                                "margin-top": "0",
                                "font-weight": "bold",
                                "text-align": "center",
                            },
                        ),
                        dcc.RadioItems(
                            id="yaxis-type",
                            options=[
                                {"label": i, "value": i}
                                for i in ["Kinh tế", "Giáo dục"]
                            ],
                            value="Kinh tế",
                            labelStyle={"display": "inline-block"},
                            style={"padding-left": "43%"},
                        ),

                        html.Div(
                            [html.Div(
                                [dcc.Graph(id="indicator-graphic")],
                                className="no_border_container twelve columns",
                            ),
                            ],
                            className="row no_border_container",
                        ),
                    ],
                    className="no_border_container twelve columns",
                ),
            ],
            className="row pretty_container",
        ),
        html.Div(
            [
                html.Div(
                    className="no_border_container two columns",
                ),
                html.Div(
                    [
                        html.P("Kết quả khảo sát lý do di cư theo tổng điều tra dân số và nhà ở năm 2019",
                               className="control_label",
                               style={"text-align": "center", "font-weight": "bold"},
                               ),
                        html.Img(
                            src=app.get_asset_url("reason.png"),
                            className="bare_container"
                        ),
                    ],
                    className="no_border_container eight columns",
                ),
                html.Div(
                    className="no_border_container two columns",
                ),
            ],
            className="row pretty_container",
        ),
        html.Div(
            [
                html.H6(
                    "Người thực hiện",
                    style={
                        "margin-top": "0",
                        "font-weight": "bold",
                        "text-align": "center",
                    },
                ),
                html.P(
                    "Đào Thị Thu Hồng (21007975)  -  Lê Kim Dũng (21007978)",
                    style={"text-align": "center", "font-size": "10pt"},
                ),
            ],
            className="pretty_container",
        ),
        html.Div(
            [
                html.H6(
                    "Số liệu \n",
                    style={
                        "margin-top": "0",
                        "font-weight": "bold",
                        "text-align": "center",
                    },
                ),
                dcc.Markdown(
                    """\
                         1.	Số liệu tổng cục thống kê: https://www.gso.gov.vn/
                         2.	Số liệu tổng điều tra dân số và nhà ở 2019: https://www.gso.gov.vn/tong-dieu-tra-dan-so-va-nha-o/
                        """,
                    style={"font-size": "10pt"},
                ),
            ],
            className="pretty_container",
        ),
    ],
    id="mainContainer",
    style={"display": "flex", "flex-direction": "column"},
)


####
@app.callback(Output("choropleth", "figure"), [Input("nutrition_types", "value")])
def display_choropleth(candi):
    midpoint = None
    colors = color_continuous_scale=["blue",'white',"red" ]
    if candi == 'nhap_cu_2021':
        colors = color_continuous_scale = ['white',"red"]
    if candi == 'xuat_cu_2021':
        colors = color_continuous_scale = ['white',"blue"]
    if candi == 'di_cu_2021':
        midpoint = 0
    fig = px.choropleth_mapbox(
        vn_tctk,
        geojson=vn_map,
        color=candi,
        locations="Name",
        featureidkey="properties.Name",
        hover_name="Name",
        opacity=0.7,  # hover_data = [],
        center={"lat": 16, "lon": 106},
        zoom=4.3,
        labels=labels_,
        color_continuous_scale=colors,
        color_continuous_midpoint  = midpoint
    )
    fig.update_layout(
        margin={"r": 0, "t": 0, "l": 0, "b": 0}, mapbox_accesstoken=mapbox_access_token,
        legend=dict(title=food_options_[candi])
    )

    # fig.update_traces(legendgrouptitle_text=food_options_[candi], selector=dict(type='choreographically'))

    return fig


# @app.callback(
#     Output("indicator-graphic2", "figure"),
#     Output("indicator-graphic3", "figure"),
#     Input("yaxis-type", "value"),
# )
# def update_graph(yaxis_type):
#
#     if yaxis_type == "Kinh tế":
#         col2 = "thu_nhap_2021"
#         col3 = "von_DTNN_2021"
#     else:
#         col2 = "SV_DH_2020"
#         col3 = "SV_nghe_2020"
#
#
#     fig2 = get_fig_bar(col2)
#     fig3 = get_fig_bar(col3)
#     fig2.update_yaxes(showticklabels = False)
#     fig3.update_yaxes(showticklabels=False)
#     return fig2, fig3

@app.callback(
    Output("indicator-graphic", "figure"),
    Input("yaxis-type", "value"),
)
def update_graph(yaxis_type):
    dat = vn_tctk.sort_values(by='di_cu_2021')
    cols = plotly.colors.DEFAULT_PLOTLY_COLORS
    if yaxis_type == "Kinh tế":
        fig = make_subplots(rows=1, cols=3, subplot_titles=("Di cư thuần", "Thu nhập", "Vốn ĐTNN"))
        fig.add_trace(
            go.Bar(y=dat["Name"], x=dat["thu_nhap_2021"] / 1000, orientation='h',marker=dict(color=cols[0])),
            row=1, col=2
        )
        fig.add_trace(
            go.Bar(y=dat["Name"], x=dat["von_DTNN_2021"] / 1000, orientation='h',marker=dict(color=cols[0])),
            row=1, col=3
        )
        fig.update_layout(xaxis1=dict(title='Tỉ lệ di cư (\u2030)'),
                          xaxis2=dict(title='Thu nhập bình quân (triệu/tháng)'), xaxis3=dict(title='Vốn ĐTNN (tỉ USD)'))

    else:
        fig = make_subplots(rows=1, cols=3, subplot_titles=("Di cư thuần", "Sinh viên đại học", "Sinh viên nghề"))

        fig.add_trace(
            go.Bar(y=dat["Name"], x=dat["SV_DH_2020"] / 1000, orientation='h', marker=dict(color=cols[0])),
            row=1, col=2
        )
        fig.add_trace(
            go.Bar(y=dat["Name"], x=dat["SV_nghe_2020"] / 1000, orientation='h',marker=dict(color=cols[0])),
            row=1, col=3
        )
        fig.update_layout(xaxis1=dict(title='Tỉ lệ di cư (\u2030)'), xaxis2=dict(title='Sinh viên đại học (nghìn)'),
                          xaxis3=dict(title='Sinh viên nghề (nghìn)'))
    colors = []

    for i in range(63):
        if dat.iloc[i, 3] < 0:
            colors.append(cols[0])
        else:
            colors.append('red')
    fig.add_trace(
        go.Bar(y=dat["Name"], x=dat['di_cu_2021'], orientation='h', marker=dict(color=colors)),
        row=1, col=1
    )
    fig.update_yaxes(tickfont_size=7, showgrid=True, gridwidth=0.5, gridcolor='LightGrey', showline=True, linewidth=1,
                     linecolor='black', mirror=True)
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True)
    fig.update_layout(height=750, yaxis2=dict(showticklabels=False), yaxis3=dict(showticklabels=False),
                      showlegend=False, paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    return fig


server = app.server

if __name__ == "__main__":
    app.run_server(debug=True)
