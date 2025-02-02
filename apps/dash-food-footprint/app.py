import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import plotly.graph_objs as go
import geojson
import os
from plotly.subplots import make_subplots
import dash_bootstrap_components as dbc

mapbox_access_token = "pk.eyJ1Ijoic3RlZmZlbmhpbGwiLCJhIjoiY2ttc3p6ODlrMG1ybzJwcG10d3hoaDZndCJ9.YE2gGNJiw6deBuFgHRHPjg"

dirname = os.path.dirname(__file__)
path = os.path.join(dirname, "data/")

with open(path + "diaphantinhenglish.geojson") as f:
    vn_map = geojson.load(f)
vn_tctk = pd.read_csv(path + 'tctk.csv', sep=',')
vn_tctk['Name_v'] = vn_tctk["Name"].map({
"An Giang":"An Giang",
"Bac Giang":"Bắc Giang",
"Bac Kan":"Bắc Cạn",
"Bac Lieu":"Bạc Liêu",
"Bac Ninh":"Bắc Ninh",
"Ben Tre":"Bến Tre",
"Ba Ria - Vung Tau":"Bà Rịa - Vũng Tàu",
"Binh Dinh":"Bình Định",
"Binh Duong":"Bình Dương",
"Binh Phuoc":"Bình Phước",
"Binh Thuan":"Bình Thuận",
"Can Tho":"Cần Thơ",
"Ca Mau":"Cà Mau",
"Cao Bang":"Cao Bằng",
"Dak Lak":"Đắc Lắc",
"Dak Nong":"Đắc Nông",
"Dong Nai":"Đồng Nai",
"Dong Thap":"Đồng Tháp",
"Da Nang":"Đà Nẵng",
"Dien Bien":"Điện Biên",
"Gia Lai":"Gia Lai",
"Hai Duong":"Hải Dương",
"Hai Phong":"Hải Phòng",
"Hau Giang":"Hậu Giang",
"Ha Giang":"Hà Giang",
"Ha Noi":"Hà Nội",
"Ha Nam":"Hà Nam",
"Ha Tinh":"Hà Tĩnh",
"Hoa Binh":"Hòa Bình",
"Hung Yen":"Hưng Yên",
"Khanh Hoa":"Khánh Hòa",
"Kien Giang":"Kiên Giang",
"Kon Tum":"Kon Tum",
"Lang Son":"Lạng Sơn",
"Lai Chau":"Lai Châu",
"Lam Dong":"Lâm Đồng",
"Lao Cai":"Lào Cai",
"Long An":"Long An",
"Nam Dinh":"Nam Định",
"Nghe An":"Nghệ An",
"Ninh Binh":"Ninh Bình",
"Ninh Thuan":"Ninh Thuận",
"Phu Tho":"Phú Thọ",
"Phu Yen":"Phú Yên",
"Quang Binh":"Quảng Bình",
"Quang Nam":"Quảng Nam",
"Quang Ngai":"Quảng Ngãi",
"Quang Ninh":"Quảng Ninh",
"Quang Tri":"Quảng Trị",
"Soc Trang":"Sóc Trăng",
"Son La":"Sơn La",
"Tay Ninh":"Tây Ninh",
"Thua Thien - Hue":"Thừa Thiên - Huế",
"Thai Binh":"Thái Bình",
"Thai Nguyen":"Thái Nguyên",
"Thanh Hoa":"Thanh Hóa",
"Tien Giang":"Tiền Giang",
"TP. Ho Chi Minh":"TP. Hồ Chí Minh",
"Tra Vinh":"Trà Vinh",
"Tuyen Quang":"Tuyên Quang",
"Vinh Long":"Vĩnh Long",
"Vinh Phuc":"Vĩnh Phúc",
"Yen Bai":"Yên Bái"
})
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
# tuổi
ages = pd.read_csv(path + 'age.csv', sep=',')
age_fig = go.Figure(
    data=[go.Bar(y=ages["percent"], x=ages['ages'],hovertemplate = '<extra>%{y}</extra>')],
    layout=go.Layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        yaxis=dict(title='%'),
        xaxis=dict(title='Tuổi')
    )
)

age_fig.update_yaxes(tickfont_size=7, showgrid=True, gridwidth=1, gridcolor='LightGrey', showline=True, linewidth=1,
                     linecolor='black')
age_fig.update_xaxes(showline=True, linewidth=1, linecolor='black', tickangle=-40)

# trình độ chuyên môn
training = pd.read_csv(path + 'training.csv', sep=',')
training['group_v'] = training['group'].map({
    "no professional knowledge": "Không có chuyên môn",
    "elementary level": "Sơ cấp",
    "intermediate level": "Trung cấp",
    "college degree": "Cao đẳng",
    "university degree and above": "Đại học và SĐH"
})

training_fig = go.Figure(
    data=[go.Bar(x=training["percent"], y=training['group_v'], orientation='h',hovertemplate = '<extra>%{x}</extra>')],
    layout=go.Layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(title='%'),
    )
)
training_fig.update_yaxes(tickfont_size=13, showline=True,
                          linewidth=1, linecolor='black')
training_fig.update_xaxes(showline=True, linewidth=1, linecolor='black', showgrid=True, gridwidth=1, gridcolor='LightGrey')

# hôn nhân
married = pd.read_csv(path + 'marriage.csv', sep=',')
married['group_v'] = married['marriage_status'].map({
    "not married": "Chưa kết hôn",
    "married": "Đã kết hôn",
    "widow/widower": "Góa chồng/vợ",
    "divorce": "Ly dị",
    "separated": "Ly thân"
})

married_fig = go.Figure(
    data=[go.Bar(x=married["percent"], y=married['group_v'], orientation='h',hovertemplate = '<extra>%{x}</extra>')],
    layout=go.Layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(title='%'),
    )
)

married_fig.update_yaxes(tickfont_size=13, showline=True,
                         linewidth=1, linecolor='black')
married_fig.update_xaxes(showline=True, linewidth=1, linecolor='black', showgrid=True, gridwidth=1, gridcolor='LightGrey')

# nguyên nhân
reason = pd.read_csv(path + 'reason.csv', sep=',')
reason['group_v'] = reason['reason'].map({
    "new job": "Tìm việc mới",
    "losse job": "Mất việc",
    "move new house with family": "Đi cùng gia đình",
    "marriage": "Kết hôn",
    "education": "Học tập",
    "others": "Khác"
})

reason_fig = go.Figure(
    data=[go.Bar(x=reason["percent"], y=reason['group_v'], orientation='h',hovertemplate = '<extra>%{x}</extra>')],
    layout=go.Layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(title='%'),
    )
)

reason_fig.update_yaxes(tickfont_size=13, showline=True, linewidth=1,
                        linecolor='black')
reason_fig.update_xaxes(showline=True, linewidth=1, linecolor='black', showgrid=True, gridwidth=1, gridcolor='LightGrey')

# Giới tính
gender = pd.read_csv(path + 'gender.csv', sep=',')
gender['group_v'] = gender['group'].map({
    "male": "Nam",
    "female": "Nữ",
})
# gender['Nam'] = gender['male']
# gender['Nữ'] = gender['female']

# gender_fig = px.bar(gender, y=["migrate"], x=["Nam", "Nữ"], labels={'value': "%", 'variable': ""})
# gender_fig.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
gender_fig = go.Figure(
    data=[go.Bar(x=gender["group_v"], y=gender['percent'], hovertemplate = '<extra>%{y}</extra>')],
    layout=go.Layout(
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        yaxis=dict(title='%'),
    )
)
gender_fig.update_xaxes(tickfont_size=13, showline=True, linewidth=1,
                        linecolor='black')
gender_fig.update_yaxes(showline=True, linewidth=1, linecolor='black', showgrid=True, gridwidth=1, gridcolor='LightGrey')

# Văn hóa
school = pd.read_csv(path + 'school.csv', sep=',')

school['Đang đi học'] = school['at shool']
school['Đã thôi học'] = school['stop learning']
school['Chưa đi học'] = school['never went to school']

school_fig = px.bar(school, x="Ages", y=['Đang đi học', 'Đã thôi học','Chưa đi học'], labels={'Ages': "Tuổi", 'value': "%", 'variable': ""})
school_fig.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
school_fig.update_xaxes(tickfont_size=13, showline=True, linewidth=1,
                        linecolor='black')
school_fig.update_yaxes(showline=True, linewidth=1, linecolor='black', showgrid=True, gridwidth=1, gridcolor='LightGrey')
school_fig.update_traces(hovertemplate='<extra>%{y}</extra>')

app = dash.Dash(__name__)

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
                    className="one-third column",
                ),
                dbc.Row([dbc.Col()]),
                html.Div(
                    [
                        html.Div(
                            [
                                html.H5(
                                    "Hiện trạng di cư giữa các tỉnh thành tại Việt Nam năm 2021",
                                    style={"font-weight": "bold"},
                                ),
                            ]
                        )
                    ],
                    className="three column",
                    id="title",
                ),
                html.Div(
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
                            html.Div(
                                [dcc.Graph(figure=age_fig)],
                                className="bare_container"
                            ),
                        ],
                        className="bare_container six columns",
                    ),
                        html.Div(
                            [
                                html.P("2. Hôn nhân",
                                       className="control_label",
                                       style={"text-align": "center", "font-weight": "bold"},
                                       ),
                                html.Div(
                                    [dcc.Graph(figure=married_fig)],
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
                            html.Div(
                                [dcc.Graph(figure=gender_fig)],
                                className="bare_container"
                            ),
                        ],
                        className="bare_container three columns",
                    ),
                        html.Div(
                            [
                                html.P("4. Trình độ chuyên môn",
                                       className="control_label",
                                       style={"text-align": "center", "font-weight": "bold"},
                                       ),
                                html.Div(
                                    [dcc.Graph(figure=training_fig)],
                                    className="bare_container"
                                ),
                            ],
                            className="bare_container five columns",
                        ),
                        html.Div(
                            [
                                html.P("5. Trình độ văn hóa",
                                       className="control_label",
                                       style={"text-align": "center", "font-weight": "bold"},
                                       ),
                                html.Div(
                                    [dcc.Graph(figure=school_fig)],
                                    className="bare_container"
                                ),
                            ],
                            className="bare_container four columns",
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
                        html.Div(
                            [dcc.Graph(figure=reason_fig)],
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
                         3. Dữ liệu địa không gian được sưu tầm tại trang web của tổ chức Sáng kiến phát triển mở Việt Nam (Open Development Vietnam – ODV): https://vietnam.opendevelopmentmekong.net/vi/about-us/.
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
    colors = ["blue", 'white', "red"]
    if candi == 'nhap_cu_2021':
        colors = ['white', "red"]
    if candi == 'xuat_cu_2021':
        colors = ['white', "blue"]
    if candi == 'di_cu_2021':
        midpoint = 0
    fig = px.choropleth_mapbox(
        vn_tctk,
        geojson=vn_map,
        color=candi,
        locations="Name",
        featureidkey="properties.Name",
        hover_name="Name_v",
        opacity=0.7,
        hover_data = {'Name':False},
        center={"lat": 16, "lon": 106},
        zoom=4.3,
        labels=labels_,
        color_continuous_scale=colors,
        color_continuous_midpoint=midpoint,
        # hoverinfo = "z"
    )
    fig.update_layout(
        margin={"r": 0, "t": 0, "l": 0, "b": 0}, mapbox_accesstoken=mapbox_access_token,
        legend=dict(title=food_options_[candi])
    )
    # fig.update_traces(hovertemplate="%{location}: %{z}", selector=dict(type='choroplethmapbox'))
    return fig

@app.callback(
    Output("indicator-graphic", "figure"),
    Input("yaxis-type", "value"),
)
def update_graph(yaxis_type):
    dat = vn_tctk.sort_values(by='di_cu_2021')
    if yaxis_type == "Kinh tế":
        fig = make_subplots(rows=1, cols=3, subplot_titles=("Di cư thuần (2021)", "Thu nhập (2021)", "Vốn đầu tư nước ngoài (2021)"))
        fig.add_trace(
            go.Bar(y=dat["Name_v"], x=dat["thu_nhap_2021"] / 1000, orientation='h', marker=dict(color='blue')),
            row=1, col=2
        )
        fig.add_trace(
            go.Bar(y=dat["Name_v"], x=dat["von_DTNN_2021"] / 1000, orientation='h', marker=dict(color='blue')),
            row=1, col=3
        )
        fig.update_layout(xaxis1=dict(title='Tỉ lệ di cư (\u2030)'),
                          xaxis2=dict(title='Thu nhập bình quân (triệu/tháng)'), xaxis3=dict(title='Vốn đầu tư nước ngoài (tỉ USD)'))

    else:
        fig = make_subplots(rows=1, cols=3, subplot_titles=("Di cư thuần", "Sinh viên đại học", "Sinh viên nghề"))

        fig.add_trace(
            go.Bar(y=dat["Name_v"], x=dat["SV_DH_2020"] / 1000, orientation='h', marker=dict(color='blue')),
            row=1, col=2
        )
        fig.add_trace(
            go.Bar(y=dat["Name_v"], x=dat["SV_nghe_2020"] / 1000, orientation='h', marker=dict(color='blue')),
            row=1, col=3
        )
        fig.update_layout(xaxis1=dict(title='Tỉ lệ di cư (\u2030)'), xaxis2=dict(title='Sinh viên đại học (nghìn)'),
                          xaxis3=dict(title='Sinh viên nghề (nghìn)'))
    colors = []

    for i in range(63):
        if dat.iloc[i, 3] < 0:
            colors.append('blue')
        else:
            colors.append('red')
    fig.add_trace(
        go.Bar(y=dat["Name_v"], x=dat['di_cu_2021'], orientation='h', marker=dict(color=colors)),
        row=1, col=1
    )
    fig.update_yaxes(tickfont_size=7, showgrid=True, gridwidth=0.5, gridcolor='LightGrey', showline=True, linewidth=1,
                     linecolor='black', mirror=True)
    fig.update_xaxes(showline=True, linewidth=1, linecolor='black', mirror=True, showgrid=True, gridwidth=1, gridcolor='LightGrey')
    fig.update_layout(height=750, yaxis2=dict(color='white'), yaxis3=dict(color='white'),
                      showlegend=False, paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    fig.update_traces(hovertemplate='<extra>(%{y}, %{x})</extra>')
    return fig


server = app.server

if __name__ == "__main__":
    app.run_server(debug=True)
