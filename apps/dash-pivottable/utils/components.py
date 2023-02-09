from dash import html

def header(app, header_color, header, subheader=None, header_background_color="transparent"):
    left_headers = html.Div(
        [
            html.Div(header, className="header-title"),
            html.Div(subheader, className="subheader-title"),
        ],
        style={"color": header_color}
    )

    logo = html.Img(src=app.get_asset_url("images/plotly-logo-light-theme.png"))
    logo_link = html.A(logo, href="https://plotly.com/get-demo/", target="_blank")
    demo_link = html.A(
        "LEARN MORE",
        href="https://plotly.com/dash/",
        target="_blank",
        className="demo-button",
    )
    right_logos = html.Div([demo_link, logo_link], className="header-logos")

    return html.Div([left_headers, right_logos], className="header", style={"background-color": header_background_color})