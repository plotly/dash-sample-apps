import dash_html_components as html
import dash_core_components as dcc


def Header():
    return html.Div([get_logo(), get_header(), html.Br([]), get_menu()])


def get_logo():
    logo = html.Div(
        [
            html.Div(
                # [html.Img(src=app.get_asset_url("Logo.png"), height="40px")],
                className="ten columns padded"
            ),
            html.Div(
                [dcc.Link("Full View   ", href="/dash-financial-report/full-view")],
                className="two columns page-view no-print",
            ),
        ],
        className="row gs-header",
    )
    return logo


def get_header():
    header = html.Div(
        [
            html.Div(
                [html.H5("Calibre Financial Index Fund Investor Shares")],
                className="twelve columns padded",
            )
        ],
        className="row gs-header gs-text-header",
    )
    return header


def get_menu():
    menu = html.Div(
        [
            dcc.Link(
                "Overview   ",
                href="/dash-financial-report/overview",
                className="tab first",
            ),
            dcc.Link(
                "Price Performance   ",
                href="/dash-financial-report/price-performance",
                className="tab",
            ),
            dcc.Link(
                "Portfolio & Management   ",
                href="/dash-financial-report/portfolio-management",
                className="tab",
            ),
            dcc.Link(
                "Fees & Minimums   ",
                href="/dash-financial-report/fees",
                className="tab",
            ),
            dcc.Link(
                "Distributions   ",
                href="/dash-financial-report/distributions",
                className="tab",
            ),
            dcc.Link(
                "News & Reviews   ",
                href="/dash-financial-report/news-and-reviews",
                className="tab",
            ),
        ],
        className="row ",
    )
    return menu


def make_dash_table(df):
    """ Return a dash definition of an HTML table for a Pandas dataframe """
    table = []
    for index, row in df.iterrows():
        html_row = []
        for i in range(len(row)):
            html_row.append(html.Td([row[i]]))
        table.append(html.Tr(html_row))
    return table
