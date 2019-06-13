import dash_html_components as html
import dash_core_components as dcc


def Header():
    return html.Div([get_header(), html.Br([]), get_menu()])


def get_header():
    header = html.Div(
        [
            html.Div(
                [html.H5("Calibre Financial Index Fund Investor Shares")],
                className="seven columns main-title",
            ),
            html.Div(
                [
                    dcc.Link(
                        "Full View",
                        href="/dash-financial-report/full-view",
                        className="full-view-link",
                    )
                ],
                className="five columns",
            ),
        ],
        className="row",
        style={"padding-top": "25px"},
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
        className="row",
        style={"margin-bottom": "25px"},
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
