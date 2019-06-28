from textwrap import dedent

import dash_core_components as dcc
import dash_html_components as html


def DemoDescriptionCard(markdown_text):
    """
    'width': '80%',
    'max-width': '1024px',
    'font-family': 'Roboto, sans-serif'
    :param markdown_text:
    :return: html.Div
    """
    return html.Div(
        className="row",
        style={
            "padding": "15px 30px 27px",
            "margin": "10px auto 45px",
            "width": "80%",
            "max-width": "1024px",
            "borderRadius": 5,
            "border": "thin lightgrey solid",
            "font-family": "Roboto, sans-serif",
        },
        children=dcc.Markdown(dedent(markdown_text)),
    )
