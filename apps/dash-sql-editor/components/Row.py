import dash_html_components as html

def Row(children=None, **kwargs):
    return html.Div(
        children,
        className="row mt-2",
        **kwargs
    )
