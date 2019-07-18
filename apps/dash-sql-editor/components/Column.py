import dash_html_components as html

def Column(children=None, width=1, **kwargs):
    number_mapping = {
        1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six',
        7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten', 11: 'eleven',
        12: 'twelve'
    }
    return html.Div(
        children,
        className="{} columns".format(number_mapping[width]),
        **kwargs
    )
