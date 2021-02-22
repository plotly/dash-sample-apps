import dash_html_components as html
import dash_core_components as dcc

MATERALIZE_CSS = (
    "https://cdnjs.cloudflare.com/ajax/libs/materialize/1.0.0/css/materialize.min.css"
)


def build_component(class_name, component=html.Div):
    def component_func(*args, className="", **kwargs):
        return component(*args, className=class_name + " " + className, **kwargs)

    return component_func


def Col(*args, width, **kwargs):
    class_name = f"col s{width}"
    return html.Div(*args, className=class_name, **kwargs)


def CustomSlider(id, min, max, label, **kwargs):
    mid = int((min + max) / 2)
    kwargs["value"] = kwargs.get("value", mid)
    return html.Div(
        [
            html.P(label),
            html.Br(),
            dcc.Slider(
                id=id,
                min=min,
                max=max,
                marks={i: str(i) for i in [min, mid, max]},
                tooltip={"always_visible": False},
                **kwargs,
            ),
        ]
    )


Row = build_component("row")
Card = build_component("card")
CardTitle = build_component("card-title")
CardContent = build_component("card-content")
CardAction = build_component("card-action")
