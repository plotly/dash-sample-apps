# -*- coding: utf-8 -*-
import os
import dash

from demo import demo_layout, demo_callbacks
from local import local_layout, local_callbacks

app = dash.Dash(__name__)
app = dash.Dash(meta_tags=[{"name": "viewport", "content": "width=device-width"}])

server = app.server
app.layout = demo_layout
demo_callbacks(app)

# Load external CSS
external_css = [
    "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/skeleton/2.0.4/skeleton.min.css",
    "//fonts.googleapis.com/css?family=Raleway:400,300,600",
    "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
]

for css in external_css:
    app.css.append_css({"external_url": css})

# Running the server
if __name__ == "__main__":
    app.run_server(debug=True)
