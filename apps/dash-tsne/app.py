# -*- coding: utf-8 -*-
import os
import dash

from demo import create_layout, demo_callbacks

# for the Local version, import local_layout and local_callbacks
# from local import local_layout, local_callbacks

app = dash.Dash(
    __name__, meta_tags=[{"name": "viewport", "content": "width=device-width"}]
)
app.title = "t-SNE Explorer"

server = app.server
app.layout = create_layout(app)
demo_callbacks(app)

# Running server
if __name__ == "__main__":
    app.run_server(debug=True)
