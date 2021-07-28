# Dash Sample Apps 

[![CircleCI](https://circleci.com/gh/plotly/dash-sample-apps.svg?style=svg)](https://circleci.com/gh/plotly/dash-sample-apps)

This repository hosts the code for over 100 open-source Dash apps written in Python or R. They can serve as a starting point for your own Dash app, as a learning tool to better understand how Dash works, as a reusable templates, and much more.

Most apps in this repository are hosted on [Dash Gallery](https://dash-gallery.plotly.host/), which is our internal server running on [Dash Enterprise Kubernetes](https://plotly.com/dash/kubernetes/). Note that you can find both open-sourced apps and demos for our [licensed products](https://plotly.com/dash/), including [Design Kit](https://plotly.com/dash/design-kit/) and [Snapshot Engine](https://plotly.com/dash/snapshot-engine/). If you are interested in learning more, don't hesitate to reach out to [get a demo](https://plotly.com/get-demo/). If you want to only see the open-sourced apps, select the ["Open Source" tag](https://dash-gallery.plotly.host/Portal/?search=[Open%20Source]).

## Downloading and running a single app

Visit the [releases page](https://github.com/plotly/dash-sample-apps/releases) and download and `unzip` the app you want. Then `cd` into the app directory and install its dependencies in a virtual environment in the following way:

```bash
python -m venv venv
source venv/bin/activate  # Windows: \venv\scripts\activate
pip install -r requirements.txt
```

then you can run the app:
```bash
python app.py
```

## Cloning this whole repository

To clone this repository, run:
```
git clone https://github.com/plotly/dash-sample-apps
```

Note this might take a long time since it copies over 100 apps available in the repo. If you just want to try one app, refer to the section above.

## Contributing

To contribute to this repository, please see the [contributing guide](CONTRIBUTING.md).
