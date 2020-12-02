## Dash Holoviews Demo
This is a simple Demo app demonstrating how Dash, HoloViews, and Datashader can be used together to interactively explore a 10 million row dataset.

## Mapbox Token
This example requires a mapbox token that can be created from a free mapbox account at  https://www.mapbox.com/.

To run this app you must either:
 - Create an environment variable named `MAPBOX_TOKEN` that is set to your token string.
 - Create a file named `.mapbox` in the top-level project directory containing your token.
 
## Setup environment
Set up the app environment in a fresh virtual environment with

```
pip install -r requirements.txt
```

## Run the App
Then run the app in development mode with

```
python app.py
```

View the app by visiting the url displayed in the console in your web browser.
