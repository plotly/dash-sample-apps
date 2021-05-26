# Peaky-Finders

[*Original published here*](https://github.com/kbaranko/peaky-finders)

Peaky Finders is a Plotly Dash application with helpful peak load visualizations and a day ahead forecasting model for five different ISOs. It does not demonstrate cutting-edge peak load forecasting methods -- there are a handful of high tech companies and millions of dollars spent trying to solve this problem -- but rather illustrate core concepts and explore how well a model can do with just historical load and temperature data.

The application has been deployed on Gallery: https://dash-gallery.plotly.host/dash-peaky-finders

## Tech Stack

- Python 
- Pandas
- Matplotlib
- Scikit-Learn
- Dash 
- Plotly

## Data

Historical load data was collected using the Pyiso python library, which provides clean API interfaces to make scraping ISO websites easy. The Darksky API was used for weather data, which provides historical temperature readings for a given latitude and longitude. For this model, I picked one central coordinate in each ISO territory to make API requests.

## Features

- Day of week (seven days)
- Holiday (yes or no)
- Hour of Day (24 hours)
- Temperature Reading (hourly)
- Previous Day’s Load (t-24)

## Results 

How well does each model perform? Depends on the ISO. Mean Absolute Error (MAE) for the month of February 2021 in Megawatts (MW):

- CAISO: 455.91
- MISO: 2,382.66 
- PJM: 2,886.66
- NYISO: 347.62
- ISONE: 522.43


