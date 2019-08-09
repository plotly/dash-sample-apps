# US opioid epidemic dataset and Dash app

Poison induced death data was downloaded from [CDC Wonder](dash_app_screencast.gif), using cause-of-death codes X40–X44 (unintentional), X60–X64 (suicide), X85 (homicide), or Y10–Y14 (undetermined intent).

[View the Dash app](https://dash-playground.plotly.host/dash-opioid-epidemic/)

## Getting Started

### Running the app locally

First create a virtual environment with conda or venv inside a temp folder, then activate it.

```
virtualenv venv

# Windows
venv\Scripts\activate
# Or Linux
source venv/bin/activate

```

Clone the git repo, then install the requirements with pip

```

git clone https://github.com/plotly/dash-sample-apps
cd dash-sample-apps/apps/dash-opioid-epidemic
pip install -r requirements.txt

```

Run the app

```

python app.py

```

![plotly-dash-screencast](assets/app_screencast.gif)

Dash app inspired by [this Tableau dashboard](https://www.cdc.gov/nchs/data-visualization/drug-poisoning-mortality/)
