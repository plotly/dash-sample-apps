# dash-clinical-analytics

## About this app

This dashboard allows you to explore clinic patient volume by time of the day, waiting time and care score. Drill down into high-volume days to view how the busy environment may have affected patient experience.

## Requirements

* Python 3

## How to run this app

We suggest you to create a virtual environment for running this app with Python 3. Clone this repository 
and open your terminal/command prompt in the root folder.

```
git clone https://github.com/plotly/dash-sample-apps
cd dash-sample-apps/apps/dash-clinical-analytics
python3 -m virtualenv venv

```
In Unix system:
```
source venv/bin/activate

```
In Windows: 

```
venv\Scripts\activate
```

Install all required packages by running:
```
pip install -r requirements.txt
```

Run this app locally with:
```
python app.py
```

## Screenshot

![screenshot](img/screencapture.png)

## Resources

* [Dash](https://dash.plot.ly/)
* Inspired by [Tableau template](https://www.tableau.com/solutions/workbook/improve-patient-satisfaction-improving-cycle-time).
