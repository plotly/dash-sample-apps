# Predictive Maintenance for Wind Turbines Dashboard

## Introduction

`dash-predictive-maintenance` is a dashboard designed to demonstrate the power of Machine Learning to predict failures (Remaining Useful Life (RUL)) in wind turbines. The data covers periods from May, 2014 to January, 2015. To predict the date when equipment will completely fail (RUL), XGBoost is used. The achieved RMSE error is `0.018534` days, which is highly accurate.

## Screenshots
![initial](screenshots/screenshot1.png)

![initial](screenshots/screenshot2.png)

## Built With
* [Dash](https://dash.plot.ly/) - Main server and interactive components.
* [Dash DAQ](https://dash.plot.ly/dash-daq) - Styled technical components for industrial applications.
* [XGBoost](https://xgboost.readthedocs.io/en/latest/) - Machine Learning model that was used to predict the RUL. The model was fine-tuned with RandomSearch.


## Requirements
Clone this repo and create a clean environment:
```
git clone https://github.com/iameminmammadov/dash-predictive-maintenance.git
cd dash-predictive-maintenance
python3 -m virtualenv venv
```
To activate the virtualenv in UNIX:
```
source venv/bin/activate
```
To activate the virtualenv in Windows:
```
venv\Scripts\activate
```
To install the libraries, needed to run this dashboard:
```
pip install -r requirements-predeploy.txt
```
To run this app:
```
python app.py
```
The app will be run on  http://127.0.0.1:8050/.
