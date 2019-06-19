# Dash Web Trader
This is a demo of the Dash interactive Python framework developed by [Plotly](https://plot.ly/).
                            
This app continually queries csv files and updates Ask and Bid prices for major currency 
pairs as well as Stock Charts. You can also virtually buy and sell stocks and see the 
profit updates.

## How to run this app

(The following instructions apply to Windows command line.)

To run this app first clone repository and then open a terminal to the app folder.

```
git clone https://github.com/plotly/dash-sample-apps.git
cd dash-sample-apps/apps/dash-web-trader
```

Create and activate a new virtual environment (recommended) by running
the following:

On Windows

```
virtualenv venv 
\venv\scripts\activate
```

Or if using linux

```bash
python3 -m venv myvenv
source myvenv/bin/activate
```

Install the requirements:

```
pip install -r requirements.txt
```
Run the app:

```
python app.py
```
You can run the app on your browser at http://127.0.0.1:8050

## Screenshots

![demo.png](demo.png)

## Resources

To learn more about Dash, please visit [documentation](https://plot.ly/dash).