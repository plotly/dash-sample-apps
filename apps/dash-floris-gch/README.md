# Visualizing GCH with FLORIS

## Instructions

To get started, first clone this repo:

```
git clone https://github.com/plotly/floris-dash-app.git
cd floris-dash-app
```

### Usual setup

Create and activate a virtual environment (make sure your `python3` is 3.6+):
```
python3 -m venv venv
source venv/bin/activate
```

Install all the requirements:

```
pip3 install -r requirements.txt
```

You can now run the app:
```
python3 app.py
```

and visit http://127.0.0.1:8050/.

### Alternative setups

#### Windows
```
python3 -m venv venv
venv\Scripts\activate.bat
pip3 install -r requirements.txt
python3 app.py
```

#### Conda
```
conda create -n floris-dash-app python=3.7
conda activate floris-dash-app
pip3 install -r requirements.txt
python3 app.py
```

#### Pipenv
```
pipenv --python 3.7
pipenv install -r requirements.txt
```




## Contact

Interested in building or deploying apps like this? [Reach out](https://plotly.com/contact-us/) or [get a demo](https://plotly.com/get-demo).