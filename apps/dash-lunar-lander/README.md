# Dash Lunar Lander - Mass Optimal Trajectory

This is a demo of the Dash interactive Python framework developed by [Plotly](https://plot.ly/).

Dash abstracts away all of the technologies and protocols required to build an interactive web-based application and is a simple and effective way to bind a user interface around your Python code. To learn more check out our [documentation](https://plot.ly/dash).

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
cd dash-sample-apps/apps/dash-lunar-lander
pip install -r requirements.txt

```

Run the app

```

python app.py

```

## About the app

This app allows for exploring the trade space of a lunar lander. You can adjust different parameters on the lander, 
like max angular velocity and inital mass, using the sliders and a new optimal trajectory will be recomputed in near 
real time using Casadi.

The lunar lander model was based off of one in [this paper](https://arxiv.org/pdf/1610.08668.pdf).

The lunar lander model does not allowed for free rotation, instead requiring the lander to gimbal it’s 
engine and thrust to impart a rotation, giving the mass optimal control case a distinctive hooked shape. 
Switching the optimizer to time optimal control results in a much smoother and more expected shape.

## Built With

- [Dash](https://dash.plot.ly/) - Main server and interactive components
- [Plotly Python](https://plot.ly/python/) - Used to create the interactive plots
