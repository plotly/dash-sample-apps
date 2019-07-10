# Dash DAQ IV Tracer

## Introduction
`dash-daq-iv-tracer` uses the graphic elements of Dash DAQ to create an interface for acquiring current-voltage I-V curves with a Keithley 2400 SourceMeter.
[Try this demo app](https://dash-gallery.plotly.host/dash-daq-iv-tracer) on the Dash Deployment Server], and read more about the code in [our blog post](https://www.dashdaq.io/build-an-i-v-curve-tracer-with-a-keithley-2400-sourcemeter-in-python).

### IV Curves
I-V curve is a good way to characterize electronic components (diode, transistor or solar cells) and extract their operating properties. It is widely used in electrical engineering and physics. 
The Keithley 2400 SourceMeter provides precision voltage and current sourcing as well as measurement. 

### dash-daq
[Dash DAQ](dash.plot.ly/dash-daq) is a data acquisition and control package built on top of Plotly's [Dash](https://plot.ly/products/dash/). 

![Animated1](img/Screencast.gif)
## Requirements
It is advisable	to create a separate virtual environment running Python 3 for the app and install all of the required packages there. To do so, run:

```
python3 -m virtualenv [your environment name]
```
In Linux: 

```
source [your environment name]/bin/activate
```
In Windows: 

```
[your environment name]\Scripts\activate
```

To install all of the required packages to this environment, simply run:

```
pip install -r requirements.txt
```

and all of the required `pip` packages, will be installed, and the app will be able to run.


## How to use the app

To control your SourceMeter, you need to set the `mock` attribute to `False` in the `app.py` file.

```
iv_generator = keithley_instruments.KT2400(mock_mode=False)
```

You can then run the app :

```
$ python app.py
```

If you already know the COM/GPIB port number, you can feed it to the SourceMeter class as:

```
iv_generator = keithley_instruments.KT2400(
  mock_mode=False,
  instr_port_name=[your instrument's COM/GPIB port]
)
```

Or you can enter it from the app display in your browser and click the button labelled 'Connect'.

## Run mock application

If you don't have the instrument connected to your computer but would still like to test the app you can run

```
$ python app_mock.py
```

You can also set the `mock` attribute to `True` in the `app.py` file.

Open http://0.0.0.0:8050/ in your browser, and a mock app interface will be displayed.

![initial](img/Screenshot_light.png)

View the dark theme layout of this app by clicking on the toggle in the header.

![initial](img/Screenshot_dark.png)

# Controls
* Sourcing toggle: Set source type as either voltage or current.
* Measure mode toggle: Set device in single measurement mode, or sweep mode.
* Knob: Adjust source value applied by turning the knob.
* Indicator: Indicator will light up when selected mode is active.
* Clear Graph button: Clear plot.
* LED Display: Measured value will be reflected on LED display.

Click the **LEARN MORE** button on the app to learn more about the interactions involved in this app.
For more detailed explanations please refer to the [blog post](https://www.dashdaq.io/build-an-i-v-curve-tracer-with-a-keithley-2400-sourcemeter-in-python)


## Resources
Manual of the Keithley [2400](http://research.physics.illinois.edu/bezryadin/labprotocol/Keithley2400Manual.pdf)
