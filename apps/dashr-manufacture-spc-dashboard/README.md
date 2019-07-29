# Manufacture SPC Dashboard

## About this app

Manufacture SPC Dashboard is a dashboard for monitoring real-time process quality along manufacture production line.

## Use this app

Click on buttons in `Parameter` column to visualize details of trendline on the bottom panel.

Click `Start` button, trends are updated every two seconds to simulate real-time measurements. The Sparkline on top panel and Control chart on bottom panel show Shewhart process control using mock data. Data falling outside of control limit are signals indicating 'Out of Control (OOC)', and will 
trigger alerts instantly for a detailed checkup. 

Operators may stop measurement by clicking `Stop` button, and edit specification parameters for selected process line(metrics) in Specification Tab.

![demo/demo.gif](demo/demo.gif)

## Run this app locally

Clone the repository:

```
$ git clone https://github.com/plotly/dash-sample-apps.git
```

Redirect to the respective app directory:

```
$ cd dash-sample-apps/apps/dashr-manufacture-spc-dashboard
```

Install the requirements:

```
$ Rscript init.R
```

Run the app:

```
$ Rscript app.R
```

View in your browser at http://127.0.0.1:8050.

## Resources/References

* [Dash R](https://dashr.plot.ly/) - Main server and interactive components 
* [Plotly R](https://plot.ly/r/) - Used to create the interactive plots
* [Dash DAQ](https://dashr.plot.ly/dash-daq) - Styled technical components for industrial applications
* [Shewhart statistical process control](https://en.wikipedia.org/wiki/Shewhart_individuals_control_chart)
