<!--
To get started, replace
{{app.name}} with your app name (e.g. Dash Super Cool App)
{{app.short}} with the short handle (e.g. dash-super-cool)

If this is in dash sample apps, uncomment the second "git clone https..." and remove the first one.
If this is in dash sample apps and you have a colab demo, uncomment the "Open in Colab" link to see the badge (make sure to create a ColabDemo.ipynb) first.

-->
# {{app.name}}
<!-- 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/plotly/dash-sample-apps/blob/master/apps/{{app.short}}/ColabDemo.ipynb)
 -->

## Instructions

To get started, first clone this repo:

```
git clone https://github.com/plotly/{{app.short}}.git
cd {{app.short}}
```

<!--
```
git clone https://github.com/plotly/dash-sample-apps.git
cd dash-sample-apps/apps/{{app.short}}
```
-->

Create and activate a conda env:
```
conda create -n {{app.short}} python=3.7.6
conda activate {{app.short}}
```

Or a venv (make sure your `python3` is 3.6+):
```
python3 -m venv venv
source venv/bin/activate  # for Windows, use venv\Scripts\activate.bat
```

Install all the requirements:

```
pip install -r requirements.txt
```

You can now run the app:
```
python app.py
```

and visit http://127.0.0.1:8050/.

## Contact

Interested in building or deploying apps like this? [Reach out](https://plotly.com/contact-us/) or [get a demo](https://plotly.com/get-demo).
