# Dash AlphaFold SARS-COV-2 Viewer

This app was cloned from @IvoLeist's excellent [Dash NGL viewer](https://github.com/IvoLeist/dash_ngl). Please take a look for more information.

This app lets you visualize the recently released *computational predictions of protein structures associated with COVID-19*, generated using AlphaFold by the DeepMind team. Please see the [official blog post](https://deepmind.com/research/open-source/computational-predictions-of-protein-structures-associated-with-COVID-19) for more details and citations.

## Instructions

To get started, first clone this repo:
<!-- 
```
git clone https://github.com/plotly/dash-alphafold.git
cd dash-alphafold
``` -->


```
git clone https://github.com/plotly/dash-sample-apps.git
cd dash-sample-apps/apps/dash-alphafold
```


Create and activate a conda env:
```
conda create -n dash-alphafold python=3.7.6
conda activate dash-alphafold
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
python usage.py
```

and visit http://127.0.0.1:8050/.

## Contact

Interested in building or deploying apps like this? [Reach out](https://plotly.com/contact-us/) or [get a demo](https://plotly.com/get-demo).