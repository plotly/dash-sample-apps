### Dash Financial Report
###### multi-page, external_css, external_scripts, PDF      

This is a demo of the [Dash](https://plot.ly/products/dash/) interactive Python framework developed by [Plotly](https://plot.ly/).

Dash abstracts away all of the technologies and protocols required to build an interactive web-based application and is a simple and effective way to bind a user interface around your Python code.

To learn more about Dash, take a look at our [documentation](https://dash.plot.ly). If you're interested in deploying this application, check out [Dash Deployment Server](https://dash.plot.ly/dash-deployment-server/) - Plotly's commercial offering for hosting and sharing Dash Apps on-premise or in the cloud. 

##### About this repo:

For more information about the application structure, see our [Dash Deployment Server Documentation](https://dash.plot.ly/dash-deployment-server/application-structure).

##### To run this app:

You can clone or download this repo:   
```
git clone https://github.com/plotly/dash-financial-report.git
```

Then cd into the repo:   
```
cd dash-financial-report
```

Now create and activate a virtualenv (noting the python runtime):   
On a mac:   
```
virtualenv -p <python version> venv
source venv/bin/activate
```

On a Windows:   
```
virtualenv -p <python version> venv
venv/Scripts/activate
```

Now that virtualenv is setup and active we can install the dependencies:   
```
pip install -r requirements.txt
```

Once the dependencies have been installed, run the application:
```
python app.py
```

Then visit http://127.0.0.1:8050/
