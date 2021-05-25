# Dash Canvas OCR

## How to Use this app

Write inside Canvas with your mouse toggling, or tablet stylus, and click **Sign** to translate your handwritten digits via optical character recognition.
The result of text recognition output will be reflected on the page.


### Running the app locally
We suggest you to create a separate virtual environment running Python 3 for this app, and install all of the required dependencies there. Run in Terminal/Command Prompt:

```
git clone https://github.com/plotly/dash-sample-apps
cd dash-sample-apps/apps/dash-canvas-ocr
python3 -m virtualenv venv
```
In UNIX system: 

```
source venv/bin/activate
```
In Windows: 

```
venv\Scripts\activate
```
This app requires installation of `tesseract` on your local device, check [this guide](https://linuxhint.com/install-tesseract-ocr-linux/) to install it according to your operating system.
To install all of the required packages to your virtual environment, simply run:

```
pip install -r requirements.txt
```

and all of the required `pip` packages, will be installed, and the app will be able to run.

```
pythonb app.py
```

![screenshot](screenshot.png)


