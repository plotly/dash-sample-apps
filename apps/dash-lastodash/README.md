# LAStoDash

'LAStoDash` is sample [Dash](https://plot.ly/products/dash) project that takes a
[Log ASCII Standard (LAS) file](http://www.cwls.org/las/) and builds a web app
to view its content and ready for printing.

## Installation

```
$ git clone https://github.com/n-riesco/lastodash.git
$ cd lastodash
$ pip3 install -r requirements.txt
```

## Usage

```
$ ./lastodash.py -h
usage: lastodash.py [-h] [--debug] lasfile

Launch a Dash app to view a LAS log.

positional arguments:
  lasfile      Log ASCII Standard (LAS) file

optional arguments:
  -h, --help   show this help message and exit
  --debug, -d  enable debug mode
```

## Example

```
$ ./lastodash.py alcor1.las 
Header section Parameter regexp=~P was not found.
 * Serving Flask app "lastodash" (lazy loading)
 * Environment: production
   WARNING: Do not use the development server in a production environment.
   Use a production WSGI server instead.
 * Debug mode: off
 * Running on http://127.0.0.1:8050/ (Press CTRL+C to quit)
```

And open http://127.0.0.1:8050/ to view the LAS file:

![Screencast](alcor1.gif)

See [here](alcor1.pdf) the report printed in PDF format.
