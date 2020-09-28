#!/bin/bash


[ -d "venv" ] && (echo "venv already exisits, aborting..." && exit 1)
[ "$VIRTUAL_ENV" == "" ] || (echo "virtual env already active, aborting..." && exit 2)

# create virtual environment
[ -d "venv" ] || virtualenv --prompt '(`basename $PWD`/`basename "$VIRTUAL_ENV"`)' venv

source venv/bin/activate

# install requirements
[ "$VIRTUAL_ENV" == "" ] \
    && echo "Error starting virtual environment." \
    || (pip install -r requirements.txt && cat <<DOC

setup successful
don't forget to 
  source venv/bin/activate
before getting to work :-)
DOC
)
