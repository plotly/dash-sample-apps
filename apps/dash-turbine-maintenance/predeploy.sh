pip install --upgrade pip # to the latest version
pip install -r requirements-predeploy.txt  # Now, we can install them
if [ -f "predeploy.py" ]; then
    python predeploy.py
fi