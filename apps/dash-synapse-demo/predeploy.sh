curl https://packages.microsoft.com/ubuntu/18.04/prod/pool/main/m/msodbcsql17/msodbcsql17_17.3.1.1-1_amd64.deb --output msodbcsql17_17.3.1.1-1_amd64.deb

ACCEPT_EULA=Y dpkg -i msodbcsql17_17.3.1.1-1_amd64.deb

if [ -f "predeploy.py" ]; then
    python predeploy.py
fi