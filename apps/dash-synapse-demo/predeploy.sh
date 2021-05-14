sudo su
curl https://packages.microsoft.com/keys/microsoft.asc | apt-key add -

#Ubuntu 18.04
curl https://packages.microsoft.com/config/ubuntu/18.04/prod.list > /etc/apt/sources.list.d/mssql-release.list

exit
sudo apt-get update
sudo ACCEPT_EULA=Y apt-get install -y msodbcsql17
# optional: for unixODBC development headers
sudo apt-get install -y unixodbc-dev

# install requirements that needed the commands above
pip install -r requirements-predeploy.txt