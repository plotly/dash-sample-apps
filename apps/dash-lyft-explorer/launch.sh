wget -O data.zip -Nq https://github.com/plotly/dash-lyft-perception-data/archive/master.zip
unzip -nq ./data.zip
mv dash-lyft-perception-data-master ./data
rm ./data.zip
gunicorn app:server --workers 4