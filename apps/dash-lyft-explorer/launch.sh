wget -O data.zip -Nq https://sampleappsdata.blob.core.windows.net/dash-sample-apps-data/dash-lyft-perception-data-master.zip
unzip -nq ./data.zip
mv dash-lyft-perception-data-master ./data
rm ./data.zip
gunicorn app:server --workers 4