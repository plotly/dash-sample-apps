
# Omnisci Test DB


## Build docker image

The docker image can be built by running the following command from the project directory:
```
docker build -t omnisc_db .
```

To push the image to Google Container Registry (GCR):

```
docker tag omnisci_db gcr.io/plotly-hosting/omnisci_db
```

Get kubernetes credentials and push to GCR:
```
gcloud auth configure-docker
```

Followed by:

```
docker push gcr.io/plotly-hosting/omnisci_db
```

## Running locally:

The container can be run locally using the following command:

```
docker run -d -v $HOME/omnisci-docker-storage:/omnisci-storage -p 6273-6280:6273-6280 omnisci_db
```

## Deploying to GKE:

Deployment requires kubectl configurations. 

To Create GKE Cluster:

1. Create the cluster (if it is not present already):

```
gcloud container clusters create omnisci-db  --zone us-central1-b --machine-type n1-highcpu-2 --num-nodes 2 --enable-autoupgrade --project plotly-hosting
```

2. Deploy:

```
kubectl apply -f deployment/
```
