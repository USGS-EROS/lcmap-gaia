# Gaia

Gaia provides Change and Classification product calculation from CCDC results.

Use it to retrieve the change and classification product values for a given area.

## Deploying Gaia

Gaia is run as a Docker container. 

```
docker run -p 9876:9876 -e NEMO_HOST=${NEMO_HOST} \
                        -it usgseros/lcmap-gaia:latest
```

Gaia is configured using these environment variables:

| ENV             | Description                          |
|-----------------|--------------------------------------|
| `NEMO_HOST`     | base url for lcmap-nemo resource     |


## Running a local Gaia

Start backing services, and add seed data 

```
make deps-up-d
```


## Running tests

```
make runtests
```
