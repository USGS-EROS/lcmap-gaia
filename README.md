[![Build Status](https://travis-ci.org/USGS-EROS/lcmap-gaia.svg?branch=develop)](https://travis-ci.org/USGS-EROS/lcmap-gaia)

# Gaia

Gaia provides Change and Classification product calculation from CCDC results.

Use it to retrieve the change and classification product values for a given area.

## Deploying Gaia

Gaia is run as a Docker container. 

```
export NEMO_HOST="http://awesomehost.org/nemo"
export SEGMENTS_PATH="/conus_segments"
export PREDICTIONS_PATH="/conus_predictions"
export REGION="cu"
export CCD_VERSION="v01"
export HTTP_PORT=9876

docker run -p 9876:${HTTP_PORT} -e NEMO_HOST=${NEMO_HOST} \
                                -e SEGMENTS_PATH=${SEGMENTS_PATH} \
                                -e PREDICTIONS_PATH=${PREDICTIONS_PATH} \
                                -e HTTP_PORT=${HTTP_PORT} \
                                -it usgseros/lcmap-gaia:latest
```

Gaia is configured using these environment variables:

| ENV                | Description                          |
|--------------------|--------------------------------------|
| `NEMO_HOST`        | base url for lcmap-nemo resource     |
| `SEGMENTS_PATH`    | resource path for segments data      |
| `PREDICTIONS_PATH` | resource path for prediction data    |
| `REGION`           | region abbreviation (cu, ak, hi)     |
| `CCD_VERSION`      | version of ccd algorithm used to     |
|                    | generate input data                  |


## Running a local Gaia

Using docker-compose, build a container to provide sample json, and a gaia instance

```
docker-compose -f resources/docker-compose.yml up

```

## Running tests

```
lein test
```

## Requesting products using HTTPie https://httpie.org
```
http GET  localhost:9876/product/time-since-change/-2115585/3119805/1996-07-01 Accept:application/json > time_since_change.json

```

## Jupyter Notebook with Clojure kernel
```
(require '[cemerick.pomegranate :as pom])
(pom/add-classpath "/workspace/lcmap-gaia/target/gaia-0.1.0-SNAPSHOT-standalone.jar")
```
