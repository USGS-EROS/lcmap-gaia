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

## Running a local Gaia

Start backing services, and add seed data 

```
make deps-up-d
```


## Running tests

```
make runtests
```

## Jupyter Notebook with Clojure kernel
```
(require '[cemerick.pomegranate :as pom])
(pom/add-classpath "/workspace/lcmap-gaia/target/gaia-0.1.0-SNAPSHOT-standalone.jar")
```
