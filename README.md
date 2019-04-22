[![Build Status](https://travis-ci.org/USGS-EROS/lcmap-gaia.svg?branch=develop)](https://travis-ci.org/USGS-EROS/lcmap-gaia)

# Gaia

Gaia calculates Change and Classification product values from CCDC results.
Use it to generate change and classification product values, and to produce maps.

Results are persisted to the configured Object Storage service.


## Deploying Gaia

Gaia is run as a Docker container. 

```
export NEMO_HOST="http://awesomehost.org/nemo"
export CHIPMUNK_HOST="http://awesomehost.org/chipmunk"
export CHIPMUNK_ACQUIRED="1999-01-01/2002-01-01"
export SEGMENTS_PATH="/conus_segments"
export PREDICTIONS_PATH="/conus_predictions"
export STORAGE_ENDPOINT="http://localhost:7480"
export STORAGE_ACCESS_KEY="9876asdrd"
export STORAGE_SECRET_KEY="13235lkjis"
export STORAGE_BUCKET="some-ceph-bucket"
export REGION="cu"
export CCD_VERSION="v01"
export HTTP_PORT=9876
export Xmx=4352m
export Xms=4352m

docker run -p 9876:${HTTP_PORT} -e NEMO_HOST=${NEMO_HOST} \
                                -e CHIPMUNK_HOST=${CHIPMUNK_HOST} \
                                -e CHIPMUNK_ACQUIRED=${CHIPMUNK_ACQUIRED} \
                                -e SEGMENTS_PATH=${SEGMENTS_PATH} \
                                -e PREDICTIONS_PATH=${PREDICTIONS_PATH} \
                                -e STORAGE_ENDPOINT=$(STORAGE_ENDPOINT) \
                                -e STORAGE_ACCESS_KEY=$(STORAGE_ACCESS_KEY) \
                                -e STORAGE_SECRET_KEY=$(STORAGE_SECRET_KEY) \
                                -e STORAGE_BUCKET=$(STORAGE_BUCKET) \
                                -e HTTP_PORT=${HTTP_PORT} \
                                -e REGION=${REGION} \
                                -e CCD_VERSION=${CCD_VERSION} \
                                -e Xmx=${Xmx} \
                                -e Xms=${Xms} \
                                -it usgseros/lcmap-gaia:latest
```

Gaia is configured using these environment variables:

| ENV                  | Description                            |
|----------------------|----------------------------------------|
| `NEMO_HOST`          | base url for lcmap-nemo resource       |
| `CHIPMUNK_HOST`      | base url for lcmap-chipmunk resource   |
| `CHIPMUNK_ACQUIRED`  | acquired value for requesting aux data |
| `SEGMENTS_PATH`      | resource path for segments data        |
| `PREDICTIONS_PATH`   | resource path for prediction data      |
| `STORAGE_ENDPOINT`   | url for object storage service         |
| `STORAGE_ACCESS_KEY` | access key for object storage service  |
| `STORAGE_SECRET_KEY` | secret key for object storage service  |
| `STORAGE_BUCKET`     | name of the object store bucket to use | 
| `REGION`             | region abbreviation (cu, ak, hi)       |
| `CCD_VERSION`        | version of ccd algorithm used to       |
|                      | generate input data                    |
| `Xmx`                | maximum JVM memory                     |
| `Xms`                | minimum JVM memory                     |
| `HTTP_PORT`          | HTTP port to expose the server at      |

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
http POST 127.0.0.1:9876/product cx="1484415" \
                                 cy="2114805" \
                                 dates:="[\"2006-07-01\"]" \
                                 tile="003008" \
                                 product="time-since-change"
```

## Requesting a map using HTTPie
```
http POST 127.0.0.1:9876/raster date="2006-07-01" \
                                tile="003008" \
                                tilex="1484415" \
                                tiley="2114805" \
                                product="time-since-change" \
                                chips:="[{\"cx\":\"1484415\", \"cy\":\"2114805\"},...]"
```

## Jupyter Notebook with Clojure kernel
```
(require '[cemerick.pomegranate :as pom])
(pom/add-classpath "/workspace/lcmap-gaia/target/gaia-0.2.7-standalone.jar")
```
