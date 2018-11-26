TAG:=`head -n 1 project.clj | grep -o '[0-9]*\.[0-9]*\.[0-9]*'`
IMAGE:=usgseros/lcmap-gaia

docker-build:
	lein uberjar
	docker build -t $(IMAGE):$(TAG) .

docker-push:
	docker push $(IMAGE):$(TAG)




