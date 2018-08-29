TAG:=`head -n 1 project.clj | grep -o '[0-9]*\.[0-9]*\.[0-9]*'`
IMAGE:=usgseros/lcmap-gaia

checkdeps:
	lein deps

runtests:
	lein test	

docker-build:
	lein uberjar
	docker build -t $(IMAGE):$(TAG) .

docker-push:
	docker push $(IMAGE):$(TAG)

deps-up-d:
	docker-compose -f resources/docker-compose.yml up -d cassandra
	sleep 20
	docker-compose -f resources/docker-compose.yml up -d chipmunk
	sleep 10
	bin/seed

deps-down:
	docker-compose -f resources/docker-compose.yml down

foo:
	export CHIPMUNK_HOST=http://127.0.0.1:5656



