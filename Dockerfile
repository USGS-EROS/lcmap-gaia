from nginx:latest
MAINTAINER USGS LCMAP http://eros.usgs.gov

RUN apt-get update
RUN apt-get install default-jdk curl vim -y

COPY startup.sh /startup.sh
COPY project.clj /project.clj
COPY target/lcmap-gaia-*-standalone.jar /

CMD java -jar lcmap-gaia-*-standalone.jar
