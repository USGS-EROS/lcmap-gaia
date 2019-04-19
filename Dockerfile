FROM ubuntu:18.04
MAINTAINER USGS LCMAP http://eros.usgs.gov

RUN apt-get update
RUN apt-get install openjdk-8-jre curl vim -y
RUN apt-get install gdal-bin libgdal-dev libgdal-java libgdal20 python-gdal -y
COPY resources/log4j.properties /log4j.properties
COPY target/gaia-*-standalone.jar /
CMD java -server -Xms$Xms -Xmx$Xmx -XX:+UseG1GC -jar gaia-*-standalone.jar

