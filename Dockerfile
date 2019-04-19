FROM ubuntu:18.04
MAINTAINER USGS LCMAP http://eros.usgs.gov

RUN apt-get update
RUN apt-get install default-jdk curl vim -y
RUN apt-get install gdal-bin libgdal-dev libgdal-java libgdal20 python-gdal -y
COPY resources/log4j.properties /log4j.properties
COPY target/gaia-*-standalone.jar /
# java.xml.bind issue related to http-kit and java 9 https://github.com/http-kit/http-kit/issues/356
ENTRYPOINT java -server -Xms$Xms -Xmx$Xmx -XX:+UseG1GC --add-modules java.xml.bind -jar gaia-*-standalone.jar

