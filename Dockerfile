from ubuntu:latest
MAINTAINER USGS LCMAP http://eros.usgs.gov

RUN apt-get update
RUN apt-get install default-jdk curl vim -y
COPY target/gaia-*-standalone.jar /
# java.xml.bind issue related to http-kit and java 9 https://github.com/http-kit/http-kit/issues/356
CMD java --add-modules java.xml.bind -jar gaia-*-standalone.jar
