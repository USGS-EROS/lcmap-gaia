from ubuntu:latest
MAINTAINER USGS LCMAP http://eros.usgs.gov

RUN apt-get update
RUN apt-get install default-jdk curl vim -y
COPY target/faux-0.1.0-SNAPSHOT-standalone.jar /
COPY cx-2115585_cy3119805_prediction_with_fake_date.json /
COPY cx-2115585_cy3119805_segment.json /
# java.xml.bind issue related to http-kit and java 9 https://github.com/http-kit/http-kit/issues/356
CMD java --add-modules java.xml.bind -jar faux-0.1.0-SNAPSHOT-standalone.jar
