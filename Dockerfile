# GCC support can be specified at major, minor, or micro version
# (e.g. 8, 8.2 or 8.2.0).
# See https://hub.docker.com/r/library/gcc/ for all supported GCC
# tags from Docker Hub.
# See https://docs.docker.com/samples/library/gcc/ for more on how to use this image
# FROM gcc:9.2
FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get -y update && apt-get install -y

RUN apt-get -y install cmake libgtest-dev libvtk7-dev libproj-dev libmedc-dev mpich libboost-test-dev && rm -rf /var/lib/apt/lists/* 

COPY . /usr/src/fluidchen

WORKDIR /usr/src/fluidchen

RUN cmake .

RUN make
