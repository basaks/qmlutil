#
# Ubuntu Dockerfile
#
# https://github.com/dockerfile/ubuntu
#

# Pull base image.
FROM ubuntu:16.04

# Install.
RUN \
  sed -i 's/# \(.*multiverse$\)/\1/g' /etc/apt/sources.list && \
  apt update && \
  apt -y upgrade && \
  apt install -y build-essential && \
  apt install -y software-properties-common && \
  apt install -y curl git htop man nano unzip vim wget && \
  apt install -y libblas-dev liblapack-dev libatlas-dev && \
  apt install -y libatlas-base-dev gfortran && \
  apt install -y python-pip && \\
  pip install -U pip numpy && \\
  pip install -U scipy obspy pytest

# Add files.
RUN git clone https://github.com/basaks/qmlutil.git
WORKDIR qmlutil
RUN python setup.py install

# Set environment variables.
ENV HOME /root

# Define working directory.
WORKDIR /root

# Define default command.
CMD ["bash"]



