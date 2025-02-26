FROM python:3.13

## ------------ ##
## How to build ##
## ------------ ##
# sudo docker build --file containers/eland-run/Dockerfile --no-cache -t rtpop/eland-run:tag --progress=plain . 2>&1 | tee containers/eland-run/logs/eland-run:tag

## --------------------------- ##
## Install system requirements ##
## --------------------------- ##

RUN mkdir -p /opt/software
WORKDIR /opt

RUN apt-get update -qq && \
        apt-get install -y --allow-downgrades \
                apt-utils

RUN apt-cache policy \
        apt-utils \
        build-essential \
        freetype2-demos \
        freetype2-doc \
        git \
        libcurl4-openssl-dev \
        libcairo2-dev \
        libfontconfig1-dev \
        libfreetype-dev \
        libfreetype6 \
        librsvg2-bin \
        librsvg2-dev \
        libudunits2-dev \
        make \
        pkg-config \
        procps \
        libssl-dev \
        libopenmpi-dev \
        libgsl-dev \
        libcurl4-gnutls-dev \
        libmysqlclient-dev \
        libpq-dev

                
RUN apt-get update -qq && \
        apt-get install -y --allow-downgrades \
                build-essential \
                cmake \
                git \
                libcurl4-openssl-dev \
                libcairo2-dev \
                libssl-dev \
                vim

## ------------ ##
## Pip installs ##
## ------------ ##
                
RUN pip install --upgrade pip \
        setuptools

## ---------------##
## Install sisana ##
## ---------------##

# try installing from github until he pushes fixes to pypi
WORKDIR /opt/software
RUN git clone https://github.com/kuijjerlab/sisana.git
WORKDIR /opt/software/sisana
RUN pip install .

#RUN pip install sisana

## ------------- ##
## Install eland ##
## ------------- ##

WORKDIR /opt/software
RUN git clone https://github.com/rtpop/ELAND.git
WORKDIR /opt/software/ELAND
RUN pip install .

## --------------- ##
## Switch to shell ##
## --------------- ##

WORKDIR ~
ENV HOME=/home

CMD ["/bin/bash"]