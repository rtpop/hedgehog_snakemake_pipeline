FROM python:3.10.16

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
        setuptools \
        community

## ---------------##
## Install sisana ##
## ---------------##

# try installing from github until he pushes fixes to pypi
WORKDIR /opt/software
RUN git clone https://github.com/kuijjerlab/sisana.git
WORKDIR /opt/software/sisana
RUN git checkout v1.2.0
RUN pip install .

#RUN pip install sisana

## ------------- ##
## Install eland ##
## ------------- ##

WORKDIR /opt/software
RUN git clone https://github.com/rtpop/ELAND.git
WORKDIR /opt/software/ELAND
RUN pip install .

## ----------------------------------------- ##
## Force older versions of numpy and scipy   ##
## Because otherwise hidef/bihidef complains ##
## Plus make sure we keep the current        ##
## version of netZooPy                       ##
## ----------------------------------------- ##

RUN pip install --no-cache-dir numpy==1.24.4 scipy==1.13.0 netZooPy==0.10.8

## -------------------------- ##
## Set matplotlib environment ##
## -------------------------- ##

# Set environment variables
ENV MPLCONFIGDIR=/tmp/matplotlib
ENV FONTCONFIG_PATH=/tmp/fontconfig
ENV FONTCONFIG_FILE=/tmp/fontconfig/fonts.conf

## --------------- ##
## Switch to shell ##
## --------------- ##

WORKDIR ~
ENV HOME=/home

CMD ["/bin/bash"]