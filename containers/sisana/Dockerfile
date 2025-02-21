FROM python:3.10.16

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
                git \
                libcurl4-openssl-dev \
                libcairo2-dev \
                libssl-dev

# install sisana
# RUN pip install sisana
# for now cloning the repo and installing it that way until bug fix gets pushed to pypi
WORKDIR /usr/local/bin
RUN git clone https://github.com/kuijjerlab/sisana.git
WORKDIR sisana
RUN pip install .
WORKDIR /opt

## --------------- ## 
## Switch to shell ##
## --------------- ##

ENV HOME=/home

CMD ["/bin/bash"]