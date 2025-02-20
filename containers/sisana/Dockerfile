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
        libnss-unknown \
        librsvg2-bin \
        librsvg2-dev \
        libudunits2-dev \
        littler \
        make \
        pkg-config \
        procps \
        libssl-dev \
        libopenmpi-dev \
        libgsl-dev \
        libcurl4-gnutls-dev \
        libmysqlclient-dev \
        libpq-dev \
        r-cran-devtools \
        r-cran-RcppParallel

                
RUN apt-get update -qq && \
        apt-get install -y --allow-downgrades \
                build-essential \
                git \
                libcurl4-openssl-dev \
                libcairo2-dev \
                littler \
                libssl-dev \
                libgsl-dev=2.8+dfsg-5

## --------------- ##
## Switch to shell ##
## --------------- ##

ENV HOME=/home

CMD ["/bin/bash"]