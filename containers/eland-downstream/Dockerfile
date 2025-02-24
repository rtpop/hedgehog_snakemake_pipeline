FROM r-base:4.4.2

## build with 
## sudo docker build --file containers/eland-downstream/Dockerfile --no-cache -t rtpop/eland-downstream:0.0.1 --progress=plain . 2>&1 | tee eland-downstream:0.0.1

## --------------------------- ##
## Install system requirements ##
## --------------------------- ##

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
        littlr \
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

## ------------------- ##
## Copy required files ##
## ------------------- ##

COPY containers/eland-downstream/requirements.R ~/requirements.R

## ---------------------------- ##
## Install package requirements ##
## ---------------------------- ##

RUN Rscript requirements.R

## --------------- ##
## Switch to shell ##
## --------------- ##

ENV HOME=/home

CMD ["/bin/bash"]