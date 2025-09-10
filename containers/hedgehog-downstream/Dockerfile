FROM r-base:4.4.2

## ------------ ##
## How to build ##
## ------------ ##

## sudo docker build --file containers/eland-downstream/Dockerfile --no-cache -t rtpop/eland-downstream:tag --progress=plain . 2>&1 | tee containers/eland-downstream/logs/eland-downstream:tag

## --------------------------- ##
## Install system requirements ##
## --------------------------- ##

RUN mkdir -p /opt/software
WORKDIR /opt

RUN apt-get update -qq && \
        apt-get install -y --allow-downgrades \
                apt-utils \
                libxml2-dev \
                libjpeg-dev \
                libpng-dev \
                libtiff-dev \
                libicu-dev \
                libgfortran5 \
                libnlopt-dev

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
                libssl-dev \
                vim

## ------------------- ##
## Copy required files ##
## ------------------- ##

COPY containers/eland-downstream/requirements.R /opt/requirements.R

## ---------------------------- ##
## Install package requirements ##
## ---------------------------- ##

RUN Rscript /opt/requirements.R

## --------------- ##
## Switch to shell ##
## --------------- ##

ENV HOME=/home

CMD ["/bin/bash"]