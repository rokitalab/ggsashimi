FROM rocker/tidyverse:4.4.0
LABEL maintainer="psullivan@childrensnational.org"
WORKDIR /rocker-build/

# update and upgrade packages
RUN apt-get -y update --fix-missing \
    && apt-get -y upgrade \
    && apt-get -y update \
    && apt-get -y autoremove

# install needed tools
RUN apt-get update \
    &&  apt-get install -y --no-install-recommends \
         ca-certificates \
         curl \
         gdebi-core \
         libcairo2-dev \
         libfontconfig1-dev \
         lsb-release \
         python3-pip \
         python-is-python3 \
         fuse \
         wget \
         tar \
         git \
         man-db \
         unzip \
         openjdk-17-jdk \
         vim

# install CAVATICA SBFS
RUN curl https://igor.sbgenomics.com/downloads/sbfs/install.sh -sSf | sudo sh

# install R
ARG R_VER=3.6.3
RUN export UBUNTU_VER=$(lsb_release -rs | tr -d '.') \
    && curl -O https://cdn.rstudio.com/r/ubuntu-${UBUNTU_VER}/pkgs/r-${R_VER}_1_amd64.deb \
    && DEBIAN_FRONTEND=noninteractve gdebi --non-interactive r-${R_VER}_1_amd64.deb

ENV PATH=/opt/R/${R_VER}/bin:$PATH

RUN echo 'options(repos = "https://cloud.r-project.org/")' > ~/.Rprofile

## Install R packages
ARG GGPLOT_VER=3.3.3
RUN R -e 'install.packages("remotes");' && \
    R -e 'remotes::install_version("rlang", version="0.4.10")' && \
    R -e 'remotes::install_version("vctrs", version="0.3.6")' && \
    R -e 'remotes::install_version("pillar", version="1.4.7")' && \
    R -e 'remotes::install_version("tibble", version="3.0.6")' && \
    R -e 'remotes::install_version("gtable", version="0.3.0")' && \
    R -e 'remotes::install_version("scales", version="1.1.1")' && \
    R -e 'remotes::install_version("ggplot2", version="'${GGPLOT_VER}'")' && \
    R -e 'remotes::install_cran(c("gridExtra", "data.table", "svglite"))'

# Install pysam
RUN pip3 install pysam


# install needed tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libblas3 \
        libgomp1 \
        liblapack3 \
        locales \
        python-is-python3

# set locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

# set environment variables
ARG R_VER=3.6.3
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    PATH=/opt/R/${R_VER}/bin:$PATH

# install samtools
RUN apt-get update && apt-get -y upgrade && \
        apt-get install -y build-essential wget \
                libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev && \
        apt-get clean && apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /usr/src

#Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2 && \
        tar jxf samtools-1.21.tar.bz2 && \
        rm samtools-1.21.tar.bz2 && \
        cd samtools-1.21 && \
        ./configure --prefix $(pwd) && \
        make

ENV PATH=${PATH}:/usr/src/samtools-1.21

# Copy ggsashimi in the docker image
# ADD ggsashimi.py /

# Run the container as an executable
# ENTRYPOINT ["/ggsashimi.py"]

