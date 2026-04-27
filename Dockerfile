FROM rocker/tidyverse:3.6.3
LABEL maintainer="psullivan@childrensnational.org"
WORKDIR /rocker-build/

RUN sed -i 's|deb.debian.org|archive.debian.org|g' /etc/apt/sources.list \
 && sed -i 's|security.debian.org|archive.debian.org|g' /etc/apt/sources.list \
 && sed -i '/buster-updates/d' /etc/apt/sources.list \
 && apt-get -o Acquire::Check-Valid-Until=false update


# update and upgrade packages
RUN apt-get -y update --fix-missing \
    && apt-get -y upgrade \
    && apt-get -y update \
    && apt-get -y autoremove

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
         ca-certificates \
         curl \
         gdebi-core \
         libcairo2-dev \
         libfontconfig1-dev \
         lsb-release \
         python3 \
         python3-pip \
         python3-dev \
	 fuse \
         wget \
         tar \
         git \
         man-db \
         unzip \
         openjdk-11-jdk \
         vim \
         libbz2-dev \
         zlib1g-dev \
         liblzma-dev \
         libcurl4-openssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install CAVATICA SBFS
RUN curl https://igor.sbgenomics.com/downloads/sbfs/install.sh -sSf | sudo sh

# Install system libs required by R packages
RUN apt-get update && apt-get install -y \
    build-essential \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    && rm -rf /var/lib/apt/lists/*

# Clean site library
RUN rm -rf /usr/local/lib/R/site-library/*

RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/2020-04-01'))" >> /usr/local/lib/R/etc/Rprofile.site

RUN R -e "install.packages(c( \
  'BiocManager', \
  'data.table', \
  'ggplot2', \
  'gridExtra', \
  'gtable', \
  'withr', \
  'gridExtra', \
  'openxlsx', \
  'optparse', \
  'pillar', \
  'rlang', \
  'scales', \
  'svglite', \
  'tibble', \
  'tidyverse', \
  'vctrs' \
), dependencies=TRUE, verbose=TRUE)"

RUN R -e 'BiocManager::install(c( \
  "AnnotationHub", \
  "biomaRt", \
  "BSgenome.Hsapiens.UCSC.hg38", \
  "data.table", \
  "ensembldb", \
  "GenomicFeatures", \
  "GenomicRanges", \
  "rtracklayer" \
))'

RUN R -e "library(ggplot2); library(withr); library(gridExtra); library(openxlsx); library(optparse)"

RUN python3 -m pip install --upgrade pip setuptools wheel

# Install pysam
RUN python3 -m pip install --no-cache-dir pysam

RUN python3 -c "import pysam; print(pysam.__version__)"

# install needed tools
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libblas3 \
        libgomp1 \
        liblapack3 \
        locales \
        python3

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
