FROM ubuntu:18.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update --fix-missing && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y \
  wget \
  ca-certificates \
  git \
  unzip \
  build-essential \
  flex \
  bison \
  autoconf \
  cmake \
  libboost-all-dev \
  libjemalloc-dev \
  zlib1g-dev \
  libbz2-dev \
  libcereal-dev \
  liblzma-dev \
  curl && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/COMBINE-lab/salmon.git --branch develop  && \
  mkdir -p salmon/build && \
  mkdir -p /opt/salmon/ && \
  cd salmon/build && \
  cmake -DCMAKE_INSTALL_PREFIX=/opt/salmon .. && \
  make && make install && \
  cd && \
  rm -rf salmon/
ENV PATH /opt/salmon/bin:$PATH
ENV LD_LIBRARY_PATH /opt/salmon/lib:$LD_LIBRARY_PATH

RUN \
  mkdir -p /opt/gffread && \
  cd /opt/gffread && \
  git clone https://github.com/gpertea/gclib && \
  git clone https://github.com/gpertea/gffread && \
  cd gffread && \
  make release
ENV PATH /opt/gffread/gffread:$PATH

RUN git clone https://github.com/lh3/bioawk.git && \
  cd bioawk && \
  make && \
  cp bioawk /usr/local/bin/ && \
  cd .. && rm -rf bioawk