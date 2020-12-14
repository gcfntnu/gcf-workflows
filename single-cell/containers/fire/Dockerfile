FROM ubuntu:18.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
# conda path
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y \
    wget \
    git \
    unzip \
    gzip \
    build-essential \
    g++>=4.8.4 \
    libboost-dev>=1.54.0 \
    sed  && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]

RUN conda update conda
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

RUN conda install --quiet --yes \
    python=2.7 \
    cython>=0.23.4 \
    numpy>=1.13.3 \
    pandas>=0.20.3 \
    statsmodels>=0.8.0 \
    scipy>=1.1.0 \
    matplotlib>=2.1.0 \
    scikit-learn>=0.19.1 \
    cmocean	

WORKDIR /opt
RUN git clone https://github.com/princethewinner/FiRE.git 
RUN cd FiRE && touch utils/__init__.py && bash ./INSTALL --py

ENV PYTHONPATH /opt/FiRE:$PYTHONPATH

