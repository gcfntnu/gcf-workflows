FROM ubuntu:18.04

# conda path
ENV PATH /opt/conda/bin:$PATH
# use C.UTF-8 locale 
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# R pre-requisites ++
RUN apt-get update --fix-missing && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    wget \
    ca-certificates \
    fonts-dejavu \
    build-essential \
    gfortran \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


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
    conda config --add channels r && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

RUN conda install --quiet --yes r-base r-argparse r-devtools r-biocmanager r-tidyverse r-ggplot2 && \
    conda clean -tipsy
 
RUN ln -s /bin/tar /bin/gtar
RUN R --slave -e "options(unzip='/usr/bin/unzip'); BiocManager::install('mikelove/tximport'); BiocManager::install('flatberg/alevinQC')"
