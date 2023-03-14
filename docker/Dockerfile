# Cellula image

################################################################

## Start Debian basic image
FROM python:3.9.16-bullseye

## Change as u want
LABEL mantainer = "andrea.cossa@ieo.it"

## Compilers, utilities, git
RUN apt-get --quiet update --yes \
    && apt-get --quiet upgrade --yes \
    && apt-get --quiet install --yes \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      build-essential \
      libz-dev \
      wget \
      procps \
      git \
      cmake \
    && rm -rf /var/lib/apt/lists/*

## Mamba
RUN cd / \
  && wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && bash Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda3 \
  && rm -f Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/miniconda3/bin/:$PATH"

# Mamba
RUN conda install -c conda-forge mamba
ENV PATH="/mambaforge/bin/:$PATH"

## Clone repo and install code
RUN echo changes env
RUN git clone https://github.com/andrecossa5/MI_TO.git
RUN cd MI_TO \
    && mamba env update -n base -f envs/environment.yml \
    && mamba develop .. \
    && mamba develop ../MI_TO/ 
 
## Default command to initialize the container
CMD [ "bash" ]

#################################################################

