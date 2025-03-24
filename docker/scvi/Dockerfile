FROM nciccbr/ccbr_ubuntu_base_20.04:v5

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# install miniconda
ENV CONDA_DIR /opt2/conda
WORKDIR /opt2
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt2/miniconda.sh && \
  /bin/bash /opt2/miniconda.sh -b -p /opt2/conda
ENV PATH=$CONDA_DIR/bin:$PATH

# install environment
RUN conda create -n scvi python=3.9
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate scvi && \
  cconda install scvi-tools -c conda-forge
RUN . /opt2/conda/etc/profile.d/conda.sh && \
  conda activate scvi && \
  conda install -c conda-forge r-base r-essentials r-reticulate
WORKDIR /opt2

# create wrapper for scvi conda env
COPY scvi /opt2
RUN chmod a+x /opt2/scvi
ENV PATH=/opt2:$PATH

# cleanup etc
# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod -R a+rx /opt2/argparse.bash
ENV PATH="/opt2/:$PATH"
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

COPY Dockerfile /opt2/Dockerfile
RUN chmod -R a+rX /opt2/Dockerfile
WORKDIR /data2
