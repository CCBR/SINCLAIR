FROM nciccbr/ccbr_ubuntu_base_20.04:v7

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# R, python, and other packages available in conda should be installed with mamba here.
# Specify all dependencies in environment.yml
COPY environment.yml /data2/
COPY install.R /data2/
ENV CONDA_ENV=scvi
RUN mamba env create -n ${CONDA_ENV} -f /data2/environment.yml && \
    echo "conda activate ${CONDA_ENV}" > ~/.bashrc
ENV PATH="/opt2/conda/envs/${CONDA_ENV}/bin:$PATH"
RUN Rscript install.R

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
