FROM nciccbr/ccbr_ubuntu_base_20.04:v7

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}

# install necessary library
COPY install_packages_v1.0.R /opt2
COPY Rpack_v1.0.yml /opt2
WORKDIR /opt2
RUN Rscript install_packages_v1.0.R

# cleanup etc
# Save Dockerfile in the docker
COPY Dockerfile /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
RUN chmod a+r /opt2/Dockerfile_${REPONAME}.${BUILD_TAG}
#COPY argparse.bash /opt2
RUN chmod -R a+rx /opt2/argparse.bash
ENV PATH="/opt2/:$PATH"
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
WORKDIR /data2
