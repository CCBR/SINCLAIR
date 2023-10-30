FROM nciccbr/ccbr_ubuntu_base_20.04:v5
FROM r-base:4.3.0

# NOTES
# https://github.com/satijalab/seurat-docker/blob/master/latest/Dockerfile
# https://github.com/CCBR/Dockers/blob/master/common/ccbr_r_4.3.0/Dockerfile.v1

# Set global R options
RUN echo "options(repos = 'http://cran.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# build time variables
ARG BUILD_DATE="000000"
ENV BUILD_DATE=${BUILD_DATE}
ARG BUILD_TAG="000000"
ENV BUILD_TAG=${BUILD_TAG}
ARG REPONAME="000000"
ENV REPONAME=${REPONAME}


# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev \
    pkg-config

#RUN apt-get install -y llvm-10

# Install system library for rgeos
RUN apt-get install -y libgeos-dev

# Install UMAP
#RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite
#RUN pip3 install numpy
#RUN pip3 install umap-learn

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

# Github
RUN R  --no-echo --no-restore --no-save -e 'install.packages("remotes",repos="http://cran.r-project.org",update=TRUE,ask=FALSE,dependencies=TRUE); \
    remotes::install_github("satijalab/seurat",version="seurat5",update=TRUE,ask=FALSE,dependencies=TRUE)'

# Install bioconductor dependencies & suggests
# RUN Rscript --no-echo --no-restore --no-save -e 'install.packages("BiocManager",repos="http://cran.r-project.org",update=TRUE,ask=FALSE,dependencies=TRUE); BiocManager::install(c("multtest", "S4Vectors", "SummarizedExperiment", "SingleCellExperiment", "MAST", "DESeq2", "BiocGenerics", "GenomicRanges", "IRanges", "rtracklayer", "monocle", "Biobase", "limma", "glmGamPoi"),repos="http://cran.r-project.org",update=TRUE,ask=FALSE,dependencies=TRUE)'

# Install CRAN suggests
# RUN R --no-echo --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools'))"

# Install spatstat
# RUN R --no-echo --no-restore --no-save -e "install.packages(c('spatstat.explore', 'spatstat.geom'))"

# Install hdf5r
# RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"

# Install latest Matrix
# RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix')"

# Install rgeos
# RUN R --no-echo --no-restore --no-save -e "install.packages('rgeos')"

# Install Seurat
# RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"

# install necessary library
WORKDIR /opt2
COPY install_packages_v1.0.R /opt2
COPY Rpack_v1.0.config /opt2
WORKDIR /opt2
RUN Rscript /opt2/install_packages_v1.0.R

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
