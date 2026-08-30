FROM bioconductor/bioconductor_docker:RELEASE_3_23-r-4.6.0@sha256:b10002b39efa30c3779ad839549806ebdbb29b3266f0d2428478b04426e55929

RUN apt-get update && \
    apt-get install -y --no-install-recommends jq && \
    rm -rf /var/lib/apt/lists/*

RUN R --vanilla -e 'BiocManager::install("methylKit", version = "3.23", ask = FALSE, update = FALSE)'
RUN curl --fail --location \
      https://bioconductor.org/packages/3.23/bioc/src/contrib/methylKit_1.38.0.tar.gz \
      --output /tmp/methylKit_1.38.0.tar.gz && \
    echo "008990fdf453e7ec66781ea6736b8f373dfae12d891c76ef3bf3456e97ff264e  /tmp/methylKit_1.38.0.tar.gz" | sha256sum --check - && \
    R CMD INSTALL /tmp/methylKit_1.38.0.tar.gz && \
    rm /tmp/methylKit_1.38.0.tar.gz

COPY check_environment.R /opt/gigas/check_environment.R
RUN Rscript --vanilla /opt/gigas/check_environment.R

WORKDIR /project-gigas_ploidy
