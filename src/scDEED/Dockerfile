FROM satijalab/seurat:5.0.0

USER root

ENV DEBIAN_FRONTEND=noninteractive \
    R_REMOTES_UPGRADE=never

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        pkg-config \
        libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e '\
    options(repos = c(CRAN = "https://cloud.r-project.org")); \
    install.packages("remotes"); \
    stopifnot(requireNamespace("remotes", quietly = TRUE)) \
'

# For a reproducible production image, replace "main" with a tested commit SHA.
ARG SCDEED_REF=main

RUN set -eux; \
    mkdir -p /tmp/scDEED; \
    curl -fsSL \
        "https://codeload.github.com/JSB-UCLA/scDEED/tar.gz/${SCDEED_REF}" \
        | tar -xz -C /tmp/scDEED --strip-components=1; \
    \
    sed -Ei \
        's/[[:space:]]*\(==[[:space:]]*[^)]*\)//g' \
        /tmp/scDEED/DESCRIPTION; \
    \
    sed -Ei \
        '/Seurat::(GetAssayData|SetAssayData)/ \
         s/slot[[:space:]]*=[[:space:]]*slot/layer = slot/g' \
        /tmp/scDEED/R/scDEED.R; \
    \
    Rscript -e '\
        options(repos = c(CRAN = "https://cloud.r-project.org")); \
        remotes::install_local( \
            "/tmp/scDEED", \
            dependencies = NA, \
            upgrade = "never" \
        ) \
    '; \
    \
    Rscript -e '\
        library(Seurat); \
        library(scDEED); \
        message("Seurat: ", packageVersion("Seurat")); \
        message("scDEED: ", packageVersion("scDEED")) \
    '; \
    \
    rm -rf /tmp/scDEED

# setup Rstudio server
ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2022.02.3+492
ENV DEFAULT_USER=rstudio
ENV PANDOC_VERSION=default
RUN apt update
RUN /bin/sh -c /rocker_scripts/install_rstudio.sh # buildkit
RUN /bin/sh -c /rocker_scripts/install_pandoc.sh # buildkit

# 
EXPOSE 8787/tcp
CMD ["/init"]
