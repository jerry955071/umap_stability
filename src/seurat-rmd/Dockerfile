FROM satijalab/seurat:5.0.0

# 
ENV S6_VERSION=v2.1.0.2
ENV RSTUDIO_VERSION=2022.02.3+492
ENV DEFAULT_USER=rstudio
ENV PANDOC_VERSION=default
RUN apt update
RUN /bin/sh -c /rocker_scripts/install_pandoc.sh # buildkit 
