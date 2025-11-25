#!/bin/bash
set -euxo pipefail

# docker run
docker run \
    -d \
    -p 8792:8787 \
    -e USERID=1102 \
    -e GROUPID=1102 \
    -e PASSWORD=0000 \
    -e ROOT=true \
    --name rstudio_umap \
    chiaenu/rstudio_seurat_spacexr:2.2.1 \
    /init