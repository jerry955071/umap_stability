FROM mambaorg/micromamba:git-35c00b7-amazon2023

# install umap-learn
RUN micromamba install -y -n base conda-forge::umap-learn pandas tqdm && \
    micromamba clean -a -y
