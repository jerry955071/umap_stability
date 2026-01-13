# This script aligns multiple embeddings using Procrustes analysis.
from scipy.linalg import orthogonal_procrustes
from itertools import combinations
import pandas as pd
import numpy as np
import argparse


def align_to_reference(Y, X):
    """
    Align embedding Y to reference X using orthogonal Procrustes.
    Returns aligned Y and RMSD.
    """
    Xc = X - X.mean(axis=0)
    Yc = Y - Y.mean(axis=0)

    R, _ = orthogonal_procrustes(Yc, Xc)
    Y_aligned = Yc @ R

    rmsd = np.sqrt(((Xc - Y_aligned) ** 2).sum(axis=1).mean())
    return Y_aligned, rmsd

def find_reference(embeddings):
    """
    embeddings: list of (n_samples, n_dims) numpy arrays
    Returns index of reference embedding and per-embedding mean RMSD.
    """
    n = len(embeddings)
    mean_rmsd = np.zeros(n)

    for i in range(n):
        Xi = embeddings[i]
        rmsds = []

        for j in range(n):
            if i == j:
                continue
            _, rmsd = align_to_reference(embeddings[j], Xi)
            rmsds.append(rmsd)

        mean_rmsd[i] = np.mean(rmsds)

    ref_idx = np.argmin(mean_rmsd)
    return ref_idx, mean_rmsd


# define parser
# parser = argparse.ArgumentParser(description="Align embeddings using Procrustes analysis.")
# parser.add_argument("--embeddings", type=str, help="Path to the source embeddings files.", nargs="+")
# parser.add_argument("--cluster", type=str, help="Path to cell clsuter csv.", required=True)
# parser.add_argument("--output", type=str, help="Path to save the aligned embeddings.", required=True)
# args = parser.parse_args()

# import embeddings
# embeddings_df = [pd.read_csv(path, index_col=0) for path in args.embeddings]
with open("/home/woodydrylab/DiskArray/b05b01002/project_umap_re/outputs/UMAP/seeds.txt", "r") as fin:
    seeds = [i.strip() for i in fin][:9]

embeddings_df = [
    pd.read_csv(
        f"/home/woodydrylab/DiskArray/b05b01002/project_umap_re/outputs/UMAP/ptr/seed{seed}/param_set97.csv",
        index_col=0
    ) for seed in seeds
] # for debugging
embeddings = [df.values for df in embeddings_df]

cluster_df = pd.read_csv(
    "/home/woodydrylab/DiskArray/b05b01002/project_umap_re/outputs/Seurat/ptr/clusters.csv", 
    index_col=0,
    header=0,
    names=["cluster"],
    dtype={"cluster": "str"}
) # for debuggin

# 
ref_idx, mean_rmsd = find_reference(embeddings)

print("Mean RMSD per seed:")
for seed, val in zip(seeds, mean_rmsd):
    print(seed, val)

print("Chosen reference seed:", seeds[ref_idx])

# 
X_ref = embeddings[ref_idx]
aligned_embeddings = []

for Y in embeddings:
    Y_aligned, _ = align_to_reference(Y, X_ref)
    aligned_embeddings.append(Y_aligned)

# 
aligned_dfs = []
for df, Y_aligned in zip(embeddings_df, aligned_embeddings):
    aligned_df = df.copy()
    aligned_df[["UMAP1", "UMAP2"]] = Y_aligned
    aligned_dfs.append(aligned_df)


# 
import matplotlib.pyplot as plt
import seaborn as sns
fig, axes = plt.subplots(3, 3, figsize=(10, 10))
for idx in range(len(aligned_dfs)):
    sns.scatterplot(
        data=pd.merge(aligned_dfs[idx], cluster_df, left_index=True, right_index=True),
        x="UMAP1",
        y="UMAP2",
        hue="cluster",
        ax=axes.flatten()[idx],
        legend=False,
        s=1
    )
    
    axes.flatten()[idx].set_title(f"seed={seeds[idx]}")
    axes.flatten()[idx].set_xlabel(None)
    axes.flatten()[idx].set_ylabel(None)


for idx in range(len(aligned_dfs)):
    aligned_dfs[idx]["seed"] = seeds[idx]
    aligned_dfs[idx] = pd.merge(
        aligned_dfs[idx],
        cluster_df,
        left_index=True,
        right_index=True
    )

concated_df = pd.DataFrame()
for idx in range(len(aligned_dfs)):
    concated_df = pd.concat([concated_df, aligned_dfs[idx]])

median_embedding = (
    concated_df
    .groupby(concated_df.index)
    .agg(
        UMAP1_median=("UMAP1", "median"),
        UMAP2_median=("UMAP2", "median"),
        cluster=("cluster", "first"),
        seed=("seed", "first")
    )
)
sns.scatterplot(
    data=median_embedding,
    x="UMAP1_median",
    y="UMAP2_median",
    hue="cluster",
    legend=False,
    s=10
)
