from scipy.spatial.distance import cosine
from itertools import combinations
from tqdm import tqdm
import pandas as pd
import argparse
from pathlib import Path
import os
import warnings
import concurrent.futures
import random
import sys

warnings.filterwarnings("ignore", category=UserWarning)

# ----------------------------
# Globals (for worker processes)
# ----------------------------
CLUSTER_DF = None

def init_worker(cluster_table_path):
    global CLUSTER_DF
    CLUSTER_DF = pd.read_csv(cluster_table_path, index_col=0)
    CLUSTER_DF.columns = ["cluster"]

# ----------------------------
# Argument parsing
# ----------------------------
argparser = argparse.ArgumentParser()
argparser.add_argument("--umap_indir", required=True)
argparser.add_argument("--cluster_table", required=True)
argparser.add_argument("--param_table", required=True)
argparser.add_argument("--seed_list")
argparser.add_argument("--output_dir", required=True)
argparser.add_argument("--n_process", type=int, default=1)
args = argparser.parse_args()

umap_indir = Path(args.umap_indir)
cluster_table = Path(args.cluster_table)
param_table = Path(args.param_table)
output_dir = Path(args.output_dir)
n_process = int(args.n_process)

os.makedirs(output_dir, exist_ok=True)

# ----------------------------
# Seeds
# ----------------------------
seeds = [int(x) for x in open(args.seed_list).read().strip().split("\n")]

# ----------------------------
# Parameter sets
# ----------------------------
params_df = pd.read_csv(param_table, index_col=0)
params_df = params_df.sample(frac=1)
param_sets = params_df.index.tolist()

# ----------------------------
# Worker functions
# ----------------------------
def calculate_inter_cluster_angle(umap_indir, seed, param_set):
    global CLUSTER_DF

    embedding = pd.read_csv(
        umap_indir / f"seed{seed}" / f"param_set{param_set}.csv",
        index_col=0
    )

    merged_df = embedding.join(CLUSTER_DF, how="inner")
    cluster_centers = merged_df.groupby("cluster").median()

    results = []
    for vertex in cluster_centers.index:
        others = cluster_centers.index.drop(vertex)
        for a, b in combinations(others, 2):
            vec_a = cluster_centers.loc[a] - cluster_centers.loc[vertex]
            vec_b = cluster_centers.loc[b] - cluster_centers.loc[vertex]
            cos_similarity = 1 - cosine(vec_a, vec_b)
            results.append({
                "param_set": param_set,
                "seed": seed,
                "vertex": vertex,
                "cluster_a": a,
                "cluster_b": b,
                "cos_similarity": cos_similarity
            })
    return results


def run_one_seed(seed):
    rows = []
    for param_set in param_sets:
        rows.extend(
            calculate_inter_cluster_angle(
                umap_indir,
                seed,
                param_set
            )
        )
    pd.DataFrame(rows).to_csv(
        output_dir / f"seed{seed}.csv",
        index=False
    )
    return

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    all_results = []

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=n_process,
        initializer=init_worker,
        initargs=(cluster_table,)
    ) as executor:

        futures = [
            executor.submit(run_one_seed, seed)
            for seed in seeds
        ]

        for fut in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures)
        ):
            fut.result()
        


sys.exit(0)
