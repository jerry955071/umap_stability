from tqdm import tqdm
import umap
import pandas as pd
import argparse
from pathlib import Path
import os
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# Create the argument parser
argparser = argparse.ArgumentParser()
argparser.add_argument("--rmd_input")
argparser.add_argument("--param_table")
argparser.add_argument("--param_set")
argparser.add_argument("--seed_list")
argparser.add_argument("--output_dir")

# Parse the arguments
args = argparser.parse_args()
rmd_input = args.rmd_input
param_table = args.param_table
param_set = int(args.param_set)
seed_list = args.seed_list
output_dir = args.output_dir

# Read input files and parameters
pca_proj = pd.read_csv(f"{rmd_input}/pca.csv", index_col=0)
params_df = pd.read_csv(param_table, index_col=0)
params = params_df.to_dict(orient="index")[param_set]
seeds = pd.read_csv(seed_list, header=None)[0].tolist()
output_dir = Path(output_dir)
os.makedirs(output_dir, exist_ok=True)

# Run UMAP for each seed and save the results
data = pca_proj.values
for seed in tqdm(seeds):
    reducer = umap.UMAP(**params, random_state=seed, n_jobs=1)
    embedding = reducer.fit_transform(data)
    embedding_df = pd.DataFrame(embedding, index=pca_proj.index, columns=["UMAP1", "UMAP2"])
    embedding_df.to_csv(output_dir / f"umap_seed{seed}.csv")