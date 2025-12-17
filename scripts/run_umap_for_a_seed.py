from tqdm import tqdm
import umap
import pandas as pd
import argparse
from pathlib import Path
import os
import warnings
import concurrent.futures

warnings.filterwarnings("ignore", category=UserWarning)

# Create the argument parser
argparser = argparse.ArgumentParser()
argparser.add_argument("--rmd_input")
argparser.add_argument("--param_table")
argparser.add_argument("--seed", type=int)
argparser.add_argument("--output_dir")
argparser.add_argument("--n_process", type=int, default=1)

# Parse the arguments
args = argparser.parse_args()
rmd_input = args.rmd_input
param_table = args.param_table
seed = args.seed
output_dir = args.output_dir
n_process = int(args.n_process)

# Read input files and parameters
pca_proj = pd.read_csv(f"{rmd_input}/pca.csv", index_col=0)
params_df = pd.read_csv(param_table, index_col=0)
params_df = params_df.sample(frac=1)
params_df = params_df.reset_index(names="param_set")
params_dict = params_df.to_dict(orient="records")
output_dir = Path(output_dir)
os.makedirs(output_dir, exist_ok=True)

# Run UMAP for each param set and save the results
data = pca_proj
def run_umap(seed, data, param, output_dir):
    values = data.values
    param_set = param.pop("param_set")
    reducer = umap.UMAP(**param, random_state=seed, n_jobs=1)
    embedding = reducer.fit_transform(values)
    embedding_df = pd.DataFrame(embedding, index=data.index, columns=["UMAP1", "UMAP2"])
    embedding_df.to_csv(output_dir / f"param_set{param_set}.csv")
    return embedding

if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_process) as executor:
        futures = [
            executor.submit(run_umap, seed, data, param, output_dir)
            for param in params_dict
        ]

        for future in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
        ):
            future.result()
