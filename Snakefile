include: "workflows/CellRanger.smk"
include: "workflows/Seurat.smk"
include: "workflows/UMAP.smk"
include: "workflows/InterClusterAngle.smk"

configfile: "config/config.json"

wildcard_constraints:
    not_lch="|".join([s["species"] for s in config["references"] if s["species"] != "lch"])

# Custom functions used by all workflows
from typing import List
def query(d:List[dict], k:str, v:str) -> dict:
    """Return the first dictionary in a list of dictionaries where the value of key k matches v."""
    return [x for x in d if x[k] == v][0]

def query_all(d:List[dict], k:str, v:str, k_out:str) -> List[str]:
    """Return a list of values from key k_out in a list of dictionaries where the value of key k matches v."""
    return [x[k_out] for x in d if x[k] == v]

def _assembly(wildcards):
    """Get genome assembly file path by sample name"""
    species = query(config["samples"], "name", wildcards.sample)["species"]
    return query(config["references"], "species", species)["assembly"]

# docker run
docker_run="docker run --rm -u $(id -u):$(id -g) -v $(pwd):/data -w /data"

# 
rule Seurat:
    input:
        expand(
            "outputs/Seurat/{sample}",
            sample=[s["name"] for s in config["samples"]]
        )