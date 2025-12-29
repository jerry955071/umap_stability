configfile: "config/config.json"

# number of seed
n_seeds = 1000
import random
random.seed(42)
seeds = random.sample(range(10**4, 10**5 -1), n_seeds)

# params to test
n_neighbors=[10, 30, 50, 100, 200]
min_dist=[0.0, 0.1, 0.3, 0.5, 0.8, 0.99]
metric=["euclidean", "cosine", "correlation"]
n_epochs=[100, 200, 500, 700]
n_combinations = len(n_neighbors) * len(min_dist) * len(metric) * len(n_epochs)
rule param_table:
    output:
        "outputs/UMAP/param_table.csv"
    params:
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        metric=metric,
        n_epochs=n_epochs
    run:
        import pandas as pd
        from itertools import product
        param_combinations = list(product(params.n_neighbors, params.min_dist, params.metric, params.n_epochs))
        df = pd.DataFrame(param_combinations, columns=["n_neighbors", "min_dist", "metric", "n_epochs"])
        df.index.name = "param_set"
        df.to_csv(output[0])

rule generate_seeds:
    output:
        "outputs/UMAP/seeds.txt"
    params:
        n_seeds=n_seeds
    log:
        "logs/UMAP/generate_seeds.log"
    run:
        import random
        random.seed(42)
        seeds = random.sample(range(10**4, 10**5 -1), params.n_seeds)
        with open(output[0], "w") as f:
            for seed in seeds:
                f.write(f"{seed}\n")


def umap_runtime(wildcards):
    runtime_dict = {
        "ath": "100h", # to be confirmed
        "ptr": "1.5h",
        "tma": "3h",
        "lch": "1.5h",
        "gar": "1.5h",
        "osa": "7h",
        "zma": "12h" # to be confirmed
    }
    return runtime_dict.get(wildcards.sample)

def umap_mem_mb_per_cpu(wildcards):
    mem_dict = {
        "ptr": 1000,
        "tma": 1000,
        "lch": 600,
        "gar": 900,
        "osa": 1500
    }
    return mem_dict.get(wildcards.sample)

rule run_umap_for_a_seed:
    container: "src/umap-learn/umap-learn_0.5.9.post2.sif"
    threads: 5
    resources:
        mem_mb_per_cpu=umap_mem_mb_per_cpu,
        runtime=umap_runtime
    input:
        rmd_output="outputs/Seurat/{sample}",
        param_table="outputs/UMAP/param_table.csv"
    output:
        dout=directory("outputs/UMAP/{sample}/seed{seed}")
    log:
        "logs/UMAP/seed{seed}/{sample}.log"
    benchmark:
        "benchmarks/UMAP/seed{seed}/{sample}.txt"
    shell:
        """
        python scripts/run_umap_for_a_seed.py \
            --rmd_input {input.rmd_output} \
            --param_table {input.param_table} \
            --seed {wildcards.seed} \
            --output_dir {output.dout} \
            --n_process {threads} \
        > {log} 2>&1
        """

rule call_umap_per_sample:
    input:
        param_table="outputs/UMAP/param_table.csv",
        seed_list="outputs/UMAP/seeds.txt",
        umaps=lambda wildcards: [
            f"outputs/UMAP/{wildcards.sample}/seed{seed}" for seed in seeds
        ]
    output:
        "outputs/UMAP/call_umap_per_sample_per_seed/{sample}/done.txt",
    log:
        "logs/UMAP/call_umap_per_sample_per_seed/{sample}.log"
    shell:
        """
        echo "UMAP calls completed for sample {wildcards.sample}" > {output[0]}
        """

# rule run_umap_test:
#     container: "src/umap-learn/umap-learn_0.5.9.post2"
#     threads: 95
#     resources:
#         mem_mb=600
#     output:
#         dout=directory("outputs/UMAP/umap_test/{repeat}")
#     log:
#         "logs/UMAP/umap_test/{repeat}.log"
#     shell:
#         """
#         python scripts/run_umap_repeat.py \
#             --n_jobs {threads} \
#         2> {log}
#         1> {log}
#         """
#
# rule run_umap_for_all_seeds:
#     container: "src/umap-learn/umap-learn_0.5.9.post2"
#     threads: 5
#     resources:
#         mem_mb=600
#     input:
#         rmd_output="outputs/Seurat/{sample}",
#         param_table="outputs/UMAP/param_table.csv",
#         seed_list="outputs/UMAP/seeds.txt"
#     output:
#         dout=directory("outputs/UMAP/umap/set{param_set}/{sample}")
#     log:
#         "logs/UMAP/umap/set{param_set}/{sample}.log"
#     shell:
#         """
#         python scripts/run_umap_for_all_params.py \
#             --rmd_input {input.rmd_output} \
#             --param_table {input.param_table} \
#             --param_set {wildcards.param_set} \
#             --seed_list {input.seed_list} \
#             --output_dir {output.dout} \
#             --n_process {threads} \
#         2> {log} 
#         1> {log}
#         """
# 
# rule call_umap_per_sample_per_param_set:
#     input:
#         param_table="outputs/UMAP/param_table.csv",
#         umaps=lambda wildcards: [
#             f"outputs/UMAP/umap/set{param_set}/{wildcards.sample}" for param_set in range(n_combinations)
#         ]
#     output:
#         "outputs/UMAP/call_umap_per_sample_per_param_set/{sample}/done.txt"
#     log:
#         "logs/UMAP/call_umap_per_sample_per_param_set/{sample}.log"
#     shell:
#         """
#         echo "UMAP calls completed for sample {wildcards.sample}" > {output[0]}
#         """
