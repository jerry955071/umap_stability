rule inter_cluster_angle_per_sample:
    container: "src/umap-learn/umap-learn_0.5.9.post2.sif"
    threads: 80
    resources:
        mem_mb_per_cpu=5000,
        runtime=240
    input:
        flag="outputs/UMAP/call_umap_per_sample_per_seed/{sample}/done.txt",
        cluster_table="outputs/Seurat/{sample}/clusters.csv",
        param_table="outputs/UMAP/param_table.csv",
        seed_list="outputs/UMAP/seeds.txt"
    output:
        directory("outputs/InterClusterAngle/{sample}")
    params:
        umap_dir="outputs/UMAP/{sample}"
    log:
        "logs/InterClusterAngle/inter_cluster_angle_per_sample/{sample}.log"
    benchmark:
        "benchmarks/InterClusterAngle/inter_cluster_angle_per_sample/{sample}.txt"
    shell:
        """
        exec python scripts/inter_cluster_angle_per_sample.py \
            --umap_indir {params.umap_dir} \
            --cluster_table {input.cluster_table} \
            --param_table {input.param_table} \
            --seed_list {input.seed_list} \
            --output_dir {output} \
            --n_process {threads} \
        > {log} 2>&1
        """
