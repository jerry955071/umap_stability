rule procrustes_alignment_per_n_embeddings:
    container: "src/umap-learn/umap-learn_0.5.9.post2.sif"
    threads: 80
    resources:
        mem_mb_per_cpu=5000,
        runtime=240
    input:
        embeddings=expand(
            "outputs/UMAP/{sample}/seed{seed}/param_set{param_set}.csv",
            sample=lambda wildcards: wildcards.sample,
            param_set=lambda wildcards: wildcards.param_set,
            seed=seeds
        )
    output:
        directory("outputs/AlignEmbeddings/procrustes_alignment/{sample}/param_set{param_set}/per_{n}_embeddings")
    params:
        n=lambda wildcards: wildcards.n
    log:
        "logs/AlignEmbeddings/procrustes_alignment/{sample}/param_set{param_set}/per_{n}_embeddings.log"
    benchmark:
        "benchmarks/AlignEmbeddings/procrustes_alignment/{sample}/param_set{param_set}/per_{n}_embeddings.log"
    shell:
        """
        exec python scripts/procrustes_alignment_per_n_embeddings.py \
            --n_process {threads} \
            --n_embeddings {params.n} \
            --output_dir {output} \
            {input.embeddings} \
        > {log} 2>&1
        """
