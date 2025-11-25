def _seurat_input(wildcards):
    # Get input data for the Seurat pipeline
    # For samples without processed data, "tma" and "osa", return CellRanger count output directories
    # For other samples, return the notebook file (always present), the processed data will be downloaded as shown in the notebook
    if wildcards.sample in ["tma", "osa"]:
        return [f"outputs/CellRanger/count/{run}" for run in query(config["samples"], "name", wildcards.sample)["runs"]]
    else:
        return f"notebooks/run-seurat/{wildcards.sample}.Rmd"

rule seurat_preprocessing:
    threads:
        8
    input:
        _seurat_input
    output:
        directory("outputs/Seurat/{sample}")
    params:
        rmd=lambda wildcards: f"notebooks/run-seurat/{wildcards.sample}.Rmd"
    log:
        "logs/Seurat/preprocessing/{sample}.log"
    shell:
        """
        {docker_run} --cpus={threads} chiaenu/seurat-rmd:5.0.0 \
            Rscript -e "rmarkdown::render('{params.rmd}', output_dir='outputs/Seurat/{wildcards.sample}', output_file='preprocessing.html', knit_root_dir='/data')" \
        2> {log} \
        1> {log}
        """
