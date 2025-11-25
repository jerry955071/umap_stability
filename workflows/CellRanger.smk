configfile: "config/config.json"

cellranger=config["executables"]["cellranger"]

rule cellranger_mkgtf:
    threads: 1
    input:
        annotation=lambda wildcards: query(config["references"], "species", wildcards.not_lch)["annotation"],
    output:
        gtf=temp("outputs/CellRanger/mkgtf/{not_lch}.gtf"),
        gtf_filtered="outputs/CellRanger/mkgtf/{not_lch}.gtf-filtered"
    log:
        "logs/CellRanger/mkgtf/{not_lch}.log"
    shell:
        """
        # convert to gtf if gff3
        if [[  "{input.annotation}" == *gff*  ]]; then
            {docker_run} quay.io/biocontainers/gffread:0.12.1--h8b12597_0 \
                gffread {input.annotation} -T -o {output.gtf} \
                    2> {log} \
                    1> {log}
        else
            cp {input.annotation} {output.gtf}
        fi
        
        # cellranger
        {cellranger} mkgtf \
            {output.gtf} {output.gtf_filtered} \
            --attribute=gene_biotype:protein_coding \
        2>> {log} \
        1>> {log}             
        """

rule cellranger_mkgtf_lch:
    threads: 1
    input:
        annotation=lambda _: query(config["references"], "species", "lch")["annotation"],
    output:
        gtf=temp("outputs/CellRanger/mkgtf/lch.gtf"),
        gtf_filtered="outputs/CellRanger/mkgtf/lch.gtf-filtered"
    log:
        "logs/CellRanger/mkgtf/lch.log"
    shell:
        """
        # convert gff to gtf
        {docker_run} quay.io/biocontainers/gffread:0.12.1--h8b12597_0 \
            gffread {input.annotation} -T -o {output.gtf} \
                2> {log} \
                1> {log}
    
        # add gene_id to attributes for Cell Ranger compatibility          
        python3 scripts/repair_lch_gtf.py {output.gtf} {output.gtf}_temp
        mv {output.gtf}_temp {output.gtf}

        # cellranger
        {cellranger} mkgtf \
            {output.gtf} {output.gtf_filtered} \
            --attribute=gene_biotype:protein_coding \
        2>> {log} \
        1>> {log}
        """

rule cellranger_mkref:
    threads: 8
    resources:
        mem_gb=20
    input:
        assembly=lambda wildcards: query(config["references"], "species", wildcards.species)["assembly"],
        gtf="outputs/CellRanger/mkgtf/{species}.gtf-filtered"
    output:
        directory("outputs/CellRanger/mkref/{species}")
    log:
        "logs/CellRanger/mkref/{species}.log"
    shell:
        """
        {cellranger} mkref \
            --output-dir={output} \
            --genome=$(basename {output}) \
            --fasta={input.assembly} \
            --genes={input.gtf} \
            --nthreads={threads} \
            --memgb={resources.mem_gb} \
        2> {log} \
        1> {log}
        """

# rule cellranger_count:
#     threads: 24
#     resources:
#         mem_gb=48
#     input:
#         transcriptome=lambda wildcards: "outputs/CellRanger/mkref/%s" % query(config["samples"], "name", wildcards.sample)["species"],
#         fastq_path=lambda wildcards: query(config["samples"], "name", wildcards.sample)["fastq_path"]
#     output:
#         outdir=directory("outputs/CellRanger/count/{sample}")
#     params:
#         create_bam="false",
#         fastq_path=lambda wildcards: ",".join(query(config["samples"], "name", wildcards.sample)["fastq_path"])
#     log: 
#         "logs/CellRanger/count/{sample}.log"
#     shell:
#         """
#         {cellranger} count \
#             --id={wildcards.sample} \
#             --output-dir={output.outdir} \
#             --transcriptome={input.transcriptome} \
#             --fastqs={params.fastq_path} \
#             --create-bam={params.create_bam} \
#             --localcores={threads} \
#             --localmem={resources.mem_gb} \
#             --nosecondary \
#             --sample=sample \
#         2> {log} \
#         1> {log}
#         """

# function to get reference path based on run
def _ref_for_run(wildcards):
    for i in config["samples"]:
        if wildcards.run in i["runs"]:
            species = i["species"]
            return "outputs/CellRanger/mkref/%s" % species

rule cellranger_count:
    threads: 24
    resources:
        mem_gb=48
    input:
        transcriptome=_ref_for_run,
        fastq_path="rawdata/fastq_paths/{run}"
    output:
        directory("outputs/CellRanger/count/{run}")
    params:
        create_bam="false"
    log: 
        "logs/CellRanger/count/{run}.log"
    shell:
        """
        {cellranger} count \
            --id={wildcards.run} \
            --output-dir={output} \
            --transcriptome={input.transcriptome} \
            --fastqs={input.fastq_path} \
            --create-bam={params.create_bam} \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
            --nosecondary \
            --sample=sample \
        2> {log} \
        1> {log}
        """