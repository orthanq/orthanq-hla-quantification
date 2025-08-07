# wrappers should be used once they are ready
rule generate_candidates:
    input:
        hla_genes="results/preparation/hla_gen.fasta",
        xml="results/preparation/hla.xml",
        genome=genome,
    output:
        vcfs=expand("results/orthanq/candidate_variants/{hla}.vcf", hla=loci),
        bam="results/orthanq/candidate_variants/alleles_alignment_sorted.bam"
    log:
        "logs/candidates/candidates.log",
    conda:
        "../envs/orthanq.yaml"
    params:
        output_folder=lambda wc, output: os.path.dirname(output.vcfs[0]),
    threads: config["threads"]
    benchmark:
        "benchmarks/orthanq_candidates/orthanq_candidates.tsv"
    shell:
        "orthanq candidates hla --alleles {input.hla_genes} --genome {input.genome} --xml {input.xml} --threads {threads} --output {params.output_folder} 2> {log}"


rule create_graph_symlink:
    input:
        "results/preparation/hprc-v1.0-mc-grch38.xg",
    output:
        "results/preparation/linked_graphs/{sample}_{hla}/hprc-v1.0-mc-grch38.xg",
    params:
        target_dir=subpath(output[0], parent=True),
    log:
        "logs/link_graph/{sample}_{hla}.log",
    shell:
        """
        mkdir -p {params.target_dir}
        ln -sfr {input} {output}
        """


rule preprocess:
    input:
        bwa_index=rules.bwa_index.output,
        candidate_variants="results/orthanq/candidate_variants/{hla}.vcf",
        genome=genome,
        genome_fai=genome_fai,
        pangenome="results/preparation/linked_graphs/{sample}_{hla}/hprc-v1.0-mc-grch38.xg",
        orthanq_input=get_orthanq_input,
    output:
        bcf="results/orthanq/preprocess/{sample}_{hla}/{sample}_{hla}.bcf",
        aux=temp(multiext("results/preparation/linked_graphs/{sample}_{hla}/hprc-v1.0-mc-grch38", ".dist", ".giraffe.gbz", ".shortread.withzip.min", ".shortread.zipcodes")),
    log:
        "logs/preprocess/{sample}_{hla}.log",
    conda:
        "../envs/orthanq.yaml"
    params:
        bwa_idx_prefix=lambda wc, input: os.path.splitext(input.bwa_index[0])[0],
        input_params=get_orthanq_input_params,
    threads: config["threads"]
    benchmark:
        "benchmarks/orthanq_preprocess/{sample}_{hla}.tsv"
    shell:
        "orthanq preprocess hla --genome {input.genome} --bwa-index {params.bwa_idx_prefix} --haplotype-variants {input.candidate_variants} --vg-index {input.pangenome} --output {output.bcf} {params.input_params} --threads {threads} 2> {log}"


# #wrappers should be used once they are ready
rule quantify:
    input:
        haplotype_variants="results/orthanq/candidate_variants/{hla}.vcf",
        haplotype_calls="results/orthanq/preprocess/{sample}_{hla}/{sample}_{hla}.bcf",
        xml="results/preparation/hla.xml",
    output:
        table="results/orthanq/calls/{sample}_{hla}/{sample}_{hla}.csv",
        three_field_solutions="results/orthanq/calls/{sample}_{hla}/3_field_solutions.json",
        two_field_solutions="results/orthanq/calls/{sample}_{hla}/2_field_solutions.json",
        final_solutions="results/orthanq/calls/{sample}_{hla}/final_solution.json",
        lp_solution="results/orthanq/calls/{sample}_{hla}/lp_solution.json",
        two_field_table="results/orthanq/calls/{sample}_{hla}/2-field.csv",
        g_groups="results/orthanq/calls/{sample}_{hla}/G_groups.csv",
    log:
        "logs/orthanq_call/{sample}_{hla}.log",
    conda:
        "../envs/orthanq.yaml"
    params:
        prior=config["orthanq_prior"],
    resources:
        mem_mb=5000,
    benchmark:
        "benchmarks/orthanq_quantify/{sample}_{hla}.tsv"
    shell:
        "orthanq call hla --haplotype-variants {input.haplotype_variants} --xml {input.xml} "
        " --haplotype-calls {input.haplotype_calls} --prior {params.prior} --output {output.table} 2> {log}"
