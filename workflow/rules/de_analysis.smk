include: "common.smk"

rule de_analysis:
    input:
        expand("{output}/Quantification/{sample}/mt.readcounts", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}/mp_normalized.tsv", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        expand("{output}/Quantification/{sample}/condition_treated_results.tsv", output=OUTPUT,
            sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}/condition_treated_results.tsv", output=OUTPUT,
            sample=set(mp_exps['Sample']))
    threads:
        15
    params:
        conditions = (lambda wildcards, input: ",".join(not_mg_exps[not_mg_exps['Sample'] == os.path.basename(
            os.path.dirname(input[0]))]['Condition'].tolist())),
        minimum_fold_change=config["minimum_differential_expression"],
        fdr=config["false_discovery_rate"],
        out_dir = (lambda wildcards, input: os.path.dirname(input[0])),
        type_of_data = (lambda wildcards, input: "rna_seq" if "Quantification" in input[0] else "proteomics")
    conda:
        "../envs/de_analysis.yaml"
    shell:
        "Rscript {SCRIPTS_DIR}/de_analysis.R --counts {input} --conditions {params.conditions} "
        "--output {params.out_dir} --foldchange {params.minimum_fold_change} "
        "--fdr {params.fdr} --datatype {params.type_of_data}"
