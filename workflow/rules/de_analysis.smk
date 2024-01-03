rule de_analysis:
    input:
        f"{OUTPUT}/Quantification/dea_input.tsv" if len(mt_exps) > 0 else f"{OUTPUT}/Metaproteomics/mp_normalized.tsv"
    output:
        f"{OUTPUT}/DE_analysis/condition_treated_results.tsv"
    threads:
        1
    params:
        conditions = (lambda wildcards, input: ",".join(not_mg_exps['Condition'].tolist())),
        foldchange = config["minimum_differential_expression"],
        fdr = config["significance_threshold"],
        output = f"{OUTPUT}/DE_analysis",
        datatype = (lambda wildcards, input: "rna_seq" if len(set(mt_exps['Sample'])) > 0 else "proteomics")
    conda:
        "../envs/de_analysis.yaml"
    script:
        "../scripts/de_analysis.R"
