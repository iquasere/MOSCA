rule normalization:
    input:
        f"{OUTPUT}/Quantification/mg_readcounts.tsv",
        f"{OUTPUT}/Quantification/mt_readcounts.tsv" if len(mt_exps) > 0 else f"{OUTPUT}/Metaproteomics/mp_spectracounts.tsv"
    output:
        f"{OUTPUT}/Quantification/mg_normalized.tsv",
        f"{OUTPUT}/Quantification/mt_normalized.tsv" if len(mt_exps) > 0 else f"{OUTPUT}/Metaproteomics/mp_normalized.tsv"
    threads:
        1
    params:
        norm_method = config["normalization_method"],
        imput_method = config["imputation_method"]
    conda:
        "../envs/normalization.yaml"
    script:
        "../scripts/normalization.R"
