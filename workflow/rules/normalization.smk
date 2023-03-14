rule normalization:
    input:
        expand("{output}/Quantification/{sample}/mt.readcounts", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}/mp.spectracounts", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        expand("{output}/Quantification/{sample}/mt_normalized.tsv",output=OUTPUT,sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}/mp_normalized.tsv",output=OUTPUT,sample=set(mp_exps['Sample']))
    threads:
        1
    params:
        nm = config["normalization_method"]
    conda:
        "../envs/normalization.yaml"
    shell:
        "Rscript ../scripts/normalization.R -c {input} -m {params.nm} -o {output}"
