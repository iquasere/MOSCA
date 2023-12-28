rule protein_report:
    input:
        expand("{output}/Annotation/{sample}/UPIMAPI_results.tsv", output=OUTPUT, sample=set(EXPS['Sample'])),
        expand("{output}/Annotation/{sample}/reCOGnizer_results.xlsx", output=OUTPUT, sample=set(EXPS["Sample"])),
        expand("{output}/Quantification/{sample}_mg.readcounts", output=OUTPUT, sample=set(mg_exps['Sample'])),
        expand("{output}/Quantification/{sample}_mg_norm.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])),
        expand("{output}/Quantification/{sample}_mt.readcounts", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Quantification/{sample}_mt_norm.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}_mp.spectracounts", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        expand("{output}/MOSCA_{sample}_General_Report.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])),
        f"{OUTPUT}/MOSCA_General_Report.xlsx",
        f"{OUTPUT}/Quantification/dea_input.tsv",
        f"{OUTPUT}/Quantification/mg_entry_quant.tsv",
        f"{OUTPUT}/Quantification/mt_entry_quant.tsv" if len(mt_exps) > 0 else f"{OUTPUT}/Metaproteomics/mp_entry_quant.tsv"
    threads:
        1
    params:
        output = OUTPUT,
        exps = f"{OUTPUT}/exps.tsv",
    conda:
        "../envs/reports.yaml"
    script:
        "../scripts/general_report.py"
