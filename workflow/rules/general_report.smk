rule general_report:
    input:
        expand("{output}/Annotation/{sample}/UPIMAPI_results.tsv", output=OUTPUT, sample=set(EXPS['Sample'])),
        expand("{output}/Annotation/{sample}/reCOGnizer_results.xlsx", output=OUTPUT, sample=set(EXPS["Sample"])),
        expand("{output}/Quantification/{sample}_mg.readcounts", output=OUTPUT, sample=set(mg_exps['Sample'])) if len(mg_exps) > 0 else [],
        expand("{output}/Quantification/{sample}_mg_norm.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])) if (len(mg_exps) > 0 and config['do_assembly']) else [],
        expand("{output}/Quantification/{sample}_mt.readcounts", output=OUTPUT, sample=set(mt_exps['Sample'])) if len(mt_exps) > 0 else [],
        expand("{output}/Quantification/{sample}_mt_norm.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])) if (len(mt_exps) > 0 and config['do_assembly']) else [],
        expand("{output}/Metaproteomics/{sample}_mp.spectracounts", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        expand("{output}/MOSCA_{sample}_General_Report.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])),
        f"{OUTPUT}/MOSCA_General_Report.xlsx",
        f"{OUTPUT}/Quantification/dea_input.tsv",
        f"{OUTPUT}/Quantification/mg_entry_quant.tsv" if len(mg_exps) > 0 else [],
        f"{OUTPUT}/Quantification/mt_entry_quant.tsv" if len(mt_exps) > 0 else [],
        f"{OUTPUT}/Metaproteomics/mp_entry_quant.tsv" if len(mp_exps) > 0 else []
    threads:
        1
    params:
        output = OUTPUT,
        exps = f"{OUTPUT}/exps.tsv",
        did_assembly = config['do_assembly']
    conda:
        "../envs/reports.yaml"
    script:
        "../scripts/general_report.py"
