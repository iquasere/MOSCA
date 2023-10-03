rule protein_report:
    input:
        expand("{output}/Annotation/{sample}/UPIMAPI_results.tsv", output=OUTPUT, sample=set(EXPS['Sample'])),
        expand("{output}/Annotation/{sample}/reCOGnizer_results.xlsx", output=OUTPUT, sample=set(EXPS["Sample"])),
        expand("{output}/Quantification/{sample}/mg.readcounts.norm", output=OUTPUT, sample=set(mg_exps['Sample'])),
        expand("{output}/Quantification/{sample}/mt.readcounts.norm", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}/spectracounts.tsv", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        f"{OUTPUT}/MOSCA_Protein_Report.xlsx",
        f"{OUTPUT}/Quantification/mg.readcounts",
        f"{OUTPUT}/Quantification/mt.readcounts" if len(mt_exps) > 0 else f"{OUTPUT}/Metaproteomics/mp.spectracounts"
    threads:
        1
    params:
        output = OUTPUT,
        exps = f"{OUTPUT}/exps.tsv",
        report = 'protein'
    conda:
        "../envs/reports.yaml"
    script:
        "../scripts/main_reports.py"
