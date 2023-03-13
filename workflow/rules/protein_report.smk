include: "common.smk"

rule protein_report:
    input:
        expand("{output}/Annotation/{sample}/UPIMAPI_results.tsv", output=OUTPUT, sample=set(EXPS['Sample'])),
        expand("{output}/Annotation/{sample}/reCOGnizer_results.xlsx", output=OUTPUT, sample=set(EXPS["Sample"])),
        expand("{output}/Quantification/{name}.readcounts", output=OUTPUT, name=set(mt_exps['Name'])),
        expand("{output}/Metaproteomics/{sample}/spectracounts.tsv", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        f"{OUTPUT}/MOSCA_Protein_Report.xlsx",
        expand("{output}/Quantification/{sample}/mt.readcounts", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}/mp.spectracounts", output=OUTPUT, sample=set(mp_exps['Sample']))
    threads:
        1
    conda:
        "../envs/reports.yaml"
    shell:
        "python {SCRIPTS_DIR}/main_reports.py -o {OUTPUT} -e {OUTPUT}/exps.tsv --protein-report"
