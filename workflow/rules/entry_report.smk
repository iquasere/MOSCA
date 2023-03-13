include: "common.smk"

rule entry_report:
    input:
        f"{OUTPUT}/MOSCA_Protein_Report.xlsx",
        expand("{output}/Quantification/{sample}/mt_normalized.tsv", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}/mp_normalized.tsv", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        f"{OUTPUT}/MOSCA_Entry_Report.xlsx"
    threads:
        1
    conda:
        "../envs/reports.yaml"
    shell:
        "python {SCRIPTS_DIR}/main_reports.py -o {OUTPUT} -e {OUTPUT}/exps.tsv --entry-report"
