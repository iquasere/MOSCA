include: "common.smk"

rule summary_report:
    input:
        f"{OUTPUT}/MOSCA_Protein_Report.xlsx",
        f"{OUTPUT}/MOSCA_Entry_Report.xlsx",
        expand("{output}/Quantification/{sample}/condition_treated_results.tsv", output=OUTPUT,
            sample=set(mt_exps['Sample'].tolist())),
        expand("{output}/Metaproteomics/{sample}/condition_treated_results.tsv", output=OUTPUT,
            sample=set(mp_exps['Sample']))
    output:
        f"{OUTPUT}/technical_report.tsv",
        f"{OUTPUT}/MOSCA_General_Report.tsv",
        f"{OUTPUT}/MOSCA_results.zip"
    threads:
        1
    conda:
        "../envs/summary.yaml"
    shell:
        "python {SCRIPTS_DIR}/summary_report.py -o {OUTPUT}"
