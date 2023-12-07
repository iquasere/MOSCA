rule summary_report:
    input:
        expand("{output}/MOSCA_{sample}_Protein_Report.tsv", output=OUTPUT, sample=set(EXPS['Sample'])),
        f"{OUTPUT}/MOSCA_Entry_Report.xlsx",
        f"{OUTPUT}/DE_analysis/condition_treated_results.tsv"
    output:
        f"{OUTPUT}/technical_report.tsv",
        f"{OUTPUT}/MOSCA_General_Report.tsv",
        f"{OUTPUT}/MOSCA_results.zip"
    threads:
        1
    params:
        output=OUTPUT,
        cutoff=config["significance_threshold"]
    conda:
        "../envs/summary.yaml"
    script:
        "../scripts/summary_report.py"
