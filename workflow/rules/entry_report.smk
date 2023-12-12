rule entry_report:
    input:
        p_reports = expand("{output}/MOSCA_{sample}_Protein_Report.tsv", output=OUTPUT, sample=set(mg_exps['Sample'])),
        norm = f"{OUTPUT}/Quantification/mt_normalized.tsv" if len(mt_exps) > 0 else f"{OUTPUT}/Metaproteomics/mp_normalized.tsv"
    output:
        f"{OUTPUT}/MOSCA_Entry_Report.xlsx",
        f"{OUTPUT}/MOSCA_Entry_Report.tsv"
    threads:
        config["threads"]
    conda:
        "../envs/reports.yaml"
    params:
        output = OUTPUT,
        exps = f"{OUTPUT}/exps.tsv",
    script:
        "../scripts/entry_report.py"
