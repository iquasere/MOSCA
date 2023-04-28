rule entry_report:
    input:
        f"{OUTPUT}/MOSCA_Protein_Report.xlsx",
        expand("{output}/Quantification/{sample}/mt_normalized.tsv", output=OUTPUT, sample=set(mt_exps['Sample'])),
        expand("{output}/Metaproteomics/{sample}/mp_normalized.tsv", output=OUTPUT, sample=set(mp_exps['Sample']))
    output:
        f"{OUTPUT}/MOSCA_Entry_Report.xlsx"
    threads:
        config["threads"]
    conda:
        "../envs/reports.yaml"
    params:
        output = OUTPUT,
        exps = f"{OUTPUT}/exps.tsv",
        report = 'entry'
    script:
        "../scripts/main_reports.py"
