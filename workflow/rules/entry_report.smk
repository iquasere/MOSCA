import os

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

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
    shell:
        "python ../scripts/main_reports.py -o {OUTPUT} -e {OUTPUT}/exps.tsv --entry-report"
