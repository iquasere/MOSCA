include: "common.smk"

rule metaproteomics:
    input:
        [directory(folder) for folder in mp_exps[mp_exps['Sample'] == (lambda wildcards: wildcards.sample)]['Files']],
        "{output}/Annotation/{sample}/UPIMAPI_results.tsv"
    output:
        "{output}/Metaproteomics/{sample}/spectracounts.tsv"
    threads:
        config["threads"]
    params:
        folders=lambda wildcards: ','.join(mp_exps[mp_exps['Sample'] == wildcards.sample]['Files'].tolist()),
        names=lambda wildcards: ','.join(mp_exps[mp_exps['Sample'] == wildcards.sample]['Name'].tolist()),
        contaminants_database=config["proteomics_contaminants_database"],
        protease=config["protease"] if config["protease_file"] == "" else config["protease_file"],
        max_memory=config["max_memory"],
        resources_directory=config["resources_directory"]
    conda:
        "../envs/metaproteomics.yaml"
    shell:
        "python {SCRIPTS_DIR}/metaproteomics_analyser.py -sf {params.folders} -ns {params.names} -t {threads} "
        "-o {OUTPUT}/Metaproteomics/{wildcards.sample} -db {OUTPUT}/Annotation/{wildcards.sample}/fgs.faa "
        "-cdb {params.contaminants_database} --protease {params.protease} "
        "-ur {OUTPUT}/Annotation/{wildcards.sample}/UPIMAPI_results.tsv -mmem {params.max_memory} "
        "-rd {params.resources_directory}"
