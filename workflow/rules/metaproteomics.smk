rule metaproteomics:
    input:
        [directory(folder) for folder in mp_exps[mp_exps['Sample'] == (lambda wildcards: wildcards.sample)]['Files']],
        "{output}/Annotation/{sample}/UPIMAPI_results.tsv"
    output:
        "{output}/Metaproteomics/{sample}_mp_spectracounts.tsv"
    threads:
        config["threads"]
    params:
        output = lambda wildcards: f'{OUTPUT}/Metaproteomics/{wildcards.sample}',
        mg_db = lambda wildcards: f'{OUTPUT}/Annotation/{wildcards.sample}/fgs.faa',
        up_res = lambda wildcards: f'{OUTPUT}/Annotation/{wildcards.sample}/UPIMAPI_results.tsv',
        spectra_folders = lambda wildcards: mp_exps[mp_exps['Sample'] == wildcards.sample]['Files'].tolist(),
        names = lambda wildcards: mp_exps[mp_exps['Sample'] == wildcards.sample]['Name'].tolist(),
        contaminants_database = config["proteomics_contaminants_database"],
        protease = config["protease"] if config["protease_file"] == "" else config["protease_file"],
        max_memory = config["max_memory"],
        resources = config["resources_directory"],
        add_reference_proteomes = config["metaproteomics_add_reference_proteomes"],
        inside_container = config["inside_container"]
    conda:
        "../envs/metaproteomics.yaml"
    script:
        "../scripts/metaproteomics.py"
