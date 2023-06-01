rule upimapi:
    input:
        expand("{output}/Annotation/{{sample}}/fgs.faa", output=OUTPUT)
    output:
        expand("{output}/Annotation/{{sample}}/UPIMAPI_results.tsv", output=OUTPUT)
    threads:
        config["threads"]
    params:
        rd = config["resources_directory"],
        upimapi_database = config["upimapi_database"],
        taxids = f' --taxids {config["upimapi_taxids"]}' if config["upimapi_database"] == 'taxids' else '',
        max_target_seqs = config["upimapi_max_target_seqs"],
        cols = '&'.join(config['uniprot_columns']),
        check_db = '' if config['upimapi_check_db'] else ' --skip-db-check',
        diamond_mode = config["upimapi_search_mode"],
        memory = config["max_memory"]
    conda:
        "../envs/upimapi.yaml"
    shell:
        'upimapi -i {input} -t {threads} -o {OUTPUT}/Annotation/{wildcards.sample}/{wildcards.part} '
        '-rd {params.rd} -db {params.upimapi_database} -mts {params.max_target_seqs}{params.taxids} '
        '-cols "{params.cols}" --max-memory {params.memory} --diamond-mode {params.diamond_mode}{params.check_db}'
