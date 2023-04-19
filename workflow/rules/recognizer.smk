rule recognizer:
    input:
        orfs = expand("{output}/Annotation/{{sample}}/fgs.faa", output=OUTPUT),
        upimapi_results = expand("{output}/Annotation/{{sample}}/UPIMAPI_results.tsv", output=OUTPUT)
    output:
        expand("{output}/Annotation/{{sample}}/reCOGnizer_results.xlsx", output=OUTPUT)
    threads:
        config["threads"]
    params:
        resources_directory = config["resources_directory"],
        recognizer_databases = ','.join(config["recognizer_databases"]),
        download_cdd_resources = '' if not config['download_cdd_resources'] else ' -dr'
    conda:
        "../envs/recognizer.yaml"
    shell:
        "recognizer -f {input.orfs} -t {threads} -o {OUTPUT}/Annotation/{wildcards.sample} "
        "-rd {params.resources_directory} -dbs {params.recognizer_databases} -sd "
        #"--tax-file {input.upimapi_results} --protein-id-col Entry --tax-col 'Organism (ID)' --species-taxids 
        "--quiet{params.download_cdd_resources}"
