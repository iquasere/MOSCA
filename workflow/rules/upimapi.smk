include: "common.smk"

rule split_gene_calling:
    input:
        expand("{output}/Annotation/{{sample}}/fgs.faa", output=OUTPUT)
    output:
        expand("{output}/Annotation/{{sample}}/fgs.faa.split/fgs.part_{{part}}.faa", output=OUTPUT)
    threads:
        config["threads"]
    params:
        parts = config["split_gene_calling"]
    conda:
        "../envs/seqkit.yaml"
    shell:
        'seqkit split --by-part {params.parts} -j {threads} {input}'

rule upimapi:
    input:
        expand("{output}/Annotation/{{sample}}/fgs.faa.split/fgs.part_{{part}}.faa", output=OUTPUT)
    output:
        expand("{output}/Annotation/{{sample}}/{{part}}/UPIMAPI_results.tsv", output=OUTPUT)
    threads:
        config["threads"]
    params:
        rd = config["resources_directory"],
        upimapi_database = config["upimapi_database"],
        taxids = f' --taxids {config["upimapi_taxids"]}' if config["upimapi_database"] == 'taxids' else '',
        max_target_seqs = config["upimapi_max_target_seqs"],
        cols = '&'.join(config['uniprot_columns']),
        check_db = '' if config['upimapi_check_db'] else ' --skip-db-check'
    conda:
        "../envs/upimapi.yaml"
    shell:
        'upimapi.py -i {input} -t {threads} -o {OUTPUT}/Annotation/{wildcards.sample}/{wildcards.part} '
        '-rd {params.rd} -db {params.upimapi_database} -mts {params.max_target_seqs}{params.taxids} '
        '-cols "{params.cols}"{params.check_db}'

rule join_upimapi:
    input:
        expand("{output}/Annotation/{sample}/{part}/UPIMAPI_results.tsv", output=OUTPUT,
            part=[f'{"0" * (3 - len(str(i + 1)))}{i + 1}' for i in range(config["split_gene_calling"])],
            sample=set(EXPS["Sample"]))
    output:
        expand("{output}/Annotation/{sample}/UPIMAPI_results.tsv", output=OUTPUT,
            sample=set(EXPS["Sample"]))
    shell:
        "awk 'FNR==1{{if(NR!=1)next;}}{{print}}' {input} > {output}"
