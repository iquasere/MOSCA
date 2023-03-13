include: "common.smk"

rule fastq2fasta:
    input:
        fastq2fasta_input
    output:
        f"{OUTPUT}/Preprocess/piled_{{sample}}.fasta"
    threads:
        1
    shell:
        "cat {input} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr '\\t' '\\n' > {output}"

rule gene_calling:
    input:
        gene_calling_input
    output:
        expand("{output}/Annotation/{{sample}}/fgs.faa", output=OUTPUT),
        expand("{output}/Annotation/{{sample}}/fgs.ffn", output=OUTPUT)
    threads:
        config["threads"]
    params:
        error_model = "complete" if config['do_assembly'] else config["error_model"],
        complete = '1' if config['do_assembly'] else '0'
    conda:
        "../envs/gene_calling.yaml"
    shell:
        'run_FragGeneScan.pl -thread={threads} -genome={input} -out={OUTPUT}/Annotation/{wildcards.sample}/fgs '
        '-complete={params.complete} -train=./{params.error_model}'