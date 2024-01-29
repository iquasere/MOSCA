rule merge_reads:
    input:
        forward_reads=f"{OUTPUT}/Preprocess/Trimmomatic/quality_trimmed_{{name}}_forward_paired.fq",
        reverse_reads=f"{OUTPUT}/Preprocess/Trimmomatic/quality_trimmed_{{name}}_reverse_paired.fq"
    output:
        expand("{output}/Preprocess/Trimmomatic/{{name}}.assembled.fastq", output=OUTPUT)
    threads:
        config["threads"]
    params:
        memory = config["max_memory"]
    conda:
        "../envs/merge_reads.yaml"
    shell:
        "pear -f {input.forward_reads} -r {input.reverse_reads} -o {OUTPUT}/Preprocess/Trimmomatic/{wildcards.name} "
        "-j {threads} -y {params.memory}G"

rule fastq2fasta:
    input:
        gene_calling_input
    output:
        expand("{output}/Preprocess/piled_{{sample}}.fasta", output=OUTPUT)
    threads:
        1
    shell:
        "cat {input} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr '\\t' '\\n' > {output}"

rule gene_calling:
    input:
        (expand("{output}/Assembly/{{sample}}/scaffolds.fasta", output=OUTPUT) if config['do_assembly'] else
        expand("{output}/Preprocess/piled_{{sample}}.fasta", output=OUTPUT))
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