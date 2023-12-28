rule fastq2fasta:
    input:
        expand("{output}/Preprocess/Trimmomatic/quality_trimmed_{name}{fr}.fq", output=OUTPUT,
            fr=(['_forward_paired', '_reverse_paired'] if EXPS["Files"].str.contains(',').tolist() else ''),
            name=lambda wildcards: wildcards.sample)
    output:
        f"{OUTPUT}/Preprocess/piled_{{sample}}.fasta"
    threads:
        1
    shell:
        "cat {input} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr '\\t' '\\n' > {output}"

rule gene_calling:
    input:
        (f"{OUTPUT}/Assembly/{{sample}}/scaffolds.fasta" if config['do_assembly'] else
        f"{OUTPUT}/Preprocess/piled_{{sample}}.fasta")
    output:
        f"{OUTPUT}/Annotation/{{sample}}/fgs.faa",
        f"{OUTPUT}/Annotation/{{sample}}/fgs.ffn"
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