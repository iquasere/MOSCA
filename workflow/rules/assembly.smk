rule join_reads:
    input:
        join_reads_input
    output:
        expand("{output}/Preprocess/{{sample}}{fr}.fastq", output=OUTPUT,
            fr=(['_forward', '_reverse'] if EXPS["Files"].str.contains(',').tolist() else ''))
    threads:
        1
    run:
        for file in input:
            if 'forward' in file:
                shell("touch {output}/Preprocess/{wildcards.sample}_forward.fastq; cat {file} >> "
                      "{output}/Preprocess/{wildcards.sample}_forward.fastq", output=OUTPUT)
            elif 'reverse' in file:
                shell("touch {output}/Preprocess/{wildcards.sample}_reverse.fastq; cat {file} >> "
                      "{output}/Preprocess/{wildcards.sample}_reverse.fastq", output=OUTPUT)
            else:
                shell("touch {output}/Preprocess/{wildcards.sample}.fastq; cat {file} >> "
                      "{output}/Preprocess/{wildcards.sample}.fastq", output=OUTPUT)

rule assembly:
    input:
        expand("{output}/Preprocess/{{sample}}{fr}.fastq", output=OUTPUT,
            fr=(['_forward', '_reverse'] if EXPS["Files"].str.contains(',').tolist() else ''))
    output:
        expand("{output}/Assembly/{{sample}}/contigs.fasta", output=OUTPUT,
            sample=set(EXPS['Sample'])),
        expand("{output}/Assembly/{{sample}}/scaffolds.fasta", output=OUTPUT,
            sample=set(EXPS['Sample']))
    threads:
        config["threads"]
    params:
        output = lambda wildcards: f'{OUTPUT}/Assembly/{wildcards.sample}',
        assembler = config["assembler"],
        reads = expand("{output}/Preprocess/{{sample}}{fr}.fastq", output=OUTPUT,
            fr=(['_forward', '_reverse'] if EXPS["Files"].str.contains(',').tolist() else '')),
        max_memory = config["max_memory"]
    conda:
        "../envs/assembly.yaml"
    script:
        "../scripts/assembly.py"
