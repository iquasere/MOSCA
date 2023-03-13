include: "common.smk"

rule binning:
    input:
        reads = expand("{output}/Preprocess/{{sample}}{fr}.fastq", output=OUTPUT,
            fr=(['_forward', '_reverse'] if EXPS["Files"].str.contains(',').tolist() else '')),
        contigs = expand("{output}/Assembly/{{sample}}/scaffolds.fasta", output=OUTPUT)
    output:
        expand("{output}/Binning/{{sample}}/checkm.tsv", output=OUTPUT, sample=set(EXPS['Sample']))
    threads:
        config["threads"]
    params:
        markerset = config["markerset"],
        iterative_binning = ' --iterative' if config['do_iterative_binning'] else '',
        reads = lambda wildcards, input: input.reads[0] if len(input.reads) == 1 else ",".join(input.reads)
    conda:
        "../envs/binning.yaml"
    shell:
        "python {SCRIPTS_DIR}/binning.py -c {input.contigs} -t {threads} -o {OUTPUT}/Binning/{wildcards.sample} "
        "-r {params.reads} -mset {params.markerset}{params.iterative_binning}"
