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
        output = lambda wildcards: f'{OUTPUT}/Binning/{wildcards.sample}',
        markerset = config["markerset"],
        iterative = config['do_iterative_binning'],
        reads = lambda wildcards, input: input.reads[0] if len(input.reads) == 1 else ",".join(input.reads)
    conda:
        "../envs/binning.yaml"
    script:
        "../scripts/binning.py"
